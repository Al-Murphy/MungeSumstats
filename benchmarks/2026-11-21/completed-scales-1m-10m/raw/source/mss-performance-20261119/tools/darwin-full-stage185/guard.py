"""Darwin process-group and scratch-space supervision for one local run.

Callers supply trusted foreground workloads and private directories. A child
that deliberately changes session/group or writes outside its admitted paths
is outside this watchdog's containment model.
"""

from contextlib import contextmanager, ExitStack
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from itertools import islice
import json
import math
import os
from pathlib import Path
import re
import selectors
import signal
import stat
import subprocess
import sys
import threading
import time


_MAX_WALL_SECONDS = 16 * 60 * 60
_MAX_RSS_BYTES = 24_000_000_000
_DIR_FLAGS = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC
_FILE_FLAGS = os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW | os.O_CLOEXEC
_PS_ENV = {"PATH": "/usr/bin:/bin", "LC_ALL": "C", "TZ": "UTC"}
_PS_OUTPUT_BYTES = 16 * 1024
_OWNED_CLEANUP_SECONDS = 0.10
_PLAN_BYTES = 16 * 1024
_RUNNING_BYTES = 4 * 1024
_TERMINAL_BYTES = 32 * 1024
_PENDING_BYTES = 4 * 1024
_ARG_BYTES = 8 * 1024
_MAX_ARGUMENTS = 256
_MAX_LABEL_BYTES = 256
_STOP_SIGNALS = (signal.SIGINT, signal.SIGTERM, signal.SIGHUP)
_STATE = re.compile(r"[IRSTUZ][+<>AELNSsVWX-]*", re.ASCII)
_DECIMAL = re.compile(r"[0-9]+", re.ASCII)
_ENV_KEYS = frozenset("""
PATH HOME TMPDIR TMP TEMP LANG LANGUAGE LC_ALL LC_CTYPE LC_COLLATE LC_MESSAGES
LC_MONETARY LC_NUMERIC LC_TIME TZ MAKEFLAGS R_LIBS R_LIBS_USER R_LIBS_SITE
R_ENVIRON_USER R_PROFILE_USER R_MAKEVARS_USER R_DATATABLE_NUM_THREADS
OMP_NUM_THREADS OMP_THREAD_LIMIT OMP_DYNAMIC OPENBLAS_NUM_THREADS MKL_NUM_THREADS
MKL_DYNAMIC VECLIB_MAXIMUM_THREADS BLIS_NUM_THREADS NUMEXPR_NUM_THREADS
""".split())


class GuardError(RuntimeError):
    def __init__(self, reason, detail):
        super().__init__(detail)
        self.reason = reason


@dataclass(frozen=True)
class Failure:
    reason: str
    detail: str

    @classmethod
    def from_error(cls, error, default):
        if isinstance(error, GuardError):
            reason, detail = error.reason, str(error)
        else:
            reason, detail = default, f"{type(error).__name__}: {error}"
        return cls(reason, detail[:512] + ("…" if len(detail) > 512 else ""))


@dataclass(frozen=True)
class GuardConfig:
    root: Path
    admitted_utc: datetime
    wall_seconds: float
    rss_bytes: int
    disk_bytes: int
    floor_bytes: int
    sample_seconds: float = 1.0
    disk_seconds: float = 5.0
    ps_timeout: float = 0.8
    disk_timeout: float = 2.0
    disk_max_entries: int = 250_000
    term_seconds: float = 2.0
    kill_seconds: float = 2.0


@dataclass(frozen=True)
class PhaseOutcome:
    status: str
    closure: str
    raw_exit: int | None
    primary: Failure | None
    cleanup_errors: tuple[Failure, ...]
    receipt_error: Failure | None


@dataclass(frozen=True)
class _Member:
    pid: int
    pgid: int
    sid: int
    ruid: int
    euid: int
    state: str
    rss_bytes: int
    sid_source: str


class _PhaseChild(subprocess.Popen):
    def __del__(self):
        # Only the supervisor may reap this child. Popen's destructor polls.
        pass


class _Leader:
    def __init__(self, child):
        self.child = child
        self.pid = child.pid
        self.exit = None
        self.reaped = False
        self.waitid_code = self.waitid_status = self.wait4_status = None
        self.wait_retry_unsafe = False
        self.verified_pids = frozenset()
        self.exiting_descendant = None

    def observe(self, *, strict=False):
        try:
            result = os.waitid(os.P_PID, self.pid, os.WEXITED | os.WNOWAIT | os.WNOHANG)
        except ChildProcessError as error:
            self.wait_retry_unsafe = True
            raise GuardError("anchor_echild", f"waitid lost owned child {self.pid} (ECHILD)") from error
        except OSError as error:
            raise GuardError("anchor_wait", f"waitid for {self.pid}: {error}") from error
        if result is None or result.si_pid == 0:
            if self.exit is not None:
                self.wait_retry_unsafe = True
                raise GuardError("anchor_identity", "Previously waitable child is no longer waitable")
            return
        if result.si_pid != self.pid:
            self.wait_retry_unsafe = True
            raise GuardError("anchor_identity", f"waitid returned unexpected PID {result.si_pid}")
        if result.si_code == os.CLD_EXITED:
            observed_exit = result.si_status
        elif result.si_code in (os.CLD_KILLED, os.CLD_DUMPED):
            observed_exit = -result.si_status
        else:
            self.wait_retry_unsafe = True
            raise GuardError("anchor_wait", f"Unexpected waitid code {result.si_code}")
        if (self.exit is not None
                and (self.exit != observed_exit or (strict and self.waitid_code != result.si_code))):
            self.wait_retry_unsafe = True
            raise GuardError("anchor_identity", "Previously observed child status changed")
        self.exit = observed_exit
        self.waitid_code, self.waitid_status = result.si_code, result.si_status

    def reap(self, *, strict=False):
        if self.exit is None or self.reaped:
            raise GuardError("anchor_reap", "Reap requires one waitable, unreaped leader")
        if strict:
            # The syscall can consume the status even if delivery is interrupted.
            self.wait_retry_unsafe = True
        try:
            pid, status, _usage = os.wait4(self.pid, os.WNOHANG)
        except ChildProcessError as error:
            raise GuardError("anchor_echild", f"wait4 lost owned child {self.pid} (ECHILD)") from error
        if pid != self.pid:
            raise GuardError("anchor_reap", "wait4 did not reap the observed leader")
        self.reaped = True
        self.wait4_status = status
        if strict:
            if os.WIFEXITED(status):
                code = os.CLD_EXITED
            elif os.WIFSIGNALED(status):
                code = os.CLD_DUMPED if os.WCOREDUMP(status) else os.CLD_KILLED
            else:
                raise GuardError("anchor_reap", "wait4 did not return terminal status")
            if code != self.waitid_code:
                raise GuardError("anchor_reap", "wait4 terminal kind disagrees with waitid")
        self.child.returncode = os.waitstatus_to_exitcode(status)
        if self.child.returncode != self.exit:
            raise GuardError("anchor_reap", "wait4 status disagrees with waitid status")


def _absolute(path):
    path = Path(path)
    if not path.is_absolute() or path.anchor != "/" or ".." in path.parts:
        raise GuardError("unsafe_path", f"Expected an absolute path without parent traversal: {path}")
    return path


def _safe_directory(info, path, uid):
    if info.st_uid not in (0, uid):
        raise GuardError("unsafe_path", f"Directory has untrusted owner UID {info.st_uid}: {path}")
    writable = info.st_mode & 0o022
    if writable and not (info.st_uid == 0 and info.st_mode & stat.S_ISVTX):
        raise GuardError("unsafe_path", f"Untrusted writable directory: {path}")


def _open_absolute(path, uid):
    with ExitStack() as handles:
        fd = os.open("/", _DIR_FLAGS)
        handles.callback(os.close, fd)
        _safe_directory(os.fstat(fd), Path("/"), uid)
        for name in path.parts[1:]:
            fd = os.open(name, _DIR_FLAGS, dir_fd=fd)
            handles.callback(os.close, fd)
            _safe_directory(os.fstat(fd), path, uid)
        return os.dup(fd)


@dataclass(frozen=True)
class OwnedCommandResult:
    command: tuple[str, ...]
    pid: int
    waitid_code: int
    waitid_status: int
    wait4_status: int
    exit_code: int
    stdout: bytes
    stderr: bytes

    def kernel_evidence(self):
        return dict(command=self.command, pid=self.pid, waitid_code=self.waitid_code,
                    waitid_status=self.waitid_status, wait4_status=self.wait4_status,
                    exit_code=self.exit_code)


class OwnedCommandError(GuardError):
    def __init__(self, primary, command, leader, cleanup_errors):
        super().__init__(primary.reason, primary.detail)
        self.primary = primary
        self.cleanup_errors = tuple(cleanup_errors)
        self.wait = dict(command=command, pid=leader.pid, waitid_code=leader.waitid_code,
                         waitid_status=leader.waitid_status, wait4_status=leader.wait4_status)


def _read_owned_pipes(streams, deadline):
    output = [bytearray(), bytearray()]
    size = 0
    with selectors.DefaultSelector() as selector:
        for index, stream in enumerate(streams):
            selector.register(stream, selectors.EVENT_READ, index)
        while selector.get_map():
            remaining = deadline - time.monotonic()
            if remaining <= 0 or not (events := selector.select(remaining)):
                raise GuardError("owned_timeout", "Owned command exceeded its pipe-read deadline")
            for key, _events in events:
                block = os.read(key.fd, min(4096, _PS_OUTPUT_BYTES - size + 1))
                if not block:
                    selector.unregister(key.fileobj)
                    continue
                size += len(block)
                if size > _PS_OUTPUT_BYTES:
                    raise GuardError("owned_output", "Owned command exceeded its pipe output budget")
                output[key.data].extend(block)
    return bytes(output[0]), bytes(output[1])


def _wait_owned(leader, deadline):
    while (remaining := deadline - time.monotonic()) > 0:
        leader.observe(strict=True)
        if leader.exit is not None:
            leader.reap(strict=True)
            return
        time.sleep(min(0.01, remaining))
    raise GuardError("owned_timeout", "Owned command lacked a terminal status before its deadline")


def _clean_owned(leader, cleanup):
    try:
        deadline = time.monotonic() + _OWNED_CLEANUP_SECONDS
        leader.observe(strict=True)
        if leader.exit is None:
            try:
                os.kill(leader.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass  # Exit still requires an exact kernel observation and reap.
            except OSError as error:
                cleanup.append(Failure("owned_cleanup", f"Signaling owned PID {leader.pid}: {error}"))
            _wait_owned(leader, deadline)
        else:
            leader.reap(strict=True)
    except BaseException as error:
        cleanup.append(Failure("owned_cleanup", f"Owned PID {leader.pid} closure unconfirmed: {error}"))


def _run_owned_command(command, *, cwd, stdout, timeout, env):
    command = tuple(command)
    deadline = time.monotonic() + timeout
    capture = stdout == subprocess.PIPE
    child = _PhaseChild(command, cwd=cwd, stdout=stdout, close_fds=True, env=env,
                        stdin=subprocess.DEVNULL if capture else None,
                        stderr=subprocess.PIPE if capture else None)
    leader = _Leader(child)
    streams = (child.stdout, child.stderr) if capture else ()
    output = (b"", b"")
    primary = None
    cleanup = []
    verified = False
    try:
        if type(leader.pid) is not int or leader.pid <= 1:
            leader.wait_retry_unsafe = True
            raise GuardError("anchor_identity", "Spawn did not yield a positive owned PID")
        if capture:
            output = _read_owned_pipes(streams, deadline)
        _wait_owned(leader, deadline)
        verified = True
    except BaseException as error:
        primary = error
    finally:
        for stream in streams:
            try:
                stream.close()
            except BaseException as error:
                if primary is None and not isinstance(error, OSError):
                    primary = error
                else:
                    cleanup.append(Failure("owned_cleanup", f"Closing owned command pipe: {error}"))
        if primary is not None and not verified:
            if not leader.wait_retry_unsafe:
                _clean_owned(leader, cleanup)
            else:
                # Lost, contradictory or possibly consumed evidence is final.
                cleanup.append(Failure("owned_cleanup", f"Owned PID {leader.pid} closure unconfirmed: {primary}"))

    if primary is not None and not isinstance(primary, (GuardError, OSError)):
        if isinstance(primary, KeyboardInterrupt):
            primary._owned_command_error = OwnedCommandError(
                Failure.from_error(primary, "owned_interrupted"), command, leader, cleanup)
        primary._guard_cleanup_errors = tuple(Failure("ps_helper_cleanup", item.detail) for item in cleanup)
        raise primary
    if primary is not None or cleanup:
        failure = Failure.from_error(primary, "owned_io") if primary is not None else cleanup.pop(0)
        raise OwnedCommandError(failure, command, leader, cleanup) from primary
    return OwnedCommandResult(command, leader.pid, leader.waitid_code, leader.waitid_status,
                              leader.wait4_status, leader.exit, *output)


def run_owned_command(command, *, cwd, stdout, timeout):
    """Run one direct child, returning its original parent's exact kernel proof."""
    env = {"PATH": "/usr/bin:/bin", "LC_ALL": "C"} if stdout == subprocess.PIPE else None
    return _run_owned_command(command, cwd=cwd, stdout=stdout, timeout=timeout, env=env)


def _query_ps(pgid, timeout):
    """Read a bounded group query and trust only an observed and reaped exit."""
    argv = ["/bin/ps", "-x", "-g", str(pgid), "-o", "pid=,pgid=,ruid=,uid=,state=,rss="]
    try:
        result = _run_owned_command(argv, cwd=None, stdout=subprocess.PIPE, timeout=timeout, env=_PS_ENV)
    except OwnedCommandError as error:
        reason = {"owned_timeout": "ps_timeout", "owned_output": "ps_truncated",
                  "owned_cleanup": "ps_helper_cleanup"}.get(error.reason, "ps_unavailable")
        failure = GuardError(reason, str(error))
        failure._guard_cleanup_errors = tuple(Failure("ps_helper_cleanup", item.detail)
                                              for item in error.cleanup_errors)
        raise failure from error
    return result.exit_code, result.stdout, result.stderr

class _ExitingLeaderRow(GuardError):
    """One exact-ID anchor row with Darwin's pre-zombie E and no usable RSS."""


class _ExitingDescendantRow(GuardError):
    """A previously verified descendant row, plus the separately strict rows."""

    def __init__(self, pid, detail):
        super().__init__("ps_parse", detail)
        self.pid = pid


def _ps_parse_detail(reason, line):
    # Only include a literal row if every character could occur in these six
    # selected fields. Keep even malformed numeric rows useful without logging
    # an arbitrary string. Capture a bounded prefix and its original length.
    excerpt = line[:120]
    if not re.fullmatch(r"[0-9 \t?+<>AELNSsVWXIRSTUZH-]*", excerpt, re.ASCII):
        fields = line.split()
        safe = []
        for i, field in enumerate(fields[:6]):
            pattern = r"[?+<>AELNSsVWXIRSTUZH-]{1,16}" if i == 4 else r"[0-9-]{1,20}"
            safe.append(field if re.fullmatch(pattern, field, re.ASCII) else "<invalid>")
        return f"{reason}; fields {safe!r} ({len(fields)} fields, {len(line)} characters)"
    return f"{reason}; ps row {excerpt!r} ({len(line)} characters)"


def _parse_ps(code, out, err, *, pgid, uid, verified_pids=()):
    if len(out) + len(err) > _PS_OUTPUT_BYTES or (out and not out.endswith(b"\n")):
        raise GuardError("ps_truncated", "Incomplete owned group output")
    if code != 0 or err or not out:
        detail = err[:256].decode("ascii", errors="backslashreplace").strip()
        raise GuardError("ps_unavailable", f"Owned group query exit {code}, stderr {detail!r}")
    try:
        lines = out.decode("ascii").splitlines()
    except UnicodeDecodeError as error:
        raise GuardError("ps_parse", "Owned group output is not ASCII") from error
    rows = {}
    transition = None
    for line in lines:
        fields = line.split()
        if len(fields) != 6:
            raise GuardError("ps_parse", _ps_parse_detail("Expected six owned group fields", line))
        for name, field in zip(("PID", "PGID", "RUID", "EUID"), fields[:4]):
            if not _DECIMAL.fullmatch(field):
                raise GuardError("ps_parse", _ps_parse_detail(f"Invalid owned group {name}", line))
        if not _STATE.fullmatch(fields[4]):
            detail = _ps_parse_detail("Invalid owned group state", line)
            if (len(lines) == 1 and fields[:4] == [str(pgid), str(pgid), str(uid), str(uid)]
                    and fields[4] in ("?E", "?Es", "HE", "HEs") and fields[5] in ("0", "-")):
                raise _ExitingLeaderRow("ps_parse", detail)
            pid = int(fields[0])
            if (fields[4:] == ["?E", "0"] and pid != pgid and pid in verified_pids
                    and fields[1:4] == [str(pgid), str(uid), str(uid)]
                    and transition is None and pid not in rows):
                transition = _ExitingDescendantRow(pid, detail)
                continue
            raise GuardError("ps_parse", detail)
        pid, group, ruid, euid = map(int, fields[:4])
        state, rss = fields[4:]
        if rss == "-" and state.startswith("Z"):
            rss = "0"
        if not _DECIMAL.fullmatch(rss):
            raise GuardError("ps_parse", _ps_parse_detail("Invalid owned group RSS field", line))
        if (pid <= 1 or pid in rows or (transition is not None and pid == transition.pid)
                or group != pgid or ruid != uid or euid != uid):
            raise GuardError("group_identity", f"Unverified PID/PGID/UID for row PID {pid}")
        rows[pid] = (group, ruid, euid, state, int(rss) * 1024)
    if transition is not None:
        # Complete parsing first, including every row after the invalid one.
        transition.strict_rows = rows
        raise transition
    return rows


def _transition_evidence(leader, captured, uid):
    """Keep obligations and attributable RSS even when the strict parser refuses this snapshot."""
    code, out, err = captured
    if (code or err or not out.endswith(b"\n") or len(out) > _PS_OUTPUT_BYTES
            or not out.isascii() or len(leader.verified_pids) < 2):
        return 0
    if leader.exiting_descendant is None and b"?E" not in out:
        return 0

    previous = {str(pid).encode(): pid for pid in leader.verified_pids}
    identity = [str(leader.pid).encode(), str(uid).encode(), str(uid).encode()]
    known = {}
    duplicate = set()
    for line in out.splitlines():
        fields = line.split()
        if not fields or not fields[0].isdigit():
            continue
        pid = previous.get(fields[0].lstrip(b"0"))
        if pid is None:
            continue
        if pid in known:
            duplicate.add(pid)
        known[pid] = line
        if (pid != leader.pid and leader.exiting_descendant is None
                and fields[1:] == [*identity, b"?E", b"0"]):
            # A conservative cleanup obligation only. Other rows can still veto
            # the entire sample. Do not depend on which line the parser sees first.
            leader.exiting_descendant = pid
    if leader.exiting_descendant is None:
        return 0
    minimum = 0
    for pid, line in known.items():
        if pid == leader.exiting_descendant or pid in duplicate:
            continue
        try:
            rows = _parse_ps(0, line + b"\n", b"", pgid=leader.pid, uid=uid)
        except (GuardError, ValueError):
            continue
        _, _, _, state, rss = rows[pid]
        if not state.startswith("Z"):
            minimum += rss
    return minimum


class Guard:
    """One admitted run with sequential, independently owned phases."""

    def __init__(self, config):
        self.config = config
        if sys.platform != "darwin" or os.getuid() != os.geteuid():
            raise GuardError("unsupported_host", "Darwin and unchanged effective UID are required")
        self.uid = os.geteuid()
        self._validate_config()
        self.root = _absolute(config.root)
        checked = time.monotonic()
        elapsed = (datetime.now(timezone.utc) - config.admitted_utc).total_seconds()
        if elapsed < 0 or elapsed >= config.wall_seconds:
            raise GuardError("stale_admission", "Admission is in the future or its wall window has elapsed")
        self.deadline = checked + config.wall_seconds - elapsed
        try:
            fd = _open_absolute(self.root, self.uid)
        except OSError as error:
            raise GuardError("unsafe_path", f"Cannot open run root without symlinks: {error}") from error
        try:
            info = os.fstat(fd)
            if self.root == Path("/") or info.st_uid != self.uid or info.st_mode & 0o077:
                raise GuardError("unsafe_path", "Run root must be a private same-user directory")
            self._root_id = (info.st_dev, info.st_ino)
        finally:
            os.close(fd)
        self._shutdown = None
        self._blocked = False
        self._run_lock = threading.Lock()

    def _validate_config(self):
        config = self.config
        if not isinstance(config.admitted_utc, datetime) or (
                config.admitted_utc.tzinfo is None
                or config.admitted_utc.utcoffset() != timezone.utc.utcoffset(None)):
            raise GuardError("invalid_config", "Admission must be a timezone-aware UTC datetime")
        seconds = ((config.wall_seconds, _MAX_WALL_SECONDS), (config.sample_seconds, 1.0),
                   (config.disk_seconds, 10.0), (config.ps_timeout, 2.0),
                   (config.disk_timeout, 3.0), (config.term_seconds, 2.0), (config.kill_seconds, 2.0))
        if any(type(value) not in (int, float) or not math.isfinite(value) or not 0 < value <= maximum
               for value, maximum in seconds):
            raise GuardError("invalid_config", "A deadline or monitoring interval is outside its bound")
        limits = (config.rss_bytes, config.disk_bytes, config.floor_bytes, config.disk_max_entries)
        if (any(type(value) is not int or value <= 0 for value in limits)
                or config.rss_bytes > _MAX_RSS_BYTES or config.disk_max_entries > 250_000):
            raise GuardError("invalid_config", "A resource limit is outside its bound")

    @contextmanager
    def shutdown_signals(self):
        """Request an owned shutdown on INT, TERM or HUP in the main thread."""
        with ExitStack() as handlers:
            for signum in _STOP_SIGNALS:
                previous = signal.signal(signum, self._signal_handler)
                handlers.callback(signal.signal, signum, previous)
            yield self

    def _signal_handler(self, signum, _frame):
        if self._shutdown is None:
            self._shutdown = signum

    def _stop_check(self):
        if self._shutdown is not None:
            raise GuardError("shutdown", f"Supervisor received signal {self._shutdown}")
        if time.monotonic() >= self.deadline:
            raise GuardError("wall", "Shared monotonic admission deadline reached")

    @contextmanager
    def _directory(self, path):
        path = _absolute(path)
        if not path.is_relative_to(self.root):
            raise GuardError("unsafe_path", f"Phase directory is outside the run root: {path}")
        with ExitStack() as handles:
            try:
                fd = _open_absolute(self.root, self.uid)
                handles.callback(os.close, fd)
                info = os.fstat(fd)
                if ((info.st_dev, info.st_ino) != self._root_id or info.st_uid != self.uid
                        or info.st_mode & 0o077):
                    raise GuardError("unsafe_path", "Admitted run root identity or privacy changed")
                for name in path.relative_to(self.root).parts:
                    fd = os.open(name, _DIR_FLAGS, dir_fd=fd)
                    handles.callback(os.close, fd)
                    info = os.fstat(fd)
                    _safe_directory(info, path, self.uid)
                    if info.st_dev != self._root_id[0]:
                        raise GuardError("wrong_volume", f"Phase directory is on another volume: {path}")
            except OSError as error:
                raise GuardError("unsafe_path", f"Cannot open admitted directory {path}: {error}") from error
            yield fd

    def _new_file(self, path):
        path = _absolute(path)
        if path == self.root or not path.is_relative_to(self.root):
            raise GuardError("unsafe_path", f"Evidence file is outside the run root: {path}")
        with self._directory(path.parent) as dir_fd:
            return os.open(path.name, _FILE_FLAGS, 0o600, dir_fd=dir_fd)

    def _environment(self, env):
        environment = dict(env)
        if any(key not in _ENV_KEYS or not isinstance(value, str) or "\0" in value
               for key, value in environment.items()):
            raise GuardError("unsafe_environment", "Phase environment has an unsupported key or value")
        if not {"HOME", "TMPDIR"} <= environment.keys():
            raise GuardError("unsafe_environment", "Phase environment needs private HOME and TMPDIR")
        for key in ("HOME", "TMPDIR", "TMP", "TEMP"):
            if key in environment:
                with self._directory(environment[key]):
                    pass
        return environment

    def _space(self, *, walk):
        with self._directory(self.root) as fd:
            if walk:
                result = self._disk_metadata(fd)
            else:
                result = {}
            volume = os.fstatvfs(fd)
            fragment = volume.f_frsize
            if fragment <= 0:
                raise GuardError("disk_monitor", "Filesystem reported an invalid allocation fragment size")
            result["free_bytes"] = volume.f_bavail * fragment
            # Darwin's preferred I/O size (f_bsize) can far exceed the
            # allocation fragment charged by this filesystem.
            result["block_bytes"] = max(512, fragment)
        return result

    def _disk_metadata(self, root_fd):
        config = self.config
        deadline = time.monotonic() + config.disk_timeout
        count = logical = allocated = 0

        def charge(info):
            nonlocal count, logical, allocated
            count += 1
            if count > config.disk_max_entries or time.monotonic() >= deadline:
                raise GuardError("disk_monitor", "Owned metadata walk exceeded its work budget")
            if info.st_dev != self._root_id[0]:
                raise GuardError("wrong_volume", "Owned metadata includes a different volume")
            logical += info.st_size
            allocated += info.st_blocks * 512

        def scan(fd, depth):
            if depth > 64:
                raise GuardError("disk_monitor", "Owned metadata walk exceeded its depth budget")
            with os.scandir(fd) as entries:
                for entry in entries:
                    if time.monotonic() >= deadline:
                        raise GuardError("disk_monitor", "Owned metadata walk exceeded its time budget")
                    try:
                        info = entry.stat(follow_symlinks=False)
                    except FileNotFoundError:
                        continue
                    charge(info)
                    if not stat.S_ISDIR(info.st_mode):
                        continue
                    try:
                        child_fd = os.open(entry.name, _DIR_FLAGS, dir_fd=fd)
                    except FileNotFoundError:
                        continue
                    try:
                        actual = os.fstat(child_fd)
                        if actual.st_dev != info.st_dev:
                            raise GuardError("wrong_volume", "Directory moved to another volume during scan")
                        if actual.st_ino != info.st_ino:
                            raise GuardError("disk_monitor", "Directory identity changed during scan")
                        scan(child_fd, depth + 1)
                    finally:
                        os.close(child_fd)

        charge(os.fstat(root_fd))
        scan(root_fd, 0)
        if time.monotonic() >= deadline:
            raise GuardError("disk_monitor", "Owned metadata walk exceeded its time budget")
        return {"logical_bytes": logical, "allocated_bytes": allocated, "entries": count}

    def _space_sample(self, audit, *, walk, evidence_bytes=0, metadata_blocks=0, new_entries=0):
        audit["last_space"] = None
        try:
            sample = self._space(walk=walk)
        except OSError as error:
            raise GuardError("disk_monitor", f"Owned filesystem observation failed: {error}") from error
        audit["last_space"] = sample
        free = sample["free_bytes"]
        low = audit["minimum_free_bytes"]
        audit["minimum_free_bytes"] = free if low is None else min(low, free)
        if walk:
            for field in ("logical_bytes", "allocated_bytes"):
                peak_field = f"peak_disk_{field}"
                audit[peak_field] = max(audit[peak_field], sample[field])
        self._check_space(sample, evidence_bytes=evidence_bytes, metadata_blocks=metadata_blocks,
                          new_entries=new_entries)
        return sample

    def _check_space(self, sample, *, evidence_bytes=0, metadata_blocks=0, new_entries=0):
        free = sample["free_bytes"]
        if free < self.config.floor_bytes:
            raise GuardError("filesystem_floor", f"Available bytes {free} fell below the floor")
        logical = sample.get("logical_bytes", 0)
        allocated = sample.get("allocated_bytes", 0)
        if max(logical, allocated) > self.config.disk_bytes:
            raise GuardError("disk", "Logical or allocated bytes exceeded the owned disk budget")
        if evidence_bytes or metadata_blocks or new_entries:
            block = sample["block_bytes"]
            # Even an empty new file can grow its parent's counted st_size.
            entry_growth = new_entries * block
            logical_charge = evidence_bytes + entry_growth
            allocation = ((evidence_bytes + block - 1) // block + metadata_blocks) * block + entry_growth
            if (free - allocation < self.config.floor_bytes
                    or max(logical + logical_charge, allocated + allocation) > self.config.disk_bytes):
                raise GuardError("evidence_reserve", f"Evidence needs {logical_charge} logical bytes and "
                                 f"{allocation} bytes with allocation and metadata allowance")

    def _check_descendant_transition(self, leader, transition):
        rows = transition.strict_rows
        # The excluded row has *unknown* RSS. A breach in the measured rows
        # alone already suffices, even if the next snapshot would be smaller.
        lower_bound = sum(rss for _, _, _, state, rss in rows.values() if not state.startswith("Z"))
        if lower_bound > self.config.rss_bytes:
            raise GuardError("rss", f"Strict rows alone show at least {lower_bound} RSS bytes "
                             f"against limit {self.config.rss_bytes}; {transition}")

        leader.observe()
        if (leader.exit is not None or leader.pid not in rows or rows[leader.pid][3].startswith("Z")
                or rows.keys() | {transition.pid} != leader.verified_pids):
            raise transition
        for pid in leader.verified_pids:
            try:
                if os.getpgid(pid) != leader.pid or os.getsid(pid) != leader.pid:
                    raise GuardError("group_identity", f"PID {pid} changed group or session during exit")
            except ProcessLookupError:
                if pid != transition.pid:
                    raise GuardError("group_identity", f"Cannot verify companion PID {pid} during exit")
            except OSError as error:
                raise GuardError("group_identity", f"Cannot verify PID {pid} during exit: {error}") from error

    @staticmethod
    def _check_disappeared(leader, rows):
        if leader.exiting_descendant in rows:
            raise GuardError("group_identity", f"Previously exiting descendant {leader.exiting_descendant} remains present")
        for pid in leader.verified_pids - rows.keys():
            try:
                pgid = os.getpgid(pid)
            except ProcessLookupError:
                continue
            except OSError as error:
                raise GuardError("group_identity", f"Cannot verify missing PID {pid}: {error}") from error
            raise GuardError("group_identity", f"Missing PID {pid} still resolves in group {pgid}")

    def _snapshot(self, leader, deadline):
        """Validate one group snapshot with at most one bounded disappearance recheck."""
        # Store only the current query's lower bound on this per-call context.
        # This also covers refusals after parsing, such as a still-visible PID.
        observed = {"minimum_rss": 0}
        try:
            return self._bounded_snapshot(leader, deadline, observed)
        except GuardError as error:
            minimum = observed["minimum_rss"]
            if minimum <= self.config.rss_bytes or error.reason == "rss":
                raise
            failure = GuardError("rss", f"Strict prior-member rows alone show at least {minimum} RSS bytes "
                                 f"against limit {self.config.rss_bytes}; {error.reason}: {error}")
            failure._guard_cleanup_errors = getattr(error, "_guard_cleanup_errors", ())
            raise failure from error

    def _bounded_snapshot(self, leader, deadline, observed):
        end = min(deadline, time.monotonic() + self.config.ps_timeout)
        for attempt in range(2):
            observed["minimum_rss"] = 0
            leader.observe()
            remaining = end - time.monotonic()
            if remaining <= 0:
                raise GuardError("ps_timeout", "No time left for a bounded owned group snapshot")
            try:
                eligible = leader.verified_pids if not attempt and leader.exiting_descendant is None else ()
                captured = _query_ps(leader.pid, remaining)
                observed["minimum_rss"] = _transition_evidence(leader, captured, self.uid)
                rows = _parse_ps(*captured, pgid=leader.pid, uid=self.uid, verified_pids=eligible)
            except _ExitingDescendantRow as transition:
                # Keep this obligation across a failed snapshot and cleanup.
                # Only a fully verified fresh snapshot can discharge it.
                leader.exiting_descendant = transition.pid
                self._check_descendant_transition(leader, transition)
                # Give the real parent a bounded opportunity to reap. The sleep
                # proves nothing. Only a new full group query can supply a sample.
                time.sleep(min(0.10, max(0, end - time.monotonic()) / 4))
                continue
            except _ExitingLeaderRow:
                if attempt:
                    raise
                # Never use a pre-zombie row as a member, sample or authority
                # to signal. Wait for the direct child without releasing its PID.
                # Spend only this snapshot's existing time and requery allowance.
                leader.observe()
                while leader.exit is None and (remaining := end - time.monotonic()) > 0:
                    time.sleep(min(0.01, remaining))
                    leader.observe()
                if leader.exit is None or time.monotonic() >= end:
                    raise
                continue
            except OSError as error:
                failure = GuardError("ps_unavailable", f"Cannot query owned group: {error}")
                failure._guard_cleanup_errors = getattr(error, "_guard_cleanup_errors", ())
                raise failure from error
            leader.observe()
            if leader.exiting_descendant is not None:
                self._check_disappeared(leader, rows)
            anchor = rows.get(leader.pid)
            if anchor is None:
                raise GuardError("group_identity", "Owned group query omitted its unreaped leader")
            if anchor[3].startswith("Z") != (leader.exit is not None):
                continue
            members = []
            for pid, (pgid, ruid, euid, state, rss) in rows.items():
                sid_source = "getsid"
                try:
                    if os.getpgid(pid) != pgid:
                        raise GuardError("group_identity", f"PID {pid} changed group")
                    sid = os.getsid(pid)
                except ProcessLookupError:
                    if pid != leader.pid or leader.exit is None:
                        break
                    # A waitable direct child cannot have its PID reused. The
                    # verified zombie row and successful setsid bind its SID.
                    sid = leader.pid
                    sid_source = "owned setsid plus waitable zombie"
                except OSError as error:
                    raise GuardError("group_identity", f"Cannot verify PID {pid} identity: {error}") from error
                if sid != leader.pid:
                    raise GuardError("group_identity", f"PID {pid} changed session")
                members.append(_Member(pid, pgid, sid, ruid, euid, state, rss, sid_source))
            else:
                if leader.exiting_descendant is not None and time.monotonic() >= end:
                    raise GuardError("ps_timeout", "Descendant exit verification exceeded its snapshot budget")
                leader.exiting_descendant = None
                leader.verified_pids = frozenset(member.pid for member in members)
                return members
            if not attempt:
                time.sleep(min(0.02, max(0, end - time.monotonic())))
        raise GuardError("group_race", "Owned process identity did not settle after one recheck")

    @staticmethod
    def _closed(leader, members):
        return leader.exit is not None and [member.pid for member in members] == [leader.pid]

    def _sample_group(self, audit, members, started):
        rss = sum(member.rss_bytes for member in members if not member.state.startswith("Z"))
        sample = dict(monotonic=started, completed_monotonic=time.monotonic(),
                      rss_bytes=rss, members=[asdict(member) for member in members])
        if audit["first_group"] is None:
            audit["first_group"] = sample
        audit["group_samples"] += 1
        audit["last_group"] = sample
        if rss >= audit["peak_sampled_group_rss_bytes"]:
            audit["peak_sampled_group_rss_bytes"] = rss
            audit["peak_group"] = sample
        if rss > self.config.rss_bytes:
            raise GuardError("rss", f"Observed group RSS {rss} exceeded its byte limit")

    def _monitor(self, leader, audit):
        previous = None
        missed = 0
        next_disk = time.monotonic() + self.config.disk_seconds
        while True:
            self._stop_check()
            leader.observe()
            started = time.monotonic()
            if (leader.exit is not None or previous is None
                    or started >= previous + self.config.sample_seconds):
                members = self._snapshot(leader, self.deadline)
                closed = self._closed(leader, members)
                if closed:
                    audit["terminal_members"] = [asdict(member) for member in members]
                self._sample_group(audit, members, started)
                if previous is not None:
                    gap = started - previous
                    audit["max_sample_gap_seconds"] = max(gap, audit["max_sample_gap_seconds"])
                    if gap > 2 * self.config.sample_seconds:
                        missed += 1
                        audit["missed_samples"] += 1
                    else:
                        missed = 0
                    if missed >= 2:
                        raise GuardError("sample_cadence", "Two consecutive owned group sample gaps exceeded tolerance")
                if closed:
                    return
                if leader.exit is not None:
                    raise GuardError("leftover_descendants", "Leader exited with owned descendants remaining")
                previous = started
                walk = started >= next_disk
                self._space_sample(audit, walk=walk, evidence_bytes=_TERMINAL_BYTES, metadata_blocks=1)
                if walk:
                    next_disk = started + self.config.disk_seconds
            time.sleep(min(0.05, self.config.sample_seconds / 4,
                           max(0, self.deadline - time.monotonic())))

    @staticmethod
    def _signal(leader, signum, audit, members=None):
        leader.observe()
        if members is None:
            if leader.exit is not None:
                return
            try:
                identity = os.getpgid(leader.pid), os.getsid(leader.pid)
            except ProcessLookupError as error:
                leader.observe()
                if leader.exit is not None:
                    return
                raise GuardError("anchor_identity", "Direct child identity vanished before signaling") from error
            if identity != (leader.pid, leader.pid):
                raise GuardError("group_identity", "Cannot verify direct child identity before signaling")
        event = {"signal": int(signum), "target": "leader" if members is None else "group",
                 "pid": leader.pid, "monotonic": time.monotonic(), "result": "sent"}
        audit["signals"].append(event)
        try:
            if members is None:
                os.kill(leader.pid, signum)
            else:
                event["members"] = [asdict(member) for member in members]
                os.killpg(leader.pid, signum)
        except ProcessLookupError:
            event["result"] = "ESRCH"
        except OSError as error:
            event["result"] = f"errno={error.errno}"
            raise GuardError("signal_error", f"Owned signal {signum} failed: {error}") from error

    def _fallback_leader(self, leader, audit, errors, deadline):
        """On lost group metadata, signal only the separately anchored direct child."""
        for signum in (signal.SIGTERM, signal.SIGKILL):
            if time.monotonic() >= deadline:
                return
            try:
                leader.observe()
                if leader.exit is not None:
                    return
                self._signal(leader, signum, audit)
                end = deadline
                if signum == signal.SIGTERM:
                    end = min(end, time.monotonic() + self.config.term_seconds / 2)
                while time.monotonic() < end:
                    leader.observe()
                    if leader.exit is not None:
                        return
                    time.sleep(min(0.05, max(0, end - time.monotonic())))
            except (GuardError, OSError) as error:
                errors.append(Failure.from_error(error, "leader_cleanup"))
                return

    def _clean_group(self, leader, audit, errors):
        end = time.monotonic() + self.config.term_seconds + self.config.kill_seconds
        for signum, grace in ((signal.SIGTERM, self.config.term_seconds),
                              (signal.SIGKILL, self.config.kill_seconds)):
            stage_end = min(end, time.monotonic() + grace)
            try:
                members = self._snapshot(leader, stage_end)
                if self._closed(leader, members):
                    audit["terminal_members"] = [asdict(member) for member in members]
                    return
                self._signal(leader, signum, audit, members)
                # Use the stage deadline both for waiting and observing. Escalate
                # from a new verified snapshot, never from the old member list.
                while time.monotonic() + min(0.10, self.config.ps_timeout) < stage_end:
                    time.sleep(min(0.05, max(0, stage_end - time.monotonic())))
                    members = self._snapshot(leader, stage_end)
                    if self._closed(leader, members):
                        audit["terminal_members"] = [asdict(member) for member in members]
                        return
            except BaseException as error:
                errors.append(Failure.from_error(error, "group_cleanup"))
                errors.extend(getattr(error, "_guard_cleanup_errors", ()))
                self._fallback_leader(leader, audit, errors, end)
                return
        errors.append(Failure("group_unclosed", "Owned group and leader were not confirmed closed in time"))

    @staticmethod
    def _encode_receipt(record):
        status = record["status"]
        limit = {"planned": _PLAN_BYTES, "running": _RUNNING_BYTES}.get(status, _TERMINAL_BYTES)
        if status == "not_started" and (record.get("primary") or {}).get("reason") == "unresolved_phase":
            limit = _PENDING_BYTES
        encoded = (json.dumps(record, sort_keys=True, allow_nan=False) + "\n").encode()
        if len(encoded) > limit:
            reason = "evidence_input" if status == "planned" else "evidence_size"
            raise GuardError(reason, f"{status} receipt exceeds its {limit}-byte limit")
        return encoded

    @staticmethod
    def _write_receipt(fd, record):
        encoded = Guard._encode_receipt(record)
        if os.write(fd, encoded) != len(encoded):
            raise OSError("Short receipt write")

    def _phase_receipt(self, fd, record):
        try:
            self._write_receipt(fd, record)
        except (OSError, ValueError, TypeError) as error:
            raise GuardError("receipt_write", f"{record['status']} receipt failed: {error}") from error

    @staticmethod
    def _phase_arguments(argv, allowed_exits):
        if isinstance(argv, (str, bytes, Path)):
            raise GuardError("invalid_phase", "Phase argv must be a sequence of arguments")
        args = []
        remaining = _ARG_BYTES
        for value in islice(argv, _MAX_ARGUMENTS + 1):
            arg = os.fspath(value)
            if not isinstance(arg, str) or "\0" in arg:
                raise GuardError("invalid_phase", "Phase argv requires text without NUL characters")
            if len(args) == _MAX_ARGUMENTS or len(arg) >= remaining:
                raise GuardError("evidence_input", "Phase argv exceeded its count or 8192-byte evidence limit")
            remaining -= len(arg.encode("utf-8", errors="surrogatepass")) + 1
            if remaining < 0:
                raise GuardError("evidence_input", "Phase argv exceeded its 8192-byte evidence limit")
            args.append(arg)
        if not args or not Path(args[0]).is_absolute():
            raise GuardError("invalid_phase", "Phase requires nonempty argv and an absolute executable")
        _absolute(args[0])
        accepted = tuple(islice(allowed_exits, 257))
        if (not accepted or len(accepted) > 256
                or any(type(code) is not int or not 0 <= code <= 255 for code in accepted)):
            raise GuardError("invalid_phase", "Allowed raw exits must be ordinary exit statuses")
        return tuple(args), accepted

    def _finish_receipt(self, fd, record, audit):
        """Write at most one funded terminal record, with a compact failure when it fits."""
        sample = audit["last_space"]
        if sample is None:
            return Failure("receipt_reserve", "Terminal evidence lacks a current space observation")
        failure = None
        try:
            encoded = self._encode_receipt(record)
            self._check_space(sample, evidence_bytes=len(encoded), metadata_blocks=1)
        except (GuardError, ValueError, TypeError) as error:
            reason = "receipt_reserve" if isinstance(error, GuardError) and error.reason in (
                "disk", "filesystem_floor", "evidence_reserve") else "terminal_receipt"
            failure = Failure(reason, Failure.from_error(error, reason).detail)
            summary = audit["trigger_sample"]
            if summary is not None and summary["group"] is not None:
                group = summary["group"]
                summary = {**summary, "group": {"monotonic": group["monotonic"], "rss_bytes": group["rss_bytes"]}}
            record = {
                key: record[key] for key in ("label", "closure", "raw_exit", "leader_reaped", "utc",
                                             "primary", "cleanup_errors", "elapsed_seconds")
            }
            record.update(status="not_started" if record["closure"] == "not_started" else "error",
                          receipt_error=asdict(failure), detailed_evidence_omitted=True, trigger_sample=summary)
            try:
                compact = self._encode_receipt(record)
                self._check_space(sample, evidence_bytes=len(compact), metadata_blocks=1)
            except (GuardError, ValueError, TypeError):
                return failure
        try:
            self._write_receipt(fd, record)
        except (GuardError, OSError, ValueError, TypeError) as error:
            write_failure = Failure("terminal_receipt", Failure.from_error(error, "terminal_receipt").detail)
            if failure is not None:
                return Failure(failure.reason, f"{failure.detail}; {write_failure.detail}")
            return write_failure
        return failure

    def run_phase(self, label, argv, *, cwd, env, stdout, stderr, receipt, allowed_exits=(0,)):
        """Run an explicit argv with exclusive evidence, returning its exact outcome."""
        if not self._run_lock.acquire(blocking=False):
            raise GuardError("overlapping_phase", "An owned phase is already running")
        try:
            return self._run_phase(label, argv, cwd, env, stdout, stderr, receipt, allowed_exits)
        finally:
            self._run_lock.release()

    def _run_phase(self, label, argv, cwd, env, stdout, stderr, receipt, allowed_exits):
        started = time.monotonic()
        audit = {
            "first_group": None, "last_group": None, "peak_group": None, "terminal_members": None,
            "prelaunch_rss_bytes": 0, "prelaunch_space": None, "last_space": None,
            "trigger_sample": None, "minimum_free_bytes": None, "peak_disk_logical_bytes": 0,
            "peak_disk_allocated_bytes": 0, "peak_sampled_group_rss_bytes": 0,
            "group_samples": 0, "max_sample_gap_seconds": 0, "missed_samples": 0, "signals": [],
        }
        leader = receipt_error = None
        pending = self._blocked
        primary = (Failure("unresolved_phase", "Earlier phase closure or evidence is unresolved")
                   if pending else None)
        errors = []
        closure = "not_started"
        try:
            # Fund the entry before opening it. For a pending phase, also
            # reserve the entire bounded record and its append metadata.
            self._space_sample(audit, walk=True, new_entries=1,
                               evidence_bytes=_PENDING_BYTES if pending else 0,
                               metadata_blocks=1 if pending else 0)
        except (GuardError, OSError) as error:
            self._blocked = True
            space_error = Failure.from_error(error, "disk_monitor")
            receipt_error = Failure("receipt_reserve", space_error.detail)
            if primary is None:
                primary = space_error
            return PhaseOutcome("not_started", closure, None, primary, (), receipt_error)
        try:
            receipt_fd = self._new_file(receipt)
        except (GuardError, OSError, TypeError, ValueError) as error:
            self._blocked = True
            receipt_error = Failure("receipt_open", Failure.from_error(error, "receipt_open").detail)
            return PhaseOutcome("not_started", closure, None, primary, (), receipt_error)

        receipt_label = "[invalid phase label]"
        try:
            try:
                if (not isinstance(label, str) or not label or len(label) > _MAX_LABEL_BYTES or "\0" in label
                        or len(label.encode("utf-8", errors="surrogatepass")) > _MAX_LABEL_BYTES):
                    raise GuardError("evidence_input", "Phase label must fit the 256-byte evidence limit")
                receipt_label = label
                if primary is not None:
                    # This phase is already blocked. Only its not_started
                    # record needs funding, with no plan, logs, or child.
                    raise GuardError(primary.reason, primary.detail)
                args, accepted = self._phase_arguments(argv, allowed_exits)
                paths = {name: os.fspath(path) for name, path in (
                    ("cwd", cwd), ("stdout", stdout), ("stderr", stderr), ("receipt", receipt),
                )}
                if any(not isinstance(path, str) or len(path) > 4096 for path in paths.values()):
                    raise GuardError("evidence_input", "Phase path text exceeded its evidence bound")
                plan = {
                    "status": "planned", "label": label, "argv": args,
                    **paths, "utc": datetime.now(timezone.utc).isoformat(),
                    "allowed_exits": accepted, "admitted_utc": self.config.admitted_utc.isoformat(),
                    "monotonic_deadline": self.deadline,
                    "limits": {key: value for key, value in asdict(self.config).items()
                               if key not in {"root", "admitted_utc"}},
                }
                planned_bytes = len(self._encode_receipt(plan))
                self._space_sample(audit, walk=True, evidence_bytes=planned_bytes + _RUNNING_BYTES + _TERMINAL_BYTES,
                                   metadata_blocks=2, new_entries=2)
                self._phase_receipt(receipt_fd, plan)
                self._stop_check()
                cwd = _absolute(cwd)
                environment = self._environment(env)
                with self._directory(cwd), ExitStack() as logs:
                    self._space_sample(audit, walk=True, new_entries=1)
                    out_fd = self._new_file(stdout)
                    logs.callback(os.close, out_fd)
                    self._space_sample(audit, walk=True, new_entries=1)
                    err_fd = self._new_file(stderr)
                    logs.callback(os.close, err_fd)
                    self._space_sample(audit, walk=True, evidence_bytes=_RUNNING_BYTES + _TERMINAL_BYTES,
                                       metadata_blocks=2)
                    audit["prelaunch_space"] = audit["last_space"]
                    self._stop_check()
                    child = _PhaseChild(
                        args, cwd=cwd, env=environment, stdin=subprocess.DEVNULL,
                        stdout=out_fd, stderr=err_fd, start_new_session=True,
                        close_fds=True, restore_signals=True,
                    )
                    leader = _Leader(child)
                if child.pid <= 1 or child.pid in (os.getpgrp(), os.getsid(0)):
                    raise GuardError("group_identity", "Spawn did not yield a distinct positive leader PID")
                self._phase_receipt(receipt_fd, {
                    "status": "running", "label": label, "pid": child.pid,
                    "expected_pgid": child.pid, "expected_sid": child.pid,
                    "uid": self.uid, "utc": datetime.now(timezone.utc).isoformat(),
                    "monotonic": time.monotonic(), "environment_keys": sorted(environment),
                })
                self._monitor(leader, audit)
                if leader.exit not in accepted:
                    raise GuardError("phase_exit", f"Phase returned raw exit {leader.exit}")
            except BaseException as error:
                failure = Failure.from_error(error, "phase_error")
                if primary is None:
                    primary = failure
                elif failure != primary:
                    errors.append(failure)
                errors.extend(getattr(error, "_guard_cleanup_errors", ()))
                if failure.reason == "receipt_write":
                    receipt_error = failure
                audit["trigger_sample"] = {"group": audit["last_group"], "space": audit["last_space"]}

            if leader is not None:
                if audit["terminal_members"] is None:
                    if leader.pid > 1 and leader.pid not in (os.getpgrp(), os.getsid(0)):
                        self._clean_group(leader, audit, errors)
                    else:
                        errors.append(Failure("group_identity", "Unsafe leader PID for cleanup"))
                closure = "unresolved"
                if audit["terminal_members"] is not None:
                    try:
                        leader.reap()
                        closure = "confirmed"
                    except (GuardError, OSError) as error:
                        errors.append(Failure.from_error(error, "leader_reap"))
            try:
                self._space_sample(audit, walk=True)
                if leader is not None:
                    self._stop_check()
            except (GuardError, OSError) as error:
                failure = Failure.from_error(error, "closure_sample")
                if primary is None:
                    primary = failure
                    audit["trigger_sample"] = {"group": audit["last_group"], "space": audit["last_space"]}
                elif failure.reason != primary.reason:
                    errors.append(failure)
            raw_exit = None if leader is None else leader.exit
            status = ("not_started" if leader is None else
                      "error" if primary is not None or errors or closure != "confirmed" else "completed")
            terminal_error = self._finish_receipt(receipt_fd, {
                "status": status, "label": receipt_label, "closure": closure, "raw_exit": raw_exit,
                "leader_reaped": False if leader is None else leader.reaped,
                "primary": None if primary is None else asdict(primary),
                "cleanup_errors": [asdict(error) for error in errors],
                "utc": datetime.now(timezone.utc).isoformat(),
                "elapsed_seconds": time.monotonic() - started, **audit,
            }, audit)
            if receipt_error is None:
                receipt_error = terminal_error
            elif terminal_error is not None:
                errors.append(terminal_error)
        finally:
            try:
                os.close(receipt_fd)
            except OSError as error:
                close_error = Failure.from_error(error, "receipt_close")
                if receipt_error is None:
                    receipt_error = close_error
                else:
                    errors.append(close_error)
        observer_cleanup_failed = (primary is not None and primary.reason == "ps_helper_cleanup"
                                   or any(error.reason == "ps_helper_cleanup" for error in errors))
        self._blocked |= closure == "unresolved" or receipt_error is not None or observer_cleanup_failed
        if status == "completed" and receipt_error is not None:
            status = "error"
        return PhaseOutcome(status, closure, raw_exit, primary, tuple(errors), receipt_error)
