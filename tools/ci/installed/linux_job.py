"""Existing accepted Linux pidfd/subreaper boundary for portable gold commands.

These ownership/drain algorithms are exact copies of the qualified Linux host
core. No provider, bootstrap, cloud, private allocation or package paths remain.
The Windows lane retains its separate hard Job commitment cap.
"""
import ctypes
from collections import defaultdict
from dataclasses import dataclass
import io
import json
import os
from pathlib import Path
import signal
import subprocess
import time

from windows_job import capacity_record
from process_metrics import child_environment
from verify import write_json

PROC = Path("/proc")
GUARD_LOG_BYTES = 16 * 1024 * 1024
GUARD_LOG_FAILURE_RESERVE = 4096
RESERVE = 8_000_000_000
SAMPLED_RSS_LIMIT = 7 * 1024**3

def require(condition, message):
    if not condition:
        raise ValueError(message)

@dataclass(frozen=True)
class Process:
    pid: int
    parent: int
    start: int
    rss_pages: int

def proc_stat(path):
    try:
        fields = (path / "stat").read_text(errors="surrogateescape").rsplit(")", 1)[1].split()
    except (ProcessLookupError, FileNotFoundError):
        return None
    return Process(int(path.name), int(fields[1]), int(fields[19]), int(fields[21]))

def owned_processes(proc, owner):
    children = defaultdict(list)
    for path in proc.iterdir():
        if not path.name.isdecimal():
            continue
        process = proc_stat(path)
        if process:
            children[process.parent].append(process)
    owned = {}
    pending = [owner]
    while pending:
        for process in children[pending.pop()]:
            if process.pid not in owned:
                owned[process.pid] = process
                pending.append(process.pid)
    return owned

def enable_subreaper():
    require(hasattr(os, "pidfd_open") and hasattr(signal, "pidfd_send_signal"),
            "Require the Linux pidfd process APIs")
    libc = ctypes.CDLL(None, use_errno=True)
    libc.prctl.argtypes = [ctypes.c_int, *([ctypes.c_ulong] * 4)]
    if libc.prctl(36, 1, 0, 0, 0) != 0:  # PR_SET_CHILD_SUBREAPER
        number = ctypes.get_errno()
        raise OSError(number, os.strerror(number))
    fd = os.pidfd_open(os.getpid())
    try:
        signal.pidfd_send_signal(fd, 0)
    finally:
        os.close(fd)

def signal_owned(process, signum):
    try:
        fd = os.pidfd_open(process.pid)
    except ProcessLookupError:
        return False
    try:
        current = proc_stat(PROC / str(process.pid))
        if not current or current.start != process.start:
            return False
        signal.pidfd_send_signal(fd, signum)
        return True
    except ProcessLookupError:
        return False
    finally:
        os.close(fd)

def quiesce(child, log):
    """Drain all adopted children after stopping any surviving owned descendants."""
    began = time.monotonic()
    signaled = False
    while True:
        if child is None or child.poll() is not None:
            try:
                while os.waitpid(-1, os.WNOHANG)[0]:
                    pass
            except ChildProcessError:  # ECHILD proves the subreaper has nothing left.
                guard_event(log, "setup_cleanup", owned_tree_quiescent=True, cleanup_signaled=signaled)
                return signaled
        elapsed = time.monotonic() - began
        processes = owned_processes(PROC, os.getpid())
        if elapsed >= 8:
            guard_event(log, "setup_cleanup", owned_tree_quiescent=False,
                        remaining_pids=sorted(processes)[:64], remaining_count=len(processes))
            raise RuntimeError(f"Setup descendants failed to quiesce: {sorted(processes)}")
        signum = signal.SIGTERM if elapsed < 2 else signal.SIGKILL
        for process in processes.values():
            if signal_owned(process, signum):
                signaled = True
        time.sleep(0.1)

def guard_event(log, event, **fields):
    line = json.dumps(dict(epoch=time.time(), event=event, **fields)) + "\n"
    limit = GUARD_LOG_BYTES - (GUARD_LOG_FAILURE_RESERVE if event == "setup_sample" else 0)
    require(log.tell() + len(line) <= limit, "Setup guard journal byte limit reached")
    log.write(line)


class Execution:
    """Serial owned gold commands; RSS is sampled, unlike Windows' hard Job cap."""
    def __init__(self, root, evidence, clock):
        require(os.uname().sysname == "Linux", "Native Linux required")
        enable_subreaper()
        self.root, self.evidence, self.clock = root, evidence, clock
        self.logs = evidence / "logs"
        self.logs.mkdir()
        self.counter = 0
        self.source_pins = None
        (root / "tmp").mkdir()
        (root / "home").mkdir()
        self.environment = child_environment(root / "tmp")
        empty = evidence / "empty-startup"
        empty.write_bytes(b"")
        self.environment.update(HOME=str(root / "home"), R_USER=str(root / "home"),
            R_ENVIRON_USER=str(empty), R_PROFILE_USER=str(empty), R_MAKEVARS_USER=str(empty),
            R_LIBS_SITE="", MAKEFLAGS="-j1", NOT_CRAN="true", TZ="UTC")
        require(not owned_processes(PROC, os.getpid()), "Gold controller has unexpected existing children")

    def capacity(self, phase, growth=0):
        self.clock.remaining("compute")
        memory = {}
        for line in Path("/proc/meminfo").read_text().splitlines():
            key, value = line.split(":", 1)
            if key in {"MemTotal", "MemAvailable"}:
                number, unit = value.strip().split()
                require(unit == "kB", "Unexpected Linux memory units")
                memory[key] = int(number) * 1024
        import shutil
        result = capacity_record(memory["MemTotal"], memory["MemAvailable"],
                                 shutil.disk_usage(self.root).free, growth, phase)
        result["aggregate_RSS_limit_mode"] = "sampled_7GiB_not_hard_Windows_Job_commitment"
        self.counter += 1
        write_json(self.logs / ("capacity-" + str(self.counter).zfill(4) + ".json"), result)
        return result

    def run(self, label, argv, timeout=1800, cwd=None, *, environment=None, provenance=None, sample=None):
        import shutil
        self.capacity(label)
        require(not owned_processes(PROC, os.getpid()), "Unclosed preceding command")
        self.counter += 1
        stem = self.logs / (str(self.counter).zfill(4) + "-" + label)
        command = [str(value) for value in argv]
        began = time.monotonic()
        child = None
        error = None
        peak = 0
        cleanup = io.StringIO()
        rc = None
        signaled = None
        with stem.with_suffix(".stdout").open("xb") as out, stem.with_suffix(".stderr").open("xb") as err:
            try:
                child = subprocess.Popen(command, stdin=subprocess.DEVNULL, stdout=out, stderr=err,
                    cwd=cwd or self.root, env=self.environment if environment is None else environment,
                    close_fds=True, start_new_session=True)
                while child.poll() is None:
                    self.clock.remaining("compute")
                    require(time.monotonic() < began + timeout, "Gold command timeout")
                    require(shutil.disk_usage(self.root).free >= RESERVE, "Gold disk reserve")
                    rss = sum(max(0, item.rss_pages) for item in owned_processes(PROC, os.getpid()).values()) * os.sysconf("SC_PAGE_SIZE")
                    peak = max(peak, rss)
                    require(rss <= SAMPLED_RSS_LIMIT, "Gold aggregate sampled RSS cap")
                    if sample is not None:
                        sample()
                    time.sleep(.1)
                rc = child.returncode
            except BaseException as caught:
                error = caught
            finally:
                try:
                    signaled = quiesce(child, cleanup)
                except BaseException as caught:
                    if error is None:
                        error = caught
                write_json(stem.with_suffix(".json"), dict(argv=command, pid=child.pid if child else None,
                    exit_code=rc, elapsed_seconds=time.monotonic()-began, phase=label,
                    source_pins=self.source_pins, provenance=provenance, maximum_sampled_group_rss=peak,
                    rss_limit_mode="sampled_not_hard_Job", cleanup_signaled=signaled,
                    cleanup_events=[json.loads(line) for line in cleanup.getvalue().splitlines()],
                    error=None if error is None else dict(type=type(error).__name__, message=str(error))))
        if error is not None:
            raise error
        require(rc == 0 and signaled is False and not owned_processes(PROC, os.getpid()),
                "Gold command failed or required cleanup signaling")
        return stem.with_suffix(".stdout")
