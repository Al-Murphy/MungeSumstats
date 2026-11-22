"""Owned native Windows execution/physical evidence primitives.

Extracted from the independently qualified producer. Single-current-package
orchestration is supplied by controller.py; there are no repository/branch,
historical package archive or paired scientific role bindings here.
"""
from contextlib import contextmanager, ExitStack
from datetime import datetime, timezone
import ctypes
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import stat
import subprocess
import sys
import tempfile
import time

PACKAGE_ROOT = Path(__file__).resolve().parents[3]
BENCH = PACKAGE_ROOT / "benchmarks/2026-11-19/performance-upgrade"
sys.path.insert(0, str(BENCH))
from process_metrics import child_environment
from verify import digest, read_json, write_json

RESERVE = 8_000_000_000
ARTIFACT_LIMIT = 512_000_000
JOB_MEMORY = 7 * 1024**3
COMMIT_RESERVE = 512 * 1024**2

JOB_CONTROL_WRITER = '''import json,os,pathlib,sys,time
root=pathlib.Path(sys.argv[1])
started=time.monotonic()
(root/"writer.json").write_text(json.dumps(dict(pid=os.getpid(),parent=os.getppid(),
 started_monotonic=started,natural_exit_monotonic=started+16,attempt=json.loads(sys.argv[2]))))
with (root/"progress.bin").open("xb",buffering=0) as stream:
 while time.monotonic()<started+16:
  stream.write(b"x")
  time.sleep(.05)
'''


class Deadlines:
    """Absolute job-start budgets, also bounded by monotonic elapsed time."""
    def __init__(self, started, now=None, qpc=None, observed=None, observed_tick=None):
        now=time.time if now is None else now
        qpc=time.perf_counter if qpc is None else qpc  # Shared Stopwatch/raw QPC epoch.
        current, tick = now(), qpc()
        if not isinstance(started,(int,float)) or not 0 <= current-started < 360*60:
            raise ValueError("Missing, future or expired GitHub job-start observation")
        self.now,self.qpc = now,qpc
        self.wall = {name:started+minutes*60 for name,minutes in
                     (("compute",300),("drain",305),("seal",316),("final",317))}
        if observed is None:
            observed,observed_tick=current,tick
        if not started<=observed<=current or not 0<=observed_tick<=tick:
            raise ValueError("Invalid first-step monotonic job observation")
        self.ticks = {name:observed_tick+end-observed for name,end in self.wall.items()}

    def remaining(self, phase):
        remaining=min(self.wall[phase]-self.now(),self.ticks[phase]-self.qpc())
        if remaining<=0:
            raise TimeoutError("Absolute GitHub job budget expired: "+phase)
        return remaining


class JobControlCutoff(TimeoutError):
    def __init__(self, phase):
        self.phase=phase
        super().__init__("Owned Job control fixture cutoff: "+phase)


class JobControlClock:
    """Short fixture caps measured once, intersected with the existing job clock."""
    def __init__(self, parent, started):
        self.parent,self.started = parent,started
        self.ends={name:started+seconds for name,seconds in
                   (("readiness",5),("compute",6),("drain",10),("outer",12))}

    def remaining(self, phase):
        parent_phase="compute" if phase=="readiness" else "drain" if phase=="outer" else phase
        original=self.parent.remaining(parent_phase)
        left=self.ends[phase]-time.monotonic()
        if left<=0:
            raise JobControlCutoff(phase)
        return min(original,left)


def receipt_epoch(value):
    # PowerShell's round-trip ISO format carries seven fractional digits. Convert
    # the whole offset-aware second and fraction directly to the clock's epoch.
    match = re.fullmatch(r"(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2})(?:\.(\d{1,7}))?(Z|[+-]\d{2}:\d{2})", value)
    if not match:
        raise ValueError("Job timestamps require an explicit ISO UTC offset")
    second, fraction, offset = match.groups()
    instant = datetime.fromisoformat(second + ("+00:00" if offset == "Z" else offset))
    return instant.timestamp() + int((fraction or "0").ljust(7, "0")) / 10_000_000


def job_clock(path):
    record=json.loads(small_text(path))
    if not all(str(record.get(key,"")).isdecimal() for key in ("job_id","run_id","run_attempt")):
        raise ValueError("Require the actual GitHub Actions job metadata")
    return Deadlines(receipt_epoch(record["started_utc"]), observed=receipt_epoch(record["observed_utc"]),
                     observed_tick=record["observed_monotonic_seconds"])


def windows_host():
    # Runtime checks here are independent of the first-step PowerShell observation.
    native = ctypes.c_ushort()
    process = ctypes.c_ushort()
    kernel = ctypes.WinDLL("kernel32",use_last_error=True)
    kernel.GetCurrentProcess.restype = ctypes.c_void_p
    fn = kernel.IsWow64Process2
    fn.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_ushort), ctypes.POINTER(ctypes.c_ushort)]
    if not fn(kernel.GetCurrentProcess(),ctypes.byref(process),ctypes.byref(native)):
        raise ctypes.WinError(ctypes.get_last_error())
    if native.value != 0x8664 or process.value != 0 or ctypes.sizeof(ctypes.c_void_p) != 8:
        raise ValueError("Python must be a native x64 Windows process")
    return dict(name=platform.node().upper(), os=platform.system(), build=str(sys.getwindowsversion().platform_version[2]),
                os_arch="X64", process_arch="X64")


class AvailableRAMShortfall(ValueError):
    """Only the six-reference admission may wait for available RAM."""


def capacity_record(total_ram, available_ram, free, growth, phase):
    if any(type(x) is not int or x < 0 for x in (total_ram, available_ram, free, growth)):
        raise ValueError("Capacity observations must be nonnegative integer bytes")
    hard_shortfall = total_ram < 7_500_000_000 or free < RESERVE + growth
    if hard_shortfall or available_ram < 5_500_000_000:
        error = AvailableRAMShortfall if phase == "six-reference-minimum" and not hard_shortfall else ValueError
        raise error(f"Insufficient native Windows capacity before {phase}: "
                         f"RAM={total_ram}, available={available_ram}, disk={free}, growth={growth}")
    return dict(phase=phase, total_ram_bytes=total_ram, available_ram_bytes=available_ram,
                free_bytes=free, budget_bytes=growth, reserve_bytes=RESERVE, free_after_budget=free-growth)


def physical_type(info,directory):
    match=stat.S_ISDIR if directory else stat.S_ISREG
    if (not match(info.st_mode) or getattr(info,"st_file_attributes",0)&stat.FILE_ATTRIBUTE_REPARSE_POINT
            or not directory and info.st_nlink!=1):
        raise ValueError("Indirect or unsupported evidence object")


def physical_info(path,parent_fd):
    if os.name=="nt" or parent_fd is None:
        return path.lstat()
    return os.stat(path.name,dir_fd=parent_fd,follow_symlinks=False)


def windows_physical_handle(path,directory,live=False):
    """Open the reparse object itself. Deny directory writes and renames while descending."""
    from ctypes import wintypes as w
    kernel=ctypes.WinDLL("kernel32",use_last_error=True)
    create=kernel.CreateFileW
    create.argtypes=[w.LPCWSTR,w.DWORD,w.DWORD,w.LPVOID,w.DWORD,w.DWORD,w.HANDLE]
    create.restype=w.HANDLE
    close=kernel.CloseHandle
    close.argtypes=[w.HANDLE];close.restype=w.BOOL
    class Info(ctypes.Structure):
        _fields_=[("attributes",w.DWORD),("created",w.FILETIME),("accessed",w.FILETIME),("modified",w.FILETIME),
                  ("volume",w.DWORD),("size_high",w.DWORD),("size_low",w.DWORD),("links",w.DWORD),
                  ("index_high",w.DWORD),("index_low",w.DWORD)]
    get=kernel.GetFileInformationByHandle
    get.argtypes=[w.HANDLE,ctypes.POINTER(Info)];get.restype=w.BOOL
    # FILE_READ_ATTRIBUTES for directories, GENERIC_READ for files. Live raw
    # log files may still have a writer. Neither kind grants share-delete.
    handle=create(str(path),0x80 if directory else 0x80000000,3 if live else 1,None,3,
                  0x00200000|0x02000000,None)  # OPEN_REPARSE_POINT | BACKUP_SEMANTICS
    if handle==ctypes.c_void_p(-1).value:
        raise ctypes.WinError(ctypes.get_last_error())
    try:
        info=Info()
        if not get(handle,ctypes.byref(info)):
            raise ctypes.WinError(ctypes.get_last_error())
        if (info.attributes&stat.FILE_ATTRIBUTE_REPARSE_POINT or bool(info.attributes&0x10)!=directory
                or not directory and info.links!=1):
            raise ValueError("Indirect Windows evidence handle")
        return handle,close
    except BaseException:
        close(handle)
        raise


@contextmanager
def physical_handle(path,parent_fd,directory,live=False):
    before=physical_info(path,parent_fd)
    physical_type(before,directory)
    if os.name=="nt":
        import msvcrt
        handle,close=windows_physical_handle(path,directory,live)
        if directory:
            try:yield None
            finally:close(handle)
            return
        try:fd=msvcrt.open_osfhandle(handle,os.O_RDONLY|os.O_BINARY)
        except BaseException:
            close(handle)
            raise
    else:
        flags=os.O_RDONLY|os.O_NOFOLLOW|(os.O_DIRECTORY if directory else os.O_NONBLOCK)
        fd=os.open(path if parent_fd is None else path.name,flags,dir_fd=parent_fd)
    try:
        after=os.fstat(fd)
        physical_type(after,directory)
        if (before.st_dev,before.st_ino)!=(after.st_dev,after.st_ino):
            raise ValueError("Evidence entry changed during open")
        yield fd
    finally:os.close(fd)


@contextmanager
def physical_directory(path):
    path=path.absolute()
    if ".." in path.parts:raise ValueError("Indirect directory expression")
    # Hold every ancestor until the last read. POSIX children are addressed
    # relative to their open directory. Windows denies rename/reparse on the
    # verified ancestors while enumerating names and opening child handles.
    with ExitStack() as stack:
        parent=Path(path.anchor)
        fd=None if os.name=="nt" else stack.enter_context(physical_handle(parent,None,True))
        for part in path.parts[1:]:
            parent=parent/part
            fd=stack.enter_context(physical_handle(parent,fd,True))
        yield fd


def transfer_regular(path,parent_fd,limit,destination=None,live=False):
    with physical_handle(path,parent_fd,False,live) as fd,ExitStack() as stack:
        before=os.fstat(fd)
        if before.st_size>limit and not live:
            raise RuntimeError("Evidence file exceeds its byte budget")
        wanted=min(before.st_size,limit)
        out=None
        if destination is not None:
            destination.parent.mkdir(parents=True,exist_ok=True)
            out=stack.enter_context(destination.open("xb"))
        sha,md5=hashlib.sha256(),hashlib.md5(usedforsecurity=False)
        remaining=wanted
        while remaining:
            block=os.read(fd,min(remaining,1<<20))
            if not block:raise ValueError("Evidence changed during read")
            if out is not None:out.write(block)
            sha.update(block);md5.update(block);remaining-=len(block)
        after=os.fstat(fd)
        physical_type(after,False)
        if not live and (before.st_size,before.st_mtime_ns)!=(after.st_size,after.st_mtime_ns):
            raise ValueError("Evidence changed during read")
        return dict(bytes=wanted,sha256=sha.hexdigest(),md5=md5.hexdigest()),max(before.st_size,after.st_size)>wanted


def evidence_inventory(evidence,clock=None,destination=None,metadata_only=False):
    from evidence import MAX_MEMBERS,relative_name
    files,directories,total={},[],0
    count=0
    def walk(directory,parent_fd):
        nonlocal total,count
        with os.scandir(directory if os.name=="nt" else parent_fd) as entries:
            for entry in entries:
                count+=1
                if count>MAX_MEMBERS:raise ValueError("Too many evidence paths")
                if clock is not None:clock.remaining("seal")
                path=directory/entry.name
                name=relative_name(path.relative_to(evidence).as_posix())
                info=physical_info(path,parent_fd)
                isdir=stat.S_ISDIR(info.st_mode)
                physical_type(info,isdir)
                target=destination/name if destination is not None else None
                if isdir:
                    with physical_handle(path,parent_fd,True) as child:
                        directories.append(name)
                        if target is not None:target.mkdir()
                        walk(path,child)
                elif not metadata_only:
                    item,_=transfer_regular(path,parent_fd,ARTIFACT_LIMIT-total,target)
                    files[name]=item;total+=item["bytes"]
    with physical_directory(evidence) as fd:
        walk(evidence,fd)
    return dict(bytes=total,files=files,directories=sorted(directories))


def seal_evidence(evidence,clock=None):
    return evidence_inventory(evidence,clock)


def prepare_upload_snapshot(attempt,area):
    """Return an ordinary isolated upload directory; live failure paths are never enumerated."""
    from evidence import OUTSIDE_LIMITS,SUCCESS_OUTSIDE
    attempt,area=attempt.absolute(),area.absolute()
    if area==attempt or attempt in area.parents:raise ValueError("Snapshot must be outside attempt")
    with physical_directory(area):
        private=Path(tempfile.mkdtemp(prefix="mss-upload-",dir=area))
    payload=private/"payload";payload.mkdir()
    pending=private/"pending";pending.mkdir()
    omitted={}
    sequence=0
    def take(name,limit,live=False):
        nonlocal sequence
        sequence+=1
        target=payload/name
        incoming=pending/str(sequence)
        path=attempt/name
        try:
            with physical_directory(path.parent) as fd:
                _,truncated=transfer_regular(path,fd,limit,incoming,live)
            target.parent.mkdir(parents=True,exist_ok=True)
            incoming.rename(target)
            if truncated:omitted[name]="only the bounded observed prefix was copied"
        except FileNotFoundError:omitted[name]="absent"
        except (OSError,ValueError,RuntimeError) as error:omitted[name]="unreadable or nonphysical: "+type(error).__name__
    for name,limit in OUTSIDE_LIMITS.items():
        if name!="q/upload-omissions.json":
            take(name,limit,name.endswith((".stdout",".stderr")))
    try:
        status=read_json(payload/"q/artifact-status.json")
        if ((payload/"q/UNSEALED.json").exists() or status["source_sealed"] is not True
                or omitted.get("q/UNSEALED.json")!="absent"
                or any(read_json(payload/f"q/{mode}-stop.json")["drained"] is not True for mode in ("job-control","qualify","seal"))):
            raise ValueError("missing source seal or confirmed drains")
        if status["qualified"] is True and (any(name in omitted for name in SUCCESS_OUTSIDE)
                or any((payload/name).exists() for name in set(OUTSIDE_LIMITS)-SUCCESS_OUTSIDE)):
            raise ValueError("incomplete required qualification receipts")
        inventory=payload/"q/artifact-inventory.json"
        if digest(inventory)["sha256"]!=status["inventory"]["sha256"]:
            raise ValueError("inconsistent source inventory pin")
        expected=read_json(inventory)
        if not 0<=expected["bytes"]<=ARTIFACT_LIMIT or shutil.disk_usage(area).free<RESERVE+expected["bytes"]:
            raise ValueError("insufficient bounded snapshot capacity")
        candidate=private/"complete-evidence";candidate.mkdir()
        actual=evidence_inventory(attempt/"q/evidence",destination=candidate)
        if actual!=expected:raise ValueError("snapshot differs from independent source seal")
        candidate.rename(payload/"q/evidence")
        return payload,True
    except (KeyError,TypeError,ValueError,OSError,RuntimeError) as error:
        omitted["q/evidence"]="full tree withheld: "+type(error).__name__
    # Only these named records may be read after a seal failure or unknown
    # drain. Pin ancestors and the file first. A growing log contributes only
    # the prefix measured on the opened file, capped by its diagnostic limit.
    diagnostics={name:1<<20 for name in ("qualification.json","qualification-candidate.json","failure.json",
        "job-control-stop.json","qualify-stop.json","scratch-cleanup.json","job-control/controller.json","job-control/writer.json")}
    diagnostics["job-control/progress.bin"]=4096
    for stem in ("worker","job-control"):
        for suffix in ("stdout","stderr"):diagnostics[stem+"."+suffix]=8<<20
    for name,limit in diagnostics.items():
        take("q/evidence/"+name,limit,name.endswith((".stdout",".stderr",".bin")))
    (payload/"q").mkdir(exist_ok=True)
    write_json(payload/"q/upload-omissions.json",dict(full_tree_uploaded=False,omitted=omitted))
    return payload,False


def memory_status():
    from ctypes import wintypes as w
    class Status(ctypes.Structure):
        _fields_ = [("length",w.DWORD),("load",w.DWORD)] + [(name,ctypes.c_ulonglong) for name in
            ("total","available","total_page","available_page","total_virtual","available_virtual","extended")]
    status = Status()
    status.length = ctypes.sizeof(status)
    function = ctypes.WinDLL("kernel32",use_last_error=True).GlobalMemoryStatusEx
    function.argtypes = [ctypes.POINTER(Status)]
    function.restype = w.BOOL
    if not function(ctypes.byref(status)):
        raise ctypes.WinError(ctypes.get_last_error())
    return int(status.total), int(status.available)


def system_commit_status():
    """Native system commitment, including the pagefile; values are bytes."""
    from ctypes import wintypes as wt
    class Performance(ctypes.Structure):
        _fields_ = [("cb",wt.DWORD)] + [(n,ctypes.c_size_t) for n in (
            "commit","limit","peak","physical","available","cache","kernel","paged",
            "nonpaged","page_size")] + [(n,wt.DWORD) for n in ("handles","processes","threads")]
    value = Performance()
    value.cb = ctypes.sizeof(value)
    kernel = ctypes.WinDLL("kernel32",use_last_error=True)
    observe = kernel.K32GetPerformanceInfo
    observe.argtypes = [ctypes.POINTER(Performance),wt.DWORD]
    observe.restype = wt.BOOL
    if not observe(ctypes.byref(value),value.cb):
        raise ctypes.WinError(ctypes.get_last_error())
    return int(value.limit*value.page_size),int(value.commit*value.page_size)


def commitment_record(limit,used,before,after,active,phase):
    """Reserve the complete Job liability against native global commitment."""
    if (any(type(x) is not int or x<0 for x in (limit,used,before,after,active))
            or not 0<used<=limit or max(before,after)>used
            or phase not in ("before-worker","before-release")
            or (phase=="before-worker" and (active!=0 or before!=0 or after!=0))
            or (phase=="before-release" and (active!=1 or min(before,after)<=0))):
        raise ValueError("Invalid native committed-memory observation")
    current=min(before,after)  # Conservative bracket; never subtract RSS.
    outside=used-current
    remaining=limit-outside-JOB_MEMORY-COMMIT_RESERVE
    return dict(phase=phase,system_commit_limit_bytes=limit,system_commit_used_bytes=used,
        job_current_before_bytes=before,job_current_after_bytes=after,job_current_bytes=current,
        active_members=active,outside_job_commit_bytes=outside,job_limit_bytes=JOB_MEMORY,
        global_commit_reserve_bytes=COMMIT_RESERVE,remaining_after_reserve_bytes=remaining,
        decision="admitted" if remaining>=0 else "rejected")


def admit_job_commitment(root,job,phase,clock):
    clock.remaining("compute")
    _,anchor=phase_anchor(clock)
    before=job.committed_bytes()
    limit,used=system_commit_status()
    after=job.committed_bytes()
    record=commitment_record(limit,used,before,after,job.active_count(),phase)
    record.update(anchor=anchor,observer_pid=os.getpid())
    clock.remaining("compute")
    with (root/"evidence/job-commit-admission.jsonl").open("x" if phase=="before-worker" else "a",encoding="ascii") as stream:
        stream.write(json.dumps(record,sort_keys=True)+"\n")
        stream.flush()
        os.fsync(stream.fileno())
    clock.remaining("compute")
    if record["decision"]!="admitted":
        raise ValueError("Insufficient global commitment for bounded qualification Job")
    return record


class OwnedJob:
    """The outside supervisor owns the handle; a gated worker joins before spawning."""
    def __init__(self):
        from ctypes import wintypes as w
        class BasicLimits(ctypes.Structure):
            _fields_ = [("process_time",ctypes.c_longlong),("job_time",ctypes.c_longlong),
                ("flags",w.DWORD),("minimum",ctypes.c_size_t),("maximum",ctypes.c_size_t),
                ("active",w.DWORD),("affinity",ctypes.c_size_t),("priority",w.DWORD),("scheduling",w.DWORD)]
        class IO(ctypes.Structure):
            _fields_ = [(x,ctypes.c_ulonglong) for x in ("read_ops","write_ops","other_ops","read","write","other")]
        class Limits(ctypes.Structure):
            _fields_ = [("basic",BasicLimits),("io",IO),("process_memory",ctypes.c_size_t),
                ("job_memory",ctypes.c_size_t),("peak_process",ctypes.c_size_t),("peak_job",ctypes.c_size_t)]
        class Accounting(ctypes.Structure):
            _fields_ = [(x,ctypes.c_longlong) for x in ("user","kernel","period_user","period_kernel")] + [
                (x,w.DWORD) for x in ("faults","total","active","terminated")]
        self.Accounting = Accounting
        self.kernel = k = ctypes.WinDLL("kernel32",use_last_error=True)
        k.CreateJobObjectW.argtypes, k.CreateJobObjectW.restype = [ctypes.c_void_p,w.LPCWSTR],w.HANDLE
        k.OpenProcess.argtypes,k.OpenProcess.restype = [w.DWORD,w.BOOL,w.DWORD],w.HANDLE
        k.CloseHandle.argtypes,k.CloseHandle.restype = [w.HANDLE],w.BOOL
        k.TerminateJobObject.argtypes,k.TerminateJobObject.restype = [w.HANDLE,w.UINT],w.BOOL
        k.SetInformationJobObject.argtypes = [w.HANDLE,ctypes.c_int,ctypes.c_void_p,w.DWORD]
        k.SetInformationJobObject.restype = w.BOOL
        k.AssignProcessToJobObject.argtypes, k.AssignProcessToJobObject.restype = [w.HANDLE,w.HANDLE],w.BOOL
        k.QueryInformationJobObject.argtypes = [w.HANDLE,ctypes.c_int,ctypes.c_void_p,w.DWORD,ctypes.c_void_p]
        k.QueryInformationJobObject.restype = w.BOOL
        self.handle = k.CreateJobObjectW(None,None)
        limits = Limits()
        limits.basic.flags = 0x2000 | 0x0200  # KILL_ON_JOB_CLOSE, JOB_MEMORY; no breakaway.
        limits.job_memory = JOB_MEMORY
        if not self.handle or not k.SetInformationJobObject(self.handle,9,ctypes.byref(limits),ctypes.sizeof(limits)):
            raise ctypes.WinError(ctypes.get_last_error())

    def assign(self, pid):
        process=self.kernel.OpenProcess(0x0100|0x0001,False,pid)
        if not process:
            raise ctypes.WinError(ctypes.get_last_error())
        try:
            if not self.kernel.AssignProcessToJobObject(self.handle,process):
                raise ctypes.WinError(ctypes.get_last_error())
        finally:
            self.kernel.CloseHandle(process)

    def committed_bytes(self):
        fields = (ctypes.c_uint64 * 2)()
        if not self.kernel.QueryInformationJobObject(self.handle,28,ctypes.byref(fields),ctypes.sizeof(fields),None):
            raise ctypes.WinError(ctypes.get_last_error())
        return int(fields[0])

    def active_count(self):
        value = self.Accounting()
        if not self.kernel.QueryInformationJobObject(self.handle,1,ctypes.byref(value),ctypes.sizeof(value),None):
            raise ctypes.WinError(ctypes.get_last_error())
        return int(value.active)

    def terminate(self):
        if not self.kernel.TerminateJobObject(self.handle,99):
            raise ctypes.WinError(ctypes.get_last_error())

    def close(self):
        if not self.kernel.CloseHandle(self.handle):
            raise ctypes.WinError(ctypes.get_last_error())


def terminate_and_drain(job, process, clock, phase):
    try:
        before=job.active_count()
    finally:
        job.terminate()
    until=time.monotonic()+min(60,clock.remaining(phase))
    while job.active_count()!=0:
        if time.monotonic()>=until:
            raise TimeoutError("Owned Windows Job did not drain; do not seal")
        clock.remaining(phase)
        time.sleep(.02)
    rc=process.wait(timeout=min(2,clock.remaining(phase)))
    clock.remaining(phase)
    return dict(drained=True,active_before=before,active_after=0,exit_code=rc)


class CompletedCommandError(RuntimeError):
    """The command exited by itself with a nonzero status."""


class Execution:
    def __init__(self, root, evidence, clock):
        self.root, self.evidence = root, evidence
        self.logs = evidence/"logs"
        self.logs.mkdir()
        self.clock = clock
        self.environment = child_environment(root/"tmp")
        (root/"tmp").mkdir()
        empty = evidence/"empty-startup"
        empty.write_bytes(b"")
        self.environment.update(R_ENVIRON_USER=str(empty),R_PROFILE_USER=str(empty),R_MAKEVARS_USER=str(empty),
            R_LIBS_SITE="",MAKEFLAGS="-j1",NOT_CRAN="true",USE_BUNDLED_LIBUV="1",TZ="UTC")
        home = root/"home"
        home.mkdir()
        self.environment.update(HOME=str(home), R_USER=str(home))
        self.counter = 0
        self.source_pins = None

    def capacity_sample(self, phase, growth, observation):
        self.clock.remaining("compute")
        total, available = memory_status()
        # All writes share the owned root's volume, including R's temporary build tree.
        volumes = {}
        observation.update(total_ram_bytes=total, available_ram_bytes=available, volumes=volumes)
        for subdir in ("tmp", "dependency-archives", "dependencies"):
            path = self.root/subdir
            observed = path if path.exists() else self.root
            if observed.resolve(strict=True) != observed or observed.stat().st_dev != self.root.stat().st_dev:
                raise ValueError("Expected a direct path on the owned work volume: " + str(path))
            volumes[subdir] = dict(path=str(observed),free_bytes=shutil.disk_usage(observed).free)
        free = min(item["free_bytes"] for item in volumes.values())
        value = capacity_record(total,available,free,growth,phase)
        value["volumes"] = volumes
        return value

    def reference_capacity(self, phase, growth):
        started = time.monotonic()

        def remaining():
            original = self.clock.remaining("compute")
            left = started + 120 - time.monotonic()
            if left <= 0:
                raise TimeoutError("Six-reference available-RAM admission expired after 120 seconds")
            return min(original, left)

        with (self.logs/"six-reference-admission.jsonl").open("x", encoding="utf-8") as log:
            while True:
                remaining()
                observation = dict(utc=datetime.now(timezone.utc).isoformat(), growth_bytes=growth,
                                   elapsed_seconds=time.monotonic()-started)
                try:
                    value = self.capacity_sample(phase, growth, observation)
                    remaining()  # A read that completed after either deadline cannot admit work.
                except AvailableRAMShortfall as error:
                    observation.update(decision="wait", reason=str(error))
                except Exception as error:
                    observation.update(decision="rejected", reason=str(error))
                    raise
                else:
                    observation["decision"] = "admitted"
                    return value
                finally:
                    log.write(json.dumps(observation) + "\n")
                    log.flush()
                time.sleep(min(5, remaining()))

    def capacity(self, phase, growth=0):
        if phase == "six-reference-minimum":
            value = self.reference_capacity(phase, growth)
        else:
            value = self.capacity_sample(phase, growth, {})
        path = self.logs/f"capacity-{self.counter:04d}.json"
        self.counter += 1
        write_json(path,value)
        return value

    def run(self, label, argv, timeout=1800, cwd=None, *, environment=None, provenance=None, sample=None):
        self.capacity(label)
        self.counter += 1
        stem = self.logs/f"{self.counter:04d}-{label}"
        started = time.monotonic()
        command = [x.as_posix() if isinstance(x,Path) else str(x) for x in argv]
        with stem.with_suffix(".stdout").open("xb") as out, stem.with_suffix(".stderr").open("xb") as err:
            process = subprocess.Popen(command,stdin=subprocess.DEVNULL,stdout=out,stderr=err,
                                       cwd=cwd or self.root,env=self.environment if environment is None else environment,
                                       close_fds=True)
            reason = None
            while process.poll() is None:
                try:
                    remaining=self.clock.remaining("compute")
                except TimeoutError:
                    remaining=0
                if remaining<=0 or time.monotonic()>=started+timeout:
                    reason = "command_or_job_deadline"
                elif shutil.disk_usage(self.root).free < RESERVE:
                    reason = "disk_reserve"
                if reason:
                    process.terminate()
                    break
                if sample is not None:
                    sample()
                time.sleep(.5)
            if sample is not None:
                sample()
            rc = process.poll()
            write_json(stem.with_suffix(".json"),dict(argv=command,pid=process.pid,exit_code=rc,
                       elapsed_seconds=time.monotonic()-started,stopped=reason,phase=label,
                       source_pins=self.source_pins,provenance=provenance))
        if reason or rc != 0:
            error = CompletedCommandError if reason is None and rc is not None else RuntimeError
            raise error(f"{label} failed; see {stem.name}, reason={reason}, exit={rc}")
        return stem.with_suffix(".stdout")


def small_text(path):
    with path.open("rb") as stream:
        data = stream.read(1_048_577)
    if len(data)>1_048_576:
        raise ValueError("Unexpectedly large metadata output")
    return data.decode("utf-8-sig").strip()


def job_control_fixture(root, clock):
    """Called only after the real assignment gate; both natural exits are later than the cutoff."""
    started=time.monotonic()
    directory=root/"evidence/job-control"
    directory.mkdir()
    write_json(directory/"controller.json",dict(pid=os.getpid(),parent=os.getppid(),
               started_monotonic=started,natural_exit_monotonic=started+16,attempt=clock.attempt))
    subprocess.Popen([sys.executable,"-I","-S","-B","-c",JOB_CONTROL_WRITER,str(directory),json.dumps(clock.attempt)],
                     stdin=subprocess.DEVNULL,env=child_environment(root),close_fds=True)
    while time.monotonic()<started+16:
        clock.remaining("compute")
        time.sleep(.05)


def job_control_progress(root, process, job, started, attempt):
    directory=root/"evidence/job-control"
    paths=[directory/name for name in ("controller.json","writer.json","progress.bin")]
    count=job.active_count()
    if count<2 or not all(path.is_file() for path in paths):
        return None
    controller,writer=(json.loads(small_text(path)) for path in paths[:2])
    for record in (controller,writer):
        begin,end=record["started_monotonic"],record["natural_exit_monotonic"]
        if record.get("attempt")!=attempt or not started<=begin<=time.monotonic() or not end>=begin+16:
            raise ValueError("Invalid Job control natural-exit observation")
    if (controller["pid"]!=process.pid or controller["parent"]!=os.getpid()
            or writer["parent"]!=process.pid or writer["pid"]==process.pid
            or writer["started_monotonic"]<controller["started_monotonic"]):
        raise ValueError("Incorrect Job control descendant identity")
    size=paths[-1].stat().st_size
    if size==0:
        return None
    if not 0<size<=4096:
        raise ValueError("Invalid bounded Job control writer progress")
    return dict(elapsed_seconds=time.monotonic()-started,active_members=count,bytes=size,
                controller_pid=controller["pid"],writer_pid=writer["pid"])


def job_control_snapshot(root, clock, started):
    clock.remaining("outer")
    path=root/"evidence/job-control/progress.bin"
    before=path.lstat()
    if not stat.S_ISREG(before.st_mode) or not 0<before.st_size<=4096:
        raise ValueError("Expected the bounded regular Job control writer")
    with path.open("rb") as stream:
        data=stream.read(4097)
    after=path.lstat()
    if (len(data)!=before.st_size or (before.st_size,before.st_mtime_ns,before.st_ino)
            !=(after.st_size,after.st_mtime_ns,after.st_ino)):
        raise ValueError("Job control writer changed while observed after drain")
    clock.remaining("outer")
    return dict(elapsed_seconds=time.monotonic()-started,bytes=len(data),sha256=hashlib.sha256(data).hexdigest())


def phase_anchor(clock):
    # Bracket the local monotonic and wall observations in the shared QPC domain.
    first=clock.qpc()
    local=time.monotonic()
    wall=clock.now()
    return local,dict(local_monotonic=local,utc_epoch=wall,qpc=first,qpc_end=clock.qpc())


def owned_phase(root, clock, mode, start_file, entry):
    """Fixed workers join a noninheritable owned Job before any descendants start."""
    entry = Path(entry).resolve(strict=True)
    if entry != Path(__file__).with_name("controller.py").resolve(strict=True):
        raise ValueError("Owned worker entry must be the current checkout controller")
    if mode not in ("qualify","seal","job-control"):
        raise ValueError("Unknown owned worker mode")
    control=mode=="job-control"
    parent_clock=clock
    started,launch=phase_anchor(parent_clock)
    anchors={"launch":launch}
    if control:clock=JobControlClock(clock,started)
    phase="compute" if mode!="seal" else "seal"
    drain_phase="drain" if mode!="seal" else "final"
    clock.remaining(phase)
    job=OwnedJob()
    process=None
    assigned=False
    reason=None
    stop=None
    observations=[]
    snapshots=[]
    trigger=None
    drained_elapsed=None
    try:
        with (root/f"{mode}.stdout").open("xb") as out,(root/f"{mode}.stderr").open("xb") as err:
            if mode=="qualify":admit_job_commitment(root,job,"before-worker",clock)
            process=subprocess.Popen([sys.executable,"-I","-S","-B",str(entry),str(root),
                "--job-start",str(start_file),"--internal-worker",mode],stdin=subprocess.DEVNULL,
                stdout=out,stderr=err,env=child_environment(root),close_fds=True)
        job.assign(process.pid)
        assigned=True
        if job.active_count()!=1:
            raise RuntimeError("Gated worker did not enter its exact owned Windows Job")
        clock.remaining("readiness" if control else phase)
        if mode=="qualify":admit_job_commitment(root,job,"before-release",clock)
        temporary=root/f"{mode}-go.new"
        write_json(temporary,dict(pid=process.pid,parent=os.getpid(),mode=mode))
        temporary.rename(root/f"{mode}-go.json")
        while process.poll() is None:
            clock.remaining(phase)
            if shutil.disk_usage(root).free<RESERVE:
                reason="disk_reserve"
                break
            if control and len(observations)<2:
                clock.remaining("readiness")
                progress=job_control_progress(root,process,job,started,parent_clock.attempt)
                clock.remaining("readiness")
                if progress is not None and (not observations or progress["bytes"]>observations[0]["bytes"]):
                    observations.append(progress)
            time.sleep(min(.1,clock.remaining(phase)))
        if control and reason is None:
            reason="job_control_natural_exit"
    except JobControlCutoff as error:
        if control and error.phase=="compute":
            reason="job_control_cutoff"
            local,anchors["cutoff"]=phase_anchor(parent_clock)
            trigger=dict(phase=error.phase,elapsed_seconds=local-started)
        else:
            reason=type(error).__name__+": "+str(error)
    except TimeoutError:
        reason="absolute_"+phase+"_cutoff"
    except BaseException as error:
        reason=type(error).__name__+": "+str(error)
    finally:
        try:
            if process is not None:
                if not assigned:
                    # No go marker exists: this worker cannot have spawned children.
                    process.terminate()
                    process.wait(timeout=min(5,clock.remaining(drain_phase)))
                stop=terminate_and_drain(job,process,clock,drain_phase)
                if control and stop["drained"]:
                    local,anchors["drain"]=phase_anchor(parent_clock)
                    drained_elapsed=local-started
            else:
                job.terminate()
                stop=dict(drained=job.active_count()==0,active_before=0,active_after=job.active_count(),exit_code=None)
        except BaseException as error:
            stop=dict(drained=False,error_type=type(error).__name__,drain_error=str(error))
        finally:
            job.close()
    report=dict(mode=mode,pid=process.pid if process else None,assigned=assigned,reason=reason,**stop)
    if control:
        passed=(assigned and reason=="job_control_cutoff" and len(observations)==2
                and stop["drained"] and stop.get("active_before",0)>=2
                and stop.get("active_after")==0 and stop.get("exit_code")==99)
        if passed:
            try:
                snapshots.append(job_control_snapshot(root,clock,started))
                time.sleep(min(.05,clock.remaining("outer")))
                snapshots.append(job_control_snapshot(root,clock,started))
                clock.remaining("outer")
                if any(snapshots[0][key]!=snapshots[1][key] for key in ("bytes","sha256")):
                    raise ValueError("Job control writer changed after confirmed drain")
            except BaseException as error:
                passed=False
                report["post_drain_error"]=type(error).__name__+": "+str(error)
        local,anchors["complete"]=phase_anchor(parent_clock)
        report.update(attempt=parent_clock.attempt,anchors=anchors,control_passed=passed,readiness_observations=observations,
                      post_drain_observations=snapshots,controller_started_monotonic=started,
                      fixture_deadlines_monotonic=clock.ends,elapsed_seconds=local-started,
                      deadline_trigger=trigger,drain_elapsed_seconds=drained_elapsed)
    if mode=="qualify":report.update(launch_anchor=launch,attempt=parent_clock.attempt)
    if not control and report["drained"]:
        report.update(drained_utc_epoch=parent_clock.now(),drained_qpc=parent_clock.qpc())
    write_json(root/f"{mode}-stop.json",report)
    if control and report["control_passed"]:clock.remaining("outer")
    return report


def wait_for_gate(root, clock, mode):
    if mode not in ("qualify","seal","job-control"):
        raise ValueError("Unknown owned worker mode")
    phase="compute" if mode!="seal" else "seal"
    limit=time.monotonic()+min(5 if mode=="job-control" else 15,clock.remaining(phase))
    gate=root/f"{mode}-go.json"
    while not gate.exists():
        clock.remaining(phase)
        if time.monotonic()>=limit:
            raise TimeoutError("Supervisor did not assign the gated worker")
        time.sleep(.02)
    if json.loads(small_text(gate))!=dict(pid=os.getpid(),parent=os.getppid(),mode=mode):
        raise ValueError("Incorrect owned worker assignment handshake")
