"""Accounting for the single fresh R process used by each benchmark run."""

import os
from pathlib import Path
import subprocess
import sys
import time


def r_command(executable, worker, *args):
    if sys.platform == "win32":
        # Select the architecture-specific Rterm.exe, so the measured PID is R.
        if Path(executable).name.lower() != "rterm.exe":
            raise ValueError("Windows requires the architecture-specific Rterm.exe")
        return [str(executable), "--vanilla", "--slave", f"--file={worker}", "--args", *args]
    return [str(executable), "--vanilla", str(worker), *args]


def child_environment(scratch):
    # Paths and locale are inherited where needed. No user R startup, library,
    # credential, package-manager, cloud, or tracing variables are propagated.
    allowed = {
        "PATH", "HOME", "USERPROFILE", "SYSTEMROOT", "WINDIR", "COMSPEC",
        "SYSTEMDRIVE", "LANG", "LANGUAGE", "TZ",
        "LC_ALL", "LC_CTYPE", "LC_COLLATE", "LC_MESSAGES", "LC_MONETARY",
        "LC_NUMERIC", "LC_TIME",
    }
    env = {name: value for name, value in os.environ.items() if name in allowed}
    for name in ("TMPDIR", "TEMP", "TMP"):
        env[name] = str(scratch)
    for name in (
        "R_DATATABLE_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS",
    ):
        env[name] = "1"
    env["OMP_DYNAMIC"] = "FALSE"
    env["MKL_DYNAMIC"] = "FALSE"
    return env


def windows_accounting(child):
    import ctypes
    from ctypes import wintypes

    class Memory(ctypes.Structure):
        _fields_ = [
            ("size", wintypes.DWORD), ("faults", wintypes.DWORD),
            ("peak_working_set", ctypes.c_size_t), ("working_set", ctypes.c_size_t),
            ("peak_paged_pool", ctypes.c_size_t), ("paged_pool", ctypes.c_size_t),
            ("peak_nonpaged_pool", ctypes.c_size_t), ("nonpaged_pool", ctypes.c_size_t),
            ("pagefile", ctypes.c_size_t), ("peak_pagefile", ctypes.c_size_t),
        ]

    kernel = ctypes.WinDLL("kernel32", use_last_error=True)
    psapi = ctypes.WinDLL("psapi", use_last_error=True)
    kernel.OpenProcess.argtypes = (wintypes.DWORD, wintypes.BOOL, wintypes.DWORD)
    kernel.OpenProcess.restype = wintypes.HANDLE
    kernel.CloseHandle.argtypes = (wintypes.HANDLE,)
    kernel.CloseHandle.restype = wintypes.BOOL
    kernel.GetProcessTimes.argtypes = (wintypes.HANDLE,) + (ctypes.POINTER(wintypes.FILETIME),) * 4
    kernel.GetProcessTimes.restype = wintypes.BOOL
    psapi.GetProcessMemoryInfo.argtypes = (wintypes.HANDLE, ctypes.POINTER(Memory), wintypes.DWORD)
    psapi.GetProcessMemoryInfo.restype = wintypes.BOOL
    handle = kernel.OpenProcess(0x0400 | 0x0010, False, child.pid)
    if not handle:
        child.wait()
        raise ctypes.WinError(ctypes.get_last_error())

    def peak():
        usage = Memory()
        usage.size = ctypes.sizeof(usage)
        if psapi.GetProcessMemoryInfo(handle, ctypes.byref(usage), usage.size):
            return int(usage.peak_working_set)
        return 0

    def seconds(filetime):
        ticks = (filetime.dwHighDateTime << 32) | filetime.dwLowDateTime
        return ticks / 10_000_000

    try:
        observed = 0
        while child.poll() is None:
            observed = max(observed, peak())
            time.sleep(0.05)
        last_peak = peak()
        created, exited, system, user = (wintypes.FILETIME() for _ in range(4))
        if not kernel.GetProcessTimes(handle, ctypes.byref(created), ctypes.byref(exited),
                                      ctypes.byref(system), ctypes.byref(user)):
            raise ctypes.WinError(ctypes.get_last_error())
        return {
            "user_seconds": seconds(user),
            "system_seconds": seconds(system),
            "peak_rss_bytes": max(observed, last_peak) or None,
            "rss_source": "Win32 peak working set",
            "rss_complete": bool(last_peak and last_peak >= observed),
        }
    finally:
        kernel.CloseHandle(handle)


def accounting(child):
    if sys.platform == "win32":
        return windows_accounting(child)
    if sys.platform not in {"linux", "darwin"}:
        raise RuntimeError("Process accounting is supported on Linux, macOS and Windows")
    _, status, usage = os.wait4(child.pid, 0)
    child.returncode = os.waitstatus_to_exitcode(status)
    # Darwin reports bytes and Linux reports KiB.
    multiplier = 1 if sys.platform == "darwin" else 1024
    return {
        "user_seconds": usage.ru_utime,
        "system_seconds": usage.ru_stime,
        "peak_rss_bytes": usage.ru_maxrss * multiplier,
        "rss_source": "wait4 ru_maxrss",
        "rss_complete": usage.ru_maxrss > 0,
    }


def run_process(command, directory, *, name="process"):
    scratch = directory / "scratch"
    scratch.mkdir(exist_ok=True)
    with (directory / f"{name}.stdout.log").open("xb") as out, \
         (directory / f"{name}.stderr.log").open("xb") as err:
        started = time.monotonic()
        child = subprocess.Popen(command, cwd=directory, env=child_environment(scratch),
                                 stdin=subprocess.DEVNULL, stdout=out, stderr=err)
        try:
            result = accounting(child)
        except BaseException:
            if child.poll() is None:
                child.terminate()
                child.wait()
            raise
    result.update(pid=child.pid, exit_code=child.returncode, command=command,
                  wall_seconds=time.monotonic() - started)
    return result
