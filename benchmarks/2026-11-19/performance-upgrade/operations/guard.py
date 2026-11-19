"""Linux wall, owned process-group RSS and scratch free-space watchdog."""
import argparse
import datetime
import json
import os
from pathlib import Path
import resource
import shutil
import signal
import subprocess
import time


def rss(group):
    pages = 0
    for path in Path("/proc").iterdir():
        if not path.name.isdecimal():
            continue
        try:
            # stat command name can itself carry whitespace or parentheses.
            fields = (path / "stat").read_text().rsplit(")", 1)[1].split()
            if int(fields[2]) == group:
                pages += int(fields[21])
        except (ProcessLookupError, FileNotFoundError):
            continue
    return pages * os.sysconf("SC_PAGE_SIZE")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--scratch", type=Path, required=True)
    ap.add_argument("--allocation-utc", required=True, help="actual allocation start, ISO8601 with offset")
    ap.add_argument("--log", type=Path, required=True)
    ap.add_argument("argv", nargs=argparse.REMAINDER)
    a = ap.parse_args()
    start = datetime.datetime.fromisoformat(a.allocation_utc.replace("Z", "+00:00"))
    if start.tzinfo is None or not a.argv or a.argv[0] != "--" or len(a.argv) < 2:
        ap.error("Need absolute allocation time and -- COMMAND [ARG ...]")
    remaining = start.timestamp() + 16 * 3600 - time.time()
    if not 0 < remaining <= 16 * 3600:
        ap.error("Outside the admitted compute window")
    deadline = time.monotonic() + remaining
    reserve = 8_000_000_000
    limit = 24_000_000_000
    def per_process():
        _, hard = resource.getrlimit(resource.RLIMIT_AS)
        cap = min(26_000_000_000, hard) if hard != resource.RLIM_INFINITY else 26_000_000_000
        resource.setrlimit(resource.RLIMIT_AS, (cap, cap))
    def interrupted(signum, frame):
        raise KeyboardInterrupt
    for sig in (signal.SIGHUP, signal.SIGTERM):
        signal.signal(sig, interrupted)
    with a.log.open("x", buffering=1) as log:
        if shutil.disk_usage(a.scratch).free <= reserve:
            raise ValueError("Scratch reserve already exhausted")
        # Only this child and its ordinary foreground descendants share the group.
        child = subprocess.Popen(a.argv[1:], start_new_session=True, preexec_fn=per_process)
        last = 0
        try:
            while child.poll() is None:
                now = time.monotonic()
                free = shutil.disk_usage(a.scratch).free
                used = rss(child.pid)
                abort = free <= reserve or used >= limit or now >= deadline
                if now - last >= 30 or abort:
                    print(json.dumps(dict(epoch=time.time(), group=child.pid, free=free,
                                          observed_group_rss=used, overrun=abort)), file=log)
                    last = now
                if abort:
                    raise RuntimeError("Owned workload reached its wall, disk or memory limit")
                time.sleep(1)
            return child.returncode
        finally:
            # Kill only members of the newly created group. Also release orphaned
            # descendants if a failing parent exits first.
            try:
                os.killpg(child.pid, signal.SIGTERM)
                time.sleep(1)
                os.killpg(child.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            child.wait()


if __name__ == "__main__":
    raise SystemExit(main())
