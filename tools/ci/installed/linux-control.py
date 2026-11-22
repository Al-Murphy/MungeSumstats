"""Genuine native Linux ownership controls; future hosted execution only."""
import io
import json
import os
from pathlib import Path
import platform
import select
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import linux_job
CHECKS = vars(linux_job)

def proc_entry(proc, pid, parent, start=100, pages=7, comm="worker"):
    path = proc / str(pid)
    path.mkdir(exist_ok=True)
    fields = ["0"] * 22
    fields[0], fields[1], fields[19], fields[21] = "S", str(parent), str(start), str(pages)
    (path / "stat").write_text(f"{pid} ({comm}) " + " ".join(fields) + "\n")
    (path / "comm").write_text(comm + "\n")
    return path

class ProcessBoundary(unittest.TestCase):
    def test_transitive_workers_and_adoptees_only(self):
        with tempfile.TemporaryDirectory() as temp:
            proc = Path(temp)
            proc_entry(proc, 1, 0, pages=100)
            proc_entry(proc, 2, 1, pages=100, comm="sshd")
            proc_entry(proc, 3, 2, pages=100)
            proc_entry(proc, 4, 3, comm="name with ) and spaces")
            proc_entry(proc, 5, 4, pages=11)
            proc_entry(proc, 6, 3, pages=15)
            proc_entry(proc, 7, 6, pages=17)
            owned = CHECKS["owned_processes"](proc, 3)
            self.assertEqual(set(owned), {4, 5, 6, 7})
            self.assertEqual(sum(process.rss_pages for process in owned.values()), 50)

    @unittest.skipUnless(platform.system() == "Linux", "requires Linux subreaping and pidfds")
    def test_new_session_adoptee_gets_killed_and_reaped(self):
        CHECKS["enable_subreaper"]()
        child_code = """\
import os, signal, time
if os.fork():
    os._exit(0)
os.setsid()
signal.signal(signal.SIGTERM, signal.SIG_IGN)
print(os.getpid(), flush=True)
time.sleep(5)
os._exit(0)
"""
        child = subprocess.Popen([sys.executable, "-I", "-B", "-c", child_code],
                                 stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                                 start_new_session=True, text=True)
        self.addCleanup(child.stdout.close)
        log = io.StringIO()
        try:
            self.assertTrue(select.select([child.stdout], [], [], 2)[0], "missing detached-child PID")
            child_pid = int(child.stdout.readline())
            self.assertEqual(child.wait(timeout=2), 0)
            self.assertIn(child_pid, CHECKS["owned_processes"](Path("/proc"), os.getpid()))
            self.assertTrue(CHECKS["quiesce"](child, log))
            self.assertFalse(CHECKS["owned_processes"](Path("/proc"), os.getpid()))
            self.assertTrue(json.loads(log.getvalue())["owned_tree_quiescent"])
        finally:
            CHECKS["quiesce"](child, io.StringIO())

if __name__ == "__main__":
    unittest.main(verbosity=2)
