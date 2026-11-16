import errno
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import guard_common as guard

class GuardTests(unittest.TestCase):
    def measure(self, code, text, error=errno.EPERM):
        def memory(pid):
            if pid == 101: raise OSError(error, 'measurement unavailable')
            return {'resident_bytes':100, 'physical_footprint_bytes':90, 'guard_bytes':100}
        with patch.object(guard.subprocess,'check_output',return_value='100 1\n101 100\n'), patch.object(guard,'process_memory',side_effect=memory), patch.object(guard.subprocess,'run',return_value=SimpleNamespace(returncode=code,stdout=text)):
            return guard.tree_usage(100)
    def test_fresh_absence_can_skip(self): self.assertEqual(self.measure(1,'' )['rss_bytes'],100)
    def test_fresh_zombie_can_skip(self): self.assertEqual(self.measure(0,'Z\n')['rss_bytes'],100)
    def test_native_esrch_can_skip(self): self.assertEqual(self.measure(99,'',errno.ESRCH)['rss_bytes'],100)
    def test_live_unmeasurable_refuses(self):
        with self.assertRaises(OSError): self.measure(0,'S\n')
    def test_failed_liveness_check_refuses(self):
        with self.assertRaises(OSError): self.measure(2,'')
    def test_empty_successful_liveness_check_refuses(self):
        with self.assertRaises(OSError): self.measure(0,'')

if __name__=='__main__':unittest.main(verbosity=2)
