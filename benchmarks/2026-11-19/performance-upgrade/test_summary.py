"""Portable benchmark-summary regressions; no R, compiler or native workload."""
import copy
import importlib.util
import json
from pathlib import Path
import sys
import unittest

SOURCE = Path(__file__).resolve().parent
sys.path.insert(0, str(SOURCE))
spec = importlib.util.spec_from_file_location('benchmark_run', SOURCE/'run.py')
run = importlib.util.module_from_spec(spec)
spec.loader.exec_module(run)

def observations(child='NA'):
    rows = []
    for pair in (1, 2):
        for role in run.ROLES:
            metrics = {f'{name}.{field}': value
                       for name in ('package_setup', 'reference_open', 'processing')
                       for field, value in [('elapsed','1.25'), ('user.self','.25'),
                                            ('sys.self','.125'), ('user.child',child), ('sys.child',child)]}
            process = dict(role=role, peak_rss_bytes=1024, rss_complete=True,
                           wall_seconds=1.5, user_seconds=.3, system_seconds=.2)
            rows.append(dict(pair=pair, metrics=metrics, process=process))
    return rows

def summarize(rows):
    return run.summarize_case({'pairs':2, 'rows':10, 'target':False}, rows, suite='diagnostic')

class Summary(unittest.TestCase):
    def test_windows_unavailable_child_cpu_is_null_without_losing_measured_times(self):
        rows = observations()
        original = copy.deepcopy(rows)
        result = summarize(rows)
        for role in run.ROLES:
            for name in ('package_setup','reference_open','processing'):
                self.assertEqual(result['medians'][role][name],
                                 {'elapsed':1.25,'user.self':.25,'sys.self':.125,
                                  'user.child':None,'sys.child':None})
        self.assertEqual(result['ratio_of_median_processing_seconds'], 1.0)
        self.assertEqual(result['baseline_over_candidate_by_pair'], [1.0,1.0])
        self.assertEqual(rows, original)
        json.dumps(result, allow_nan=False)

    def test_available_zero_and_nonzero_child_cpu_remain_measured_values(self):
        for token, expected in [('0',0.0), ('.125',.125)]:
            with self.subTest(token=token):
                result = summarize(observations(token))
                for role in run.ROLES:
                    self.assertEqual(result['medians'][role]['processing']['user.child'], expected)
                    self.assertEqual(result['medians'][role]['processing']['sys.child'], expected)

    def test_nonfinite_invalid_or_incomplete_metrics_still_fail(self):
        for field in ('elapsed','user.self','sys.self','user.child','sys.child'):
            for value in ('NaN','Inf','-Inf','','bad'):
                rows=observations('0');rows[0]['metrics']['processing.'+field]=value
                with self.subTest(field=field,value=value), self.assertRaises(ValueError):
                    summarize(rows)
        for field in ('elapsed','user.self','sys.self'):
            rows=observations('0');rows[0]['metrics']['processing.'+field]='NA'
            with self.subTest(required=field), self.assertRaises(ValueError):summarize(rows)
        rows=observations('0');rows[0]['metrics']['processing.user.child']='NA'
        with self.assertRaises(ValueError):summarize(rows)
        rows=observations('0');del rows[0]['metrics']['processing.elapsed']
        with self.assertRaises(KeyError):summarize(rows)

if __name__ == '__main__':
    unittest.main(verbosity=2)
