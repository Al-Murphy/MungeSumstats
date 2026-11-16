"""Adversarial selfcheck of the benchmark oracle; no large benchmark runs."""
import copy
import csv
import hashlib
import io
import json
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from verify import FIELDS, Gold, verify
from bench import approved


class OracleTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.fasta = self.root / 'ref.fa'
        self.fasta.write_text('>22\nACGT\n')
        self.fai = self.root / 'ref.fai'
        self.fai.write_text('22\t4\t4\t4\t5\n')
        self.reference = self.root / 'reference.tsv'
        self.reference.write_text('chr\tpos\tref\talt\tid\taliases\n22\t1\tA\tC\trs1\told1\n')
        self.source = self.root / 'input.tsv'
        self.source.write_text('SNP\tCHR\tBP\tA1\tA2\tBETA\tFRQ\tSE\tP\nrs1\t22\t1\tA\tC\t.125\t.25\t.05\t.01\nrs1\t22\t1\tC\tA\t-.125\t.75\t.05\t.01\n')
        self.output = self.root / 'output'
        self.output.mkdir()
        for name in ('source.tsv', 'source-artifact'):
            (self.output / name).write_bytes(self.source.read_bytes())
        self.rows = []
        for i, status in enumerate(('match', 'swap'), 1):
            row = {key: 'NA' for key in FIELDS}
            row.update(status=status, chr='22', pos='1', ref='A', alt='C', id='rs1', effect='0.125', eaf='0.25', swapped='FALSE' if i == 1 else 'TRUE', complemented='FALSE', candidate_count='1', entity_type='small_variant', entity_id='rs1', normalized_build='GRCh37', source_row=str(i))
            self.rows.append(row)
        self.profile = self.root / 'profile.json'
        self.profile.write_text(json.dumps({'stages': {'tiny': {'input': {'path': str(self.source), 'sha256': hashlib.sha256(self.source.read_bytes()).hexdigest()}, 'rows': 2}}, 'reference': {'source': str(self.reference), 'fasta': str(self.fasta), 'fai': str(self.fai), 'build': 'GRCh37'}}))

    def tearDown(self):
        self.temporary.cleanup()

    def write_rows(self, rows=None, fields=None):
        with (self.output / 'normalized.tsv').open('w') as handle:
            writer = csv.DictWriter(handle, fieldnames=fields or FIELDS, delimiter='\t')
            writer.writeheader()
            writer.writerows(self.rows if rows is None else rows)

    def check(self, label='result'):
        with redirect_stdout(io.StringIO()):
            return verify(self.profile, 'tiny', {'mss7': str(self.output)}, self.root / (label + '.json'))

    def test_exact_match_and_swap_retention(self):
        self.write_rows()
        result = self.check()
        self.assertTrue(result['passed'])
        self.assertEqual(result['versions']['mss7']['counts']['all_19_fields_correct_rows'], 2)
        self.assertEqual(result['versions']['mss7']['statuses'], {'match': 1, 'swap': 1})

    def test_every_normalized_field_corruption_is_detected(self):
        for field in FIELDS:
            with self.subTest(field=field):
                rows = copy.deepcopy(self.rows)
                rows[0][field] = '0.5' if field in ('effect', 'eaf') else '2' if field in ('pos', 'source_row', 'candidate_count') else 'CORRUPT'
                self.write_rows(rows)
                with self.assertRaises(AssertionError):
                    self.check(field)

    def test_short_and_extra_rows_are_detected(self):
        for rows in (self.rows[:1], self.rows + self.rows[:1]):
            self.write_rows(rows)
            with self.assertRaises(AssertionError):
                self.check()

    def test_reordered_schema_is_detected(self):
        self.write_rows(fields=list(reversed(FIELDS)))
        with self.assertRaises(AssertionError):
            self.check()

    def test_each_original_copy_is_verified(self):
        self.write_rows()
        for name in ('source.tsv', 'source-artifact'):
            original = (self.output / name).read_bytes()
            (self.output / name).write_bytes(original + b'corrupt')
            with self.assertRaises(AssertionError):
                self.check()
            (self.output / name).write_bytes(original)

    def test_missing_frequency_and_direction_corruption_is_detected(self):
        for field, value in (('eaf', 'NA'), ('effect', '-0.125'), ('swapped', 'FALSE')):
            rows = copy.deepcopy(self.rows)
            rows[1][field] = value
            self.write_rows(rows)
            with self.assertRaises(AssertionError):
                self.check()

    def test_source_generation_contract_is_checked(self):
        self.write_rows()
        content = self.source.read_text().replace('\t.05\t', '\t.06\t')
        self.source.write_text(content)
        for name in ('source.tsv', 'source-artifact'):
            (self.output / name).write_bytes(self.source.read_bytes())
        profile = json.loads(self.profile.read_text())
        profile['stages']['tiny']['input']['sha256'] = hashlib.sha256(self.source.read_bytes()).hexdigest()
        self.profile.write_text(json.dumps(profile))
        with self.assertRaises(AssertionError):
            self.check()

    def test_complex_multiallelic_duplicate_id_or_palindromic_reference_is_refused(self):
        rows = ['22\t1\tA\tC,G\trs1\t', '22\t1\tA\tAC\trs1\t', '22\t1\tA\tT\trs1\t', '22\t1\tA\tC\trs1\t\n22\t2\tC\tG\trs1\t']
        for row in rows:
            self.reference.write_text('chr\tpos\tref\talt\tid\taliases\n' + row + '\n')
            with self.assertRaises(AssertionError):
                Gold(self.reference, self.fasta, self.fai, 'GRCh37')

    def test_extra_undeclared_output_column_is_refused(self):
        self.write_rows()
        p = self.output / 'normalized.tsv'
        lines = p.read_text().splitlines()
        lines[1] += '\textra'
        p.write_text('\n'.join(lines) + '\n')
        with self.assertRaises(AssertionError):
            self.check()

    def test_duplicate_identical_reference_rows_are_refused(self):
        self.reference.write_text(self.reference.read_text() + '22\t1\tA\tC\trs1\told1\n')
        with self.assertRaises(AssertionError):
            Gold(self.reference, self.fasta, self.fai, 'GRCh37')

    def approval_case(self, *, red=False, stale=False):
        revision = '5e5e5bfd1a8ccb4bddfb9520e0e0997d4691708e'
        core = 'a' * 64
        dll = 'b' * 64
        profile = {'sources': {'mss7': {'production_revision': revision, 'files': {'src/preservation_core.h': {'sha256': core}}}}, 'libraries': {'mss7': {'dll': {'sha256': dll}}}}
        profile_path = self.root / 'approval-profile.json'
        profile_path.write_text(json.dumps(profile))
        review = self.root / 'review.md'
        review.write_text(('**Verdict: RED. NOT GREEN.**' if red else '**Verdict: GREEN for the documented scientific contract of frozen `' + revision + '`.**') + '\nCore ' + core + '\n')
        suite = self.root / 'suite.json'
        suite.write_text(json.dumps({'failed_assertions': 0, 'errors': 0, 'passed_assertions': 1, 'installed_dll_sha256': dll, 'source_revision': 'old-r-only-revision' if stale else revision}))
        pin = lambda p: {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}
        approval = self.root / 'approval.json'
        approval.write_text(json.dumps({'approved': True, 'status': 'SCIENCE_GREEN', 'profile_sha256': hashlib.sha256(profile_path.read_bytes()).hexdigest(), 'authorized_stages': ['1M'], 'candidate_core_sha256': core, 'scientific_review': pin(review), 'installed_suite': pin(suite)}))
        return profile_path, profile, approval

    def test_exact_scientific_green_approval_is_accepted(self):
        path, profile, approval = self.approval_case()
        self.assertEqual(approved(path, profile, approval, '1M')['path'], str(approval.resolve()))

    def test_red_not_green_review_is_refused(self):
        path, profile, approval = self.approval_case(red=True)
        with self.assertRaises(RuntimeError):
            approved(path, profile, approval, '1M')

    def test_stale_installed_suite_source_revision_is_refused(self):
        path, profile, approval = self.approval_case(stale=True)
        with self.assertRaises(RuntimeError):
            approved(path, profile, approval, '1M')

    def test_reference_fasta_disagreement_is_refused(self):
        self.reference.write_text('chr\tpos\tref\talt\tid\taliases\n22\t1\tG\tT\trs1\t\n')
        with self.assertRaises(AssertionError):
            Gold(self.reference, self.fasta, self.fai, 'GRCh37')


if __name__ == '__main__':
    unittest.main(verbosity=2)
