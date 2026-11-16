"""Operator-only, sealed two-build count and logical planning. Never builds an index.

Trust boundary: the local operator supplies an independently reviewed, digest-pinned
profile and (for full scope) a separately digest-pinned admission. HTTP JSON is not
an authority. Profile construction is deliberately absent: expected digests must
come from reviewed producer/plan context, never from convenient replacement bytes.
The worker is a CLI subprocess enclosed by the unchanged process-tree guard.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import stat
import subprocess
import sys
import time

GIB = 1024**3
U64 = 2**64 - 1
ADDRESS_MAX = min(sys.maxsize, 2**63 - 1)
MAX_JSON_BYTES = 16*1024**2
ROLES = {'coordinator', 'counter_source', 'native', 'header', 'guard', 'harness', 'producer'}
OUTPUTS = ('reference.fa', 'variants.tsv', 'contig_aliases.tsv')
FIELDS = '''format_version source_rows source_bytes record_rows site_records allele_records
invalid_positions raw_alt_tokens rejected_alts raw_key_rows raw_sort_bytes final_key_rows
keys_bytes max_source_line_bytes max_key_row_bytes max_record_bytes max_raw_keys_per_record
records_bytes rejected_bytes fai_upper_bytes offsets_bytes'''.split()


class Refusal(RuntimeError):
    pass


def require(condition, message):
    if not condition:
        raise Refusal(message)


def integer(value, maximum=U64):
    require(type(value) is int and 0 <= value <= maximum, 'Nonexact, negative or overflowing integer')
    return value


def add(*values):
    total = sum(integer(value) for value in values)
    return integer(total)


def real_path(value, directory=False):
    """Reject symlinks at every component, including ancestors and dangling links."""
    path = Path(value)
    require(path.is_absolute() and '..' not in path.parts, 'Absolute nontraversing path required')
    for component in [path, *path.parents]:
        require(not component.is_symlink(), 'Symlink path refused')
    try:
        state = path.stat()
    except OSError as exc:
        raise Refusal('Required regular file/publication unavailable') from exc
    require(stat.S_ISDIR(state.st_mode) if directory else stat.S_ISREG(state.st_mode), 'Wrong file type')
    return path


def trusted_publications(expected):
    """Independent operator context, never inferred from unparsed/parsed profile JSON."""
    require(type(expected) is dict and set(expected) == {'GRCh37', 'GRCh38'},
            'Independently trusted expected publication roots required')
    roots = {build: real_path(expected[build], directory=True) for build in ('GRCh37', 'GRCh38')}
    a, b = roots.values()
    require(not a.is_relative_to(b) and not b.is_relative_to(a) and not a.samefile(b),
            'Distinct disjoint expected publication roots required')
    return roots


def outside_publications(value, roots, control=False):
    """Metadata-only pre-open/pre-mkdir check, including aliases and symlink parents."""
    path = Path(value)
    require(path.is_absolute() and '..' not in path.parts, 'Absolute nontraversing placement required')
    for component in [path, *path.parents]:
        require(not component.is_symlink(), 'Symlink placement/parent refused')
        if component.exists():
            require(not any(component.samefile(root) for root in roots.values()),
                    'Publication alias in control/case placement')
    require(not any(path.is_relative_to(root) for root in roots.values()),
            'Control/case must remain outside both expected publications')
    if control:
        real_path(path)
        require(path.stat().st_nlink == 1, 'Hardlinked external control refused')
    return path


def external_placement(expected, profile_path, admission_path=None, case=None):
    roots = trusted_publications(expected)
    outside_publications(profile_path, roots, control=True)
    if admission_path:
        outside_publications(admission_path, roots, control=True)
    if case is not None:
        outside_publications(case, roots)
    return roots


def identity(state):
    return {name: getattr(state, name) for name in
            ('st_dev', 'st_ino', 'st_size', 'st_mtime_ns', 'st_ctime_ns', 'st_mode')}


def pin(path, *, max_bytes=None):
    """Stable digest; a supplied bound is checked before any content open/hash."""
    path = real_path(path)
    before = identity(path.stat())
    if max_bytes is not None:
        require(before['st_size'] <= integer(max_bytes), 'Observed JSON exceeds finite bound')
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        require(identity(os.fstat(handle.fileno())) == before, 'File replaced before hash')
        remaining = before['st_size']
        while remaining:
            chunk = handle.read(min(1024*1024, remaining))
            require(bool(chunk), 'File truncated while hashing')
            digest.update(chunk)
            remaining -= len(chunk)
        require(identity(os.fstat(handle.fileno())) == before, 'File changed while hashing')
    require(identity(path.stat()) == before, 'File changed/replaced while hashing')
    return {'path': str(path), 'size_bytes': before['st_size'], 'sha256': digest.hexdigest(),
            'identity': before}


def expected_pin(specification, pinner=pin, *, max_bytes=None):
    integer(specification['size_bytes'], min(ADDRESS_MAX, max_bytes) if max_bytes is not None else ADDRESS_MAX)
    digest = specification['sha256']
    require(type(digest) is str and len(digest) == 64 and all(c in '0123456789abcdef' for c in digest),
            'Expected SHA256 required from reviewed context')
    actual = pinner(specification['path'], max_bytes=max_bytes) if max_bytes is not None else pinner(specification['path'])
    require(actual['size_bytes'] == specification['size_bytes'] and actual['sha256'] == digest,
            'Reviewed file identity mismatch: '+specification['path'])
    return actual


def bounded_json_bytes(path):
    """Check observed cap before parse read as well as before its enclosing hashes."""
    path = real_path(path)
    before = identity(path.stat())
    require(before['st_size'] <= MAX_JSON_BYTES, 'Observed JSON exceeds finite bound')
    with path.open('rb') as handle:
        require(identity(os.fstat(handle.fileno())) == before, 'JSON replaced before parse read')
        raw = handle.read(before['st_size'])
        require(len(raw) == before['st_size'], 'JSON truncated during parse read')
        require(identity(os.fstat(handle.fileno())) == before, 'JSON changed during parse read')
    require(identity(path.stat()) == before, 'JSON replaced during parse read')
    return raw


def json_load(path):
    def unique(pairs):
        result = {}
        for key, value in pairs:
            require(key not in result, 'Duplicate JSON key')
            result[key] = value
        return result
    try:
        return json.loads(bounded_json_bytes(path).decode('utf-8'), object_pairs_hook=unique,
                          parse_constant=lambda x: (_ for _ in ()).throw(Refusal('Nonfinite JSON')))
    except (ValueError, UnicodeError) as exc:
        raise Refusal('Invalid/truncated JSON') from exc


def pinned_json(specification):
    before = expected_pin(specification, max_bytes=MAX_JSON_BYTES)
    value = json_load(specification['path'])
    require(expected_pin(specification, max_bytes=MAX_JSON_BYTES) == before, 'Control changed during parse')
    return value


def check_guard(guard):
    require(guard.get('status') == 'success' and type(guard.get('exit_code')) is int
            and guard['exit_code'] == 0 and guard.get('reason') is None
            and guard.get('natural_oom') is False, 'Producer/count guard not successful')


def controls(profile):
    return [*profile['artifacts'].values(),
            *(b[k] for b in profile['builds'] for k in ('source_plan', 'receipt', 'producer_guard'))]


def control_snapshots(profile):
    return ([expected_pin(spec) for spec in profile['artifacts'].values()]
            + [expected_pin(build[role], max_bytes=MAX_JSON_BYTES)
               for build in profile['builds'] for role in ('source_plan', 'receipt', 'producer_guard')])


def validate_seals(profile, admission=None, *, expected_publications=None):
    """Validate BOTH producers and small trusted controls before any output-file open."""
    roots = trusted_publications(expected_publications)
    require(type(profile.get('builds')) is list and len(profile['builds']) == 2
            and [b.get('build') for b in profile['builds']] == ['GRCh37', 'GRCh38']
            and all(b.get('final_dir') == str(roots[b['build']]) for b in profile['builds']),
            'Parsed publication roots differ from independent expected context')
    require(type(profile.get('schema_version')) is int and profile['schema_version'] == 1, 'Wrong profile schema')
    context = profile.get('context')
    require(context in ('tiny_fixture', 'full_reference'), 'Explicit context required')
    if context == 'full_reference':
        require(admission is not None and admission.get('action') == 'sealed-full-count-only'
                and admission.get('independent_review') == 'GREEN'
                and admission.get('source_commit') and admission.get('profile_sha256'),
                'Separate explicit pinned root admission required')
        require(type(profile.get('coordinator_source_commit')) is str
                and len(profile['coordinator_source_commit']) == 40
                and all(c in '0123456789abcdef' for c in profile['coordinator_source_commit'])
                and admission['source_commit'] == profile['coordinator_source_commit'],
                'Admission/frozen coordinator source mismatch')
    require(type(profile.get('format_version')) is int and profile['format_version'] == 3, 'Only reviewed format 3 planning')
    require(set(profile['artifacts']) == ROLES, 'Missing/extra reviewed artifact roles')
    limits = profile['limits']
    require(limits['tree_memory_bytes'] == 8*GIB and type(limits['tree_memory_bytes']) is int
            and limits['minimum_free_bytes'] == 160*GIB and type(limits['minimum_free_bytes']) is int,
            '8GiB/160GiB limits must remain unchanged')
    require(0 < integer(limits['owned_case_bytes']) <= GIB, 'Finite owned output cap required')
    require(type(limits['deadline_seconds']) in (int, float)
            and math.isfinite(limits['deadline_seconds']) and 0 < limits['deadline_seconds'] <= 7*86400,
            'Finite positive deadline required')
    integer(profile['manifest_upper_bytes'], ADDRESS_MAX)
    require(profile['manifest_upper_bytes'] > 0 and bool(profile['manifest_bound_basis']), 'Explicit bounded RDS allowance required')
    integer(profile['allocation_allowance_bytes']); integer(profile['concurrent_reservation_bytes'])
    builds = profile['builds']
    require(len(builds) == 2 and [b['build'] for b in builds] == ['GRCh37', 'GRCh38'], 'Both ordered build seals required')
    require(len({b['final_dir'] for b in builds}) == 2, 'Distinct final directories required')
    # A malformed local profile must not disguise a prepared input as an artifact
    # or plan/guard control and open it before both producers are sealed.
    finals = [real_path(b['final_dir'], directory=True) for b in builds]
    # Small pre-seal controls must have independent identities. Counted data is
    # deliberately exempt: legitimate scientific hardlinks require no copies.
    controls_with_identity_rule = [(spec, False) for spec in profile['artifacts'].values()]
    controls_with_identity_rule += [(b[k], True) for b in builds for k in ('source_plan', 'producer_guard')]
    for control, single_link in controls_with_identity_rule:
        path = real_path(control['path'])
        require(not any(path.is_relative_to(final) for final in finals),
                'Source/artifact control inside a prepared publication')
        if single_link:
            require(path.stat().st_nlink == 1, 'Hardlinked small source-plan/producer-guard control refused')
    require(Path(profile['artifacts']['coordinator']['path']) == Path(__file__).resolve(), 'Coordinator self pin mismatch')
    # Authenticate BOTH guards before opening any plan, producer receipt or
    # artifact content. Metadata-only publication/receipt placement checks remain
    # permitted here; scientific input hashing stays after the complete seal gate.
    for build in builds:
        final = real_path(build['final_dir'], directory=True)
        require(all(not p.startswith('.') and '-pending-' not in p for p in final.parts[1:]), 'Pending final directory refused')
        require(Path(build['receipt']['path']) == final/'receipts.json', 'Receipt outside published final directory')
        guard = pinned_json(build['producer_guard'])
        check_guard(guard)
        require(guard.get('command') == build['producer_command'], 'Producer command mismatch')
        require(build['source_plan']['path'] in guard['command'] and str(final) in guard['command'],
                'Guard does not bind reviewed plan and publication')
    result = []
    for build in builds:
        final = real_path(build['final_dir'], directory=True)
        plan = pinned_json(build['source_plan'])
        receipt = pinned_json(build['receipt'])
        for name in ('build', 'release'):
            require(receipt.get(name) == plan.get(name) == build[name], 'Build/release mismatch')
        for name in ('scope', 'assembly', 'vcf_headers', 'sources'):
            require(name in plan, 'Incomplete reviewed source plan: '+name)
        for name, expected in plan.items():
            require(receipt.get(name) == expected, 'Reviewed source plan mismatch: '+name)
        require((plan['scope'] == 'tiny_fixture') == (context == 'tiny_fixture'), 'Fixture/full context mismatch')
        require(receipt.get('preparation_version') == build['preparation_version'] == 'biodb-reference-stream-4',
                'Producer version mismatch')
        require(receipt.get('reference_format') == 'mss-preservation-2', 'Unexpected preparer TSV contract')
        require(receipt.get('alias_history_status') == 'complete_supplied_history', 'Incomplete history')
        sources = receipt['sources']
        require(set(sources) == {'fasta', 'vcf', 'assembly_report', 'merged', 'withdrawn'}, 'Complete reviewed source roles required')
        for role, source in sources.items():
            require(source.get('build') == build['build'], 'Source build mismatch')
            integer(source['size_bytes'], ADDRESS_MAX)
            require(type(source.get('sha256')) is str and len(source['sha256']) == 64, 'Source digest missing')
        for role in ('merged', 'withdrawn'):
            require(sources[role].get('completeness') == 'complete_source'
                    and sources[role]['release'] == sources['vcf']['release'], 'Incomplete/current history mismatch')
        rows = integer(receipt['counts']['emitted_rows'], ADDRESS_MAX)
        require(rows == integer(build['expected_emitted_rows'], ADDRESS_MAX), 'Expected emitted row mismatch')
        require(receipt.get('reference_fasta') == str(final/'reference.fa')
                and receipt.get('contig_aliases') == str(final/'contig_aliases.tsv'), 'Derived output paths mismatch')
        inputs = []
        for name in OUTPUTS:
            require(name in receipt['outputs'], 'Missing prepared output: '+name)
            output = receipt['outputs'][name]
            integer(output['size_bytes'], ADDRESS_MAX)
            require(type(output.get('sha256')) is str and len(output['sha256']) == 64, 'Output digest missing')
            path = real_path(final/name)
            require(path.stat().st_size == output['size_bytes'], 'Output size mismatch')
            inputs.append(dict(path=str(path), size_bytes=output['size_bytes'], sha256=output['sha256']))
        result.append({'build': build['build'], 'inputs': inputs, 'expected_emitted_rows': rows})
    for artifact in profile['artifacts'].values():
        expected_pin(artifact)
    return result


def alias_upper(path):
    """R quote=TRUE/qmethod=double UTF8 bytes, using conservative CRLF framing."""
    csv.field_size_limit(ADDRESS_MAX)
    total = 0
    with Path(path).open(encoding='utf-8', newline='') as handle:
        reader = csv.reader(handle, delimiter='\t', strict=True)
        require(next(reader, None) == ['alias', 'canonical_contig'], 'Invalid producer alias header')
        for fields in [['alias', 'canonical_contig']]:
            total = add(total, sum(len(s.replace('"', '""').encode('utf-8'))+2 for s in fields), 3)
        for fields in reader:
            require(len(fields) == 2 and all(fields), 'Invalid producer alias row')
            total = add(total, sum(len(s.replace('"', '""').encode('utf-8'))+2 for s in fields), 3)
    return total


def logical_size(counts, aliases, manifest):
    t,r,u,k,n,j,f = [integer(counts[name], ADDRESS_MAX) for name in
                   ('source_bytes','records_bytes','raw_sort_bytes','keys_bytes','final_key_rows','rejected_bytes','fai_upper_bytes')]
    offsets = integer(8*n, ADDRESS_MAX)
    c = add(f, aliases, manifest, 19)
    retained = add(t,r,k,offsets,j,c)
    phase = max(retained, add(r,u,u,j,c), add(r,u,k,offsets,j,c))
    return {'small_files_upper_bytes': c, 'retained_bytes': retained, 'phase_bytes': phase}


def capacity(sizes, free, allocation, concurrent):
    a,b = sizes
    free,allocation,concurrent = [integer(x) for x in (free,allocation,concurrent)]
    peaks = {'GRCh37_then_GRCh38': max(a['phase_bytes'],add(a['retained_bytes'],b['phase_bytes']),add(a['retained_bytes'],b['retained_bytes'])),
             'GRCh38_then_GRCh37': max(b['phase_bytes'],add(b['retained_bytes'],a['phase_bytes']),add(a['retained_bytes'],b['retained_bytes']))}
    reserved = add(170*GIB, allocation, concurrent)
    return {'order_peak_bytes': peaks, 'filesystem_free_after_both_inputs_bytes': free,
            'planning_reserve_bytes': 170*GIB, 'allocation_allowance_bytes': allocation,
            'concurrent_reservation_bytes': concurrent,
            'orders_fit_with_reserves': {order: add(peak,reserved) <= free for order,peak in peaks.items()},
            'full_reference_admitted': False, 'index_builds_executed': 0,
            'decision': 'no_admission', 'future_build_gate': 'Separate reviewed native build/memory profile, maximum-line/sort-buffer assessment, fresh free measurement and post-first-build remeasurement required'}


def validate_counts(value, seal):
    require(set(value) == set(FIELDS), 'Missing/extra native count fields')
    for name in FIELDS:
        integer(value[name], ADDRESS_MAX)
    require(value['format_version'] == 3 and value['source_rows'] == seal['expected_emitted_rows'], 'Native format/emitted row mismatch')
    require(value['source_bytes'] == seal['inputs'][1]['size_bytes'], 'Native source byte mismatch')
    require(value['offsets_bytes'] == integer(8*value['final_key_rows'], ADDRESS_MAX), 'Offset count mismatch')
    require(value['site_records'] + value['invalid_positions'] == value['source_rows']
            and value['site_records'] + value['allele_records'] == value['record_rows']
            and value['raw_key_rows'] >= value['final_key_rows']
            and value['raw_sort_bytes'] >= value['keys_bytes'], 'Inconsistent native counts')


    # Necessary structural invariants of the frozen typed writer/counter. These
    # validate malformed derived responses; they never filter scientific records.
    rows, sites, alleles, invalid = [value[k] for k in
                                    ('source_rows', 'site_records', 'allele_records', 'invalid_positions')]
    raw_alts, rejected = value['raw_alt_tokens'], value['rejected_alts']
    require(alleles + rejected <= raw_alts <= value['source_bytes'],
            'ALT disposition/raw token count inconsistent with source bytes')
    require(invalid != 0 or raw_alts == alleles + rejected,
            'ALT disposition must be exact when no positions are invalid')
    require(sites != 0 or (alleles == 0 and rejected == 0), 'Allele dispositions require site records')
    require(rows != 0 or raw_alts == 0, 'Empty source has ALT tokens')
    # split("") has zero ALT tokens: empty ALT rows and invalid empty-ALT rows
    # legitimately violate a tempting but incorrect raw_alts >= source_rows rule.
    header_bytes = len(b'chr\tpos\tref\talt\tid\taliases')
    line_max = value['max_source_line_bytes']
    require(value['source_bytes'] >= header_bytes + 6*rows, 'Source bytes below six-field framing minimum')
    if rows:
        require(5 <= line_max <= value['source_bytes']-header_bytes-1-6*(rows-1),
                'Source line maximum inconsistent with source bytes')
    else:
        require(line_max == 0 and header_bytes <= value['source_bytes'] <= header_bytes+2,
                'Empty source framing/maxima mismatch')

    def population(label, count, total, maximum, minimum):
        if count == 0:
            require(total == 0 and maximum == 0, label+' zero population has bytes/maxima')
        else:
            require(minimum <= maximum and maximum+minimum*(count-1) <= total <= maximum*count,
                    label+' population bytes/maxima inconsistent')

    population('Record', value['record_rows'], value['records_bytes'], value['max_record_bytes'], 7)
    require(value['records_bytes'] >= 7*sites + 13*alleles, 'Typed record framing below minimum')
    require(alleles == 0 or value['max_record_bytes'] >= 13, 'Allele record maximum below framing minimum')
    raw, final = value['raw_key_rows'], value['final_key_rows']
    key_max = value['max_key_row_bytes']
    population('Raw key', raw, value['raw_sort_bytes'], key_max, 6)
    population('Final key', final, value['keys_bytes'], key_max, 6)
    require(final >= value['record_rows'] and raw >= sites+2*alleles,
            'Each record requires keys; each allele emits two raw coordinate keys')
    raw_max = value['max_raw_keys_per_record']
    if alleles:
        require(2 <= raw_max <= raw-sites-2*(alleles-1)
                and raw <= sites+alleles*raw_max, 'Raw per-record key maximum inconsistent')
    else:
        require(raw == final == sites and raw_max == (1 if sites else 0),
                'Site-only/empty key population mismatch')
    duplicates = raw-final
    require(6*duplicates <= value['raw_sort_bytes']-value['keys_bytes'] <= key_max*duplicates,
            'Raw/final byte difference inconsistent with duplicate rows')

    newline = len(os.linesep.encode())
    rejected_header = len(b'row\tstatus')+newline
    digits = len(str(rows+1))
    invalid_status = len('invalid_position')
    alt_statuses = [len('unsupported_allele'), len('reference_mismatch')]
    rejection_min = rejected_header + invalid*(2+invalid_status+newline) + rejected*(2+min(alt_statuses)+newline)
    rejection_max = rejected_header + invalid*(digits+1+invalid_status+newline) + rejected*(digits+1+max(alt_statuses)+newline)
    require(rejection_min <= value['rejected_bytes'] <= rejection_max,
            'Rejected header/disposition byte framing mismatch')
    require(value['fai_upper_bytes'] >= 85+newline, 'Nonempty FASTA requires positive FAI framing bound')


def write_once(path, value, limit_bytes=None):
    encoded = (json.dumps(value, indent=2, allow_nan=False)+'\n').encode('utf-8')
    if limit_bytes is not None:
        require(len(encoded) <= integer(limit_bytes), 'Final receipt exceeds reserved output allowance')
    with Path(path).open('xb') as handle:
        handle.write(encoded)


def count_worker(profile, case, admission=None, input_pinner=pin, *, expected_publications=None):
    roots = trusted_publications(expected_publications)
    case = outside_publications(case, roots)
    seals = validate_seals(profile, admission, expected_publications=expected_publications)
    case.mkdir()  # Exclusive, never overwrite/retry a case.
    control_before = control_snapshots(profile)
    before = [[expected_pin(s, input_pinner) for s in seal['inputs']] for seal in seals]
    # Both derived FASTA closures and prepared TSVs now exist and are authenticated.
    free = shutil.disk_usage(case).free
    require(free >= 160*GIB, '160GiB filesystem hard floor unavailable')
    results = []
    for seal in seals:
        stem = case/seal['build']
        stdout, stderr = Path(str(stem)+'.stdout'), Path(str(stem)+'.stderr')
        with stdout.open('xb') as out, stderr.open('xb') as err:
            child = subprocess.run([profile['artifacts']['native']['path'], seal['inputs'][0]['path'],
                                    seal['inputs'][1]['path'], '3'], stdout=out, stderr=err)
        require(child.returncode == 0 and stderr.stat().st_size == 0, 'Native failure or unexpected successful stderr')
        require(0 < stdout.stat().st_size <= 65536, 'Native output missing/over-limit')
        counts = json_load(stdout)
        validate_counts(counts, seal)
        aliases = alias_upper(seal['inputs'][2]['path'])
        sizes = logical_size(counts, aliases, profile['manifest_upper_bytes'])
        results.append({'build': seal['build'], 'counts': counts, 'sizes': sizes,
                        'quoted_alias_upper_bytes': aliases})
    after = [[expected_pin(s, input_pinner) for s in seal['inputs']] for seal in seals]
    require(before == after, 'Counted input changed identity/digest')
    require(control_before == control_snapshots(profile), 'Pinned control/source/binary changed')
    planning = capacity([r['sizes'] for r in results], free, profile['allocation_allowance_bytes'], profile['concurrent_reservation_bytes'])
    result = dict(planning, status='counted_planning_only', context=profile['context'], builds=results,
                  independently_expected_publications={b: str(root) for b, root in roots.items()},
                  inputs_before=before, inputs_after=after, controls_before=control_before,
                  controls_after=control_before, manifest_bound_basis=profile['manifest_bound_basis'],
                  small_file_framing='FAI conservative native bound; aliases quote-all/doubled-quotes UTF8 with CRLF upper bound; marker 19 bytes; explicit operator RDS allowance')
    write_once(case/'worker.json', result)
    return result


def load_profile(path, expected_sha, admission_path=None, admission_sha=None, *, expected_publications=None):
    external_placement(expected_publications, path, admission_path)
    require(bool(expected_sha), 'Reviewed profile digest required')
    profile_spec = {'path': str(Path(path).absolute()), 'size_bytes': Path(path).stat().st_size, 'sha256': expected_sha}
    profile = pinned_json(profile_spec)
    admission = None
    admission_spec = None
    if admission_path:
        require(bool(admission_sha), 'Independent admission digest required')
        admission_spec = {'path': str(Path(admission_path).absolute()), 'size_bytes': Path(admission_path).stat().st_size, 'sha256': admission_sha}
        admission = pinned_json(admission_spec)
        require(admission.get('profile_sha256') == expected_sha, 'Admission/profile digest mismatch')
        require(profile['context'] == 'full_reference', 'Full admission cannot promote fixture context')
    validate_seals(profile, admission, expected_publications=expected_publications)
    return profile, admission, [profile_spec] + ([admission_spec] if admission_spec else [])


def launch(profile_path, expected_sha, case, admission_path=None, admission_sha=None, *, expected_publications=None):
    """Monitor complete worker: authentication recheck, hashes, native, post-hashes."""
    roots = external_placement(expected_publications, profile_path, admission_path, case)
    profile, admission, external = load_profile(profile_path, expected_sha, admission_path, admission_sha,
                                                expected_publications=expected_publications)
    case = outside_publications(case, roots)
    real_path(case.parent, directory=True)
    case.mkdir()
    external_before = [expected_pin(s, max_bytes=MAX_JSON_BYTES) for s in external]
    artifacts_before = control_snapshots(profile)
    # Guard and memory harness are pinned before importing the reviewed implementation.
    guard_path = Path(profile['artifacts']['guard']['path'])
    require(Path(profile['artifacts']['harness']['path']).parent == guard_path.parent, 'Guard/harness directory mismatch')
    sys.path.insert(0, str(guard_path.parent))
    module_spec = importlib.util.spec_from_file_location('sealed_reviewed_guard', guard_path)
    guard = importlib.util.module_from_spec(module_spec)
    module_spec.loader.exec_module(guard)
    # The public CLI has no unmonitored worker switch. This internal entry is only
    # composed here as the owned child of the reviewed guard.
    command = [sys.executable, '-c',
               "import runpy,sys; n=runpy.run_path(sys.argv[1]); n['_worker_entry'](*sys.argv[2:])",
               str(Path(__file__).resolve()), str(Path(profile_path).absolute()), expected_sha,
               str(case/'worker'), str(Path(admission_path).absolute()) if admission_path else '', admission_sha or '',
               str(roots['GRCh37']), str(roots['GRCh38'])]
    limits = profile['limits']; started = time.monotonic()
    receipt = guard.guarded(command, case, case/'guard.json', case/'worker.log',
                            {'max_tree_memory_bytes': 8*GIB, 'min_filesystem_free_bytes': 160*GIB,
                             'poll_seconds': .05, 'max_case_disk_bytes': limits['owned_case_bytes']},
                            started+limits['deadline_seconds'], os.environ.copy())
    check_guard(receipt)
    require(time.monotonic()-started <= limits['deadline_seconds'], 'Post-exit deadline censor')
    require(guard.disk_bytes(case)+65536 <= limits['owned_case_bytes'], 'Post-exit output cap (including receipt allowance)')
    require(shutil.disk_usage(case).free >= 160*GIB, 'Post-exit hard free floor')
    require(external_before == [expected_pin(s, max_bytes=MAX_JSON_BYTES) for s in external], 'Profile/admission changed')
    require(artifacts_before == control_snapshots(profile), 'Reviewed controls changed')
    summary = json_load(case/'worker'/'worker.json')
    summary.update(guard=receipt, external_control_pins=external_before,
                   qualification='Tiny local count only' if profile['context']=='tiny_fixture' else 'Explicitly admitted full count only; no index admission')
    write_once(case/'summary.json', summary, limit_bytes=65536)
    return summary


def _worker_entry(profile_path, expected_sha, case, admission_path, admission_sha, grch37=None, grch38=None):
    expected = {'GRCh37': grch37, 'GRCh38': grch38} if grch37 is not None and grch38 is not None else None
    external_placement(expected, profile_path, admission_path or None, case)
    profile, admission, external = load_profile(profile_path, expected_sha, admission_path or None, admission_sha or None,
                                                expected_publications=expected)
    before = [expected_pin(s, max_bytes=MAX_JSON_BYTES) for s in external]
    result = count_worker(profile, case, admission, expected_publications=expected)
    require(before == [expected_pin(s, max_bytes=MAX_JSON_BYTES) for s in external], 'Profile/admission changed during worker')
    print(json.dumps({'status': result['status'], 'decision': 'no_admission',
                      'full_reference_admitted': False, 'index_builds_executed': 0}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profile', required=True)
    parser.add_argument('--expected-profile-sha256', required=True)
    parser.add_argument('--case', required=True)
    parser.add_argument('--expected-grch37-publication', required=True, help='Independent trusted root, outside profile JSON')
    parser.add_argument('--expected-grch38-publication', required=True, help='Independent trusted root, outside profile JSON')
    parser.add_argument('--full-admission')
    parser.add_argument('--expected-admission-sha256')
    args = parser.parse_args()
    os.umask(0o077)
    try:
        result = launch(args.profile, args.expected_profile_sha256, args.case, args.full_admission, args.expected_admission_sha256,
                        expected_publications={'GRCh37': args.expected_grch37_publication,
                                               'GRCh38': args.expected_grch38_publication})
        print(json.dumps({'status': result['status'], 'decision': 'no_admission',
                          'full_reference_admitted': False, 'index_builds_executed': 0}))
    except Exception as exc:
        print(json.dumps({'status': 'refused', 'decision': 'no_admission', 'full_reference_admitted': False,
                          'index_builds_executed': 0, 'reason': type(exc).__name__+': '+str(exc)}), file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
