"""Same-handle companion using frozen worker settings and verification."""
import csv
import json
from pathlib import Path

import driver


ANCHOR = '# Match the prior profile\'s fresh-process preparation for both revisions.\nheap_before <- gc(reset = TRUE)\n'
ROLES = driver.ROLES


def render_worker(source, prime):
    driver.require(source.count(ANCHOR) == 1, 'Frozen worker heap-reset anchor differs')
    return source.replace(ANCHOR, prime + ANCHOR, 1)


def edition_roots(primary, warm, size):
    case = f'matched-{size}'
    return [*(primary / 'matrix' / case / f'pair-01-{role}' for role in ROLES),
            *(warm / case / role / suffix for role in ROLES for suffix in ('prime', ''))]


def settings(pin):
    data = driver.check_control(pin).decode()
    reader = csv.DictReader(data.splitlines(keepends=True), delimiter='\t', strict=True)
    driver.require(reader.fieldnames == ['key', 'value'], 'Invalid cold worker settings')
    values = {}
    for row in reader:
        driver.require(set(row) == {'key', 'value'} and row['key'] not in values
                       and row['value'] is not None, 'Duplicate or unfinished cold worker settings')
        values[row['key']] = row['value']
    return values


def warm_inputs(runner, primary, size):
    name = f'matched-{size}'
    config = json.loads(driver.check_control(primary['matrix_config']))
    config['cases'] = [case for case in config['cases'] if case['id'] == name]
    case, = config['cases']
    reference_name = case['reference']
    config['references'] = {reference_name: config['references'][reference_name]}
    hashes, packages, builds = runner.resolve_configuration(config, 'diagnostic')
    assets = runner.asset_inventory(hashes)
    installed, toolchain = runner.read_provenance(config, packages, builds, assets)
    for key in ('matrix_pins', 'matrix_post_pins'):
        cold = json.loads(driver.check_control(primary[key]))
        driver.require(all(cold['assets'].get(path) == digest for path, digest in assets.items())
                       and cold['installed_packages'] == installed, 'Warm assets or packages differ from sealed cold run')
    identity = {'size': size, 'assets': assets, 'toolchain': toolchain,
                'packages': {role: {'library': config[role]['library'], 'commit': config[role]['commit'],
                                    'build_receipt': config[role]['build_receipt'],
                                    'tree_sha256': installed[role]['tree_sha256']} for role in ROLES},
                'cold_configuration': primary['matrix_config'], 'cold_pins': primary['matrix_pins'],
                'cold_post_pins': primary['matrix_post_pins']}
    return config, case, assets, toolchain, identity


def measure(runner, primary, root, primary_root, size, worker_source, prime_source):
    name = f'matched-{size}'
    scale = root / name
    scale.mkdir()
    # Pre-read and retain the selected pinned input/reference and matched
    # installation identity before starting either role's R process.
    config, _, _, _, identity = warm_inputs(runner, primary, size)
    runner.write_json(scale / 'inputs-before.json', identity)
    worker = scale / 'warm-worker.R'
    source = driver.control_bytes(worker_source).decode()
    prime = driver.control_bytes(prime_source).decode()
    with worker.open('x', encoding='utf-8', newline='\n') as stream:
        stream.write(render_worker(source, prime))
    worker_pin = driver.pin_control(worker)
    raw_failure = False
    for role in ROLES:
        dest = scale / role
        dest.mkdir()
        cold = primary['scale_controls'][size][role]
        values = settings(cold['settings.tsv'])
        cold_dir = primary_root / 'matrix' / name / f'pair-01-{role}'
        driver.require(values['mode'] == 'matrix' and values['profile'] == 'none'
                       and Path(values['output']) == cold_dir
                       and int(values['rows']) == driver.SCALES[size][0], 'Wrong cold comparator settings')
        values['output'] = str(dest)
        ini = dest / 'settings.tsv'
        runner.write_settings(ini, values)
        command = runner.r_command(config['r_executable']['path'], worker, 'run', str(ini))
        started = runner.utc_now()
        process = runner.run_process(command, dest)
        process.update(started_utc=started, finished_utc=runner.utc_now(), role=role, case=name,
                       application_state='same_verified_handle_after_complete_input_prime',
                       process_peak_rss_includes_prime=True, os_page_cache='uncontrolled', worker=worker_pin)
        runner.write_json(dest / 'process.json', process)
        raw_failure |= process['exit_code'] != 0
    return int(raw_failure)


def verify(runner, primary, root, primary_root, size, worker_source, prime_source):
    name = f'matched-{size}'
    scale = root / name
    opening = driver.pin_control(scale / 'inputs-before.json')
    config, case, assets, toolchain, identity = warm_inputs(runner, primary, size)
    runner.write_json(scale / 'inputs-after.json', identity)
    driver.require(json.loads(driver.check_control(opening)) == identity, 'Warm opening and closing identities differ')
    ref = config['references'][case['reference']]
    roles = {}
    for role in ROLES:
        dest = scale / role
        prefix = primary['scale_controls'][size][role]
        old = settings(prefix['settings.tsv'])
        driver.require(settings(driver.pin_control(dest / 'settings.tsv')) == {**old, 'output': str(dest)},
                       'Warm settings differ from cold comparator')
        old_dir = Path(old['output'])
        driver.require(old_dir == primary_root / 'matrix' / name / f'pair-01-{role}',
                       'Warm comparator directory differs')
        process = json.loads(driver.control_bytes(dest / 'process.json'))
        driver.require(process['exit_code'] == 0 and process['role'] == role and process['case'] == name
                       and process['process_peak_rss_includes_prime'] is True, 'Warm process failed or is incomplete')
        generated = driver.check_control(process['worker']).decode()
        driver.require(Path(process['worker']['path']) == scale / 'warm-worker.R'
                       and generated == render_worker(driver.control_bytes(worker_source).decode(),
                                                      driver.control_bytes(prime_source).decode()),
                       'Warm worker derivation differs')
        for filename in ('metrics.tsv', 'namespaces.tsv', 'process.json'):
            driver.check_control(prefix[filename])
        metrics = runner.tsv_pairs(dest / 'metrics.tsv', ('key', 'value'))
        cold_metrics = runner.tsv_pairs(old_dir / 'metrics.tsv', ('key', 'value'))
        driver.require(runner.verify_runtime(config, metrics, toolchain)
                       == runner.verify_runtime(config, cold_metrics, toolchain)
                       and runner.dependency_state(dest, config, role) == runner.dependency_state(old_dir, config, role),
                       'Warm runtime or namespace differed from cold comparator')
        prime = runner.tsv_pairs(dest / 'prime/metrics.tsv', ('key', 'value'))
        driver.require(all(prime[key] == str(case['rows']) for key in ('source_rows', 'normalized_rows'))
                       and prime['scope'] == 'complete identical input on the same verified handle',
                       'Warm prime was incomplete')
        roles[role] = {'process': process, 'warm_metrics': metrics, 'prime_metrics': prime,
                       'cold_comparator': prefix}

    # Read all six editions, compare all stable artifacts and scan all 19 fields.
    # Byte parity transfers the independently checked SNP oracle from cold pair-01.
    cold_oracle = primary['snp_oracles'][size]
    oracle = json.loads(driver.check_control(cold_oracle))
    driver.require(oracle['passed'] is True and oracle['all_rows_examined'] is True
                   and oracle['rows_expected'] == case['rows']
                   and oracle['versions']['mss7']['counts']['all_19_fields_correct_rows'] == case['rows']
                   and oracle['source_bytes_sha256'] == case['input']['sha256'], 'Wrong independent cold oracle')
    directories = edition_roots(primary_root, root, size)
    equality = runner.verify_editions(case, ref, directories, scale, assets, config.get('edition_encoding', 'utf-8'))
    recorded = json.loads(driver.check_control(primary['scale_controls'][size]['artifact_digests']))
    current = json.loads(driver.control_bytes(scale / 'artifact-digests.json'))
    driver.require(all(current[str(cold)] == recorded[str(cold)] for cold in directories[:2]),
                   'Warm comparator bytes differ from independently sealed primary evidence')
    runner.verify_manifests(config, scale, directories)
    equality['stable_rds_manifest_equality'] = True
    runner.write_json(scale / 'verification.json', equality)
    for observed in roles.values():
        observed['rss'] = driver.summarize_rss([observed['process']])
    roles['baseline']['rss'].pop('limit_status')
    print(driver.encode({'size': size, 'pairs': 1, 'roles': roles, 'verification': equality,
                         'inputs_before': opening, 'inputs_after': driver.pin_control(scale / 'inputs-after.json'),
                         'independent_oracle_transferred_by_byte_parity': cold_oracle,
                         'candidate_8GiB_status': roles['candidate']['rss']['limit_status']}).decode(), end='')
    return 0
