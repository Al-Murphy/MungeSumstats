"""Sequential matched reference-preparation measurements, shared host/cache uncontrolled."""
import os
from pathlib import Path
import subprocess

BASE = Path('/private/tmp/gwas-harvest-stage2')
OUT = BASE / 'results/reference-index-capacity-v1'
COMMON = ':'.join(str(p) for p in [BASE / 'Rlib/reference-contigs-v1', BASE / 'Rlib',
    Path('/private/tmp/harvest-Rlib'), Path('/Users/bschilder/Library/R/arm64/4.5/library')])


def main():
    for label in ('baseline', 'candidate'):
        target = 'reference-index-' + label + '-matched'
        env = os.environ.copy()
        env['R_LIBS_USER'] = (str(BASE / 'Rlib/fasta-offsets-v1') + ':' if label == 'candidate' else '') + COMMON
        with (OUT / (label + '-matched.log')).open('x') as log:
            subprocess.run(['python3', str(BASE / 'results/full-reference/guard.py'),
                '--receipt', str(OUT / (label + '-matched-guard.json')), '--filesystem', '/private/tmp',
                '--', '/usr/bin/time', '-l', 'Rscript', '--vanilla',
                str(OUT / 'build_sample_matched.R'), target], env=env, stdout=log,
                stderr=subprocess.STDOUT, check=True)
        print(label + ' complete', flush=True)


if __name__ == '__main__':
    main()
