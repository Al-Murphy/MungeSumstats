#!/usr/bin/env python3
"""Optional bcftools normalization oracle; never used by production APIs."""
import pathlib, random, shutil, subprocess, tempfile
root = pathlib.Path(__file__).resolve().parents[2]
if not shutil.which('bcftools'):
    raise SystemExit('Optional oracle requires bcftools; runtime MSS does not.')
rng = random.Random(20261115)
with tempfile.TemporaryDirectory(prefix='mss-oracle-') as work:
    work = pathlib.Path(work)
    sequence = 'CAAAAAACACACACACGTGTGTGT' + ''.join(rng.choice('ACGT') for _ in range(200))
    fasta = work / 'tiny.fa'
    fasta.write_text('>1\n' + sequence + '\n')
    records = []
    for index in range(1500):
        pos = rng.randrange(1, len(sequence)-8)
        ref = sequence[pos-1:pos-1+rng.randrange(1,7)]
        if index % 3 == 0:
            alt = ref[:rng.randrange(1,len(ref)+1)]
        elif index % 3 == 1:
            alt = ref + ''.join(rng.choice('ACGT') for _ in range(rng.randrange(1,5)))
        else:
            alt = ''.join(rng.choice('ACGT') for _ in range(rng.randrange(1,7)))
        if ref != alt:
            records.append(('1',str(pos),ref,alt,f'v{index}'))
    vcf = work / 'source.vcf'
    vcf.write_text('##fileformat=VCFv4.2\n##contig=<ID=1,length='+str(len(sequence))+'>\n'
                   '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' +
                   ''.join(f'{chrom}\t{pos}\t{name}\t{ref}\t{alt}\t.\t.\t.\n' for chrom,pos,ref,alt,name in records))
    native = work / 'normalize'
    subprocess.run(['c++','-std=c++17','-O2',str(root/'tests/native/normalize.cpp'),'-o',str(native)],check=True)
    observed = subprocess.run([str(native),str(fasta)],input=''.join('\t'.join(row)+'\n' for row in records),text=True,capture_output=True,check=True)
    actual = {}
    for line in observed.stdout.splitlines():
        name,status,pos,ref,alt = line.split('\t')
        assert status == 'ok', (name,status)
        actual[name] = (pos,ref,alt)
    oracle = subprocess.run(['bcftools','norm','-f',str(fasta),'-Ov',str(vcf)],text=True,capture_output=True,check=True)
    expected = {v[2]:(v[1],v[3],v[4]) for line in oracle.stdout.splitlines() if not line.startswith('#') for v in [line.split('\t')]}
    assert actual == expected, [(name,actual[name],expected.get(name)) for name in actual if actual[name] != expected.get(name)][:10]
    print(f'{len(actual)} seeded SNP/MNV/indel normalization contrasts agree with bcftools')
    print(subprocess.run(['bcftools','--version'],capture_output=True,text=True,check=True).stdout.splitlines()[0])
