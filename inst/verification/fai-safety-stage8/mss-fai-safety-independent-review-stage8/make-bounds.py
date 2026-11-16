from pathlib import Path
import itertools,random,re,json

root=Path(__file__).resolve().parent
maximum=2**64-1
file_bytes=len(b'>x\n'+b'ACGT'*15+b'\n')
assert file_bytes==64
cases=[]
def add(name,values,expected=None):
    if expected is None:
        nums=[]
        for value in values:
            if not re.fullmatch(rb'[0-9]+',value) or int(value)>maximum:
                expected=False;break
            nums.append(int(value))
        else:
            size,offset,bases,width=nums
            if size==0:
                expected=bases==width==0 and offset<=file_bytes
            elif bases==0 or bases>size or width<bases or offset>=file_bytes:
                expected=False
            else:
                # Independent arbitrary-precision sum; no division-before-multiply implementation.
                last=offset+((size-1)//bases)*width+(size-1)%bases
                expected=last<file_bytes
    row=b'x\t'+b'\t'.join(values)+b'\n'
    cases.append((name,row,bool(expected)))

small=[0,1,2,3,4,8,16,64]
for i,values in enumerate(itertools.product(small,repeat=4)):
    add(f'small{i}',[str(v).encode() for v in values])
large=[0,1,2,64,2**63-1,2**63,maximum-1,maximum]
for i,values in enumerate(itertools.product(large,repeat=4)):
    add(f'extreme{i}',[str(v).encode() for v in values])
rng=random.Random(80171)
for i in range(3000):
    add(f'random{i}',[str(rng.randrange(0,85)).encode() for _ in range(4)])
for field in range(4):
    for i,value in enumerate([b'',b'-0',b'+1',b' 1',b'1 ',b'0x1',b'1e1',b'1\0',b'1.0',b'18446744073709551616',b'00018446744073709551616',b'\xff',b'\xc2\xb2',b'\xd9\xa1']):
        values=[b'1',b'3',b'1',b'2'];values[field]=value
        add(f'token{field}_{i}',values)
    for i,value in enumerate([b'00',b'0001',b'00000000000000000000000000000000000000000000000000000000000001']):
        values=[b'1',b'3',b'1',b'2'];values[field]=value
        add(f'leadingzero{field}_{i}',values)
for name,row in [
    ('blankrow',b'\n'),('embeddedCR',b'x\rname\t1\t3\t1\t2\n'),
    ('NULname',b'x\0\t1\t3\t1\t2\n'),('tabname',b'x\tname\t1\t3\t1\t2\n'),
    ('duplicate',b'x\t1\t3\t1\t2\nx\t1\t3\t1\t2\n'),
    ('blanktail',b'x\t1\t3\t1\t2\n\n'),('extraCR',b'x\t1\t3\t1\t2\r\r\n')]:
    cases.append((name,row,False))
for name,row in [('CRLF',b'x\t1\t3\t1\t2\r\n'),('unterminated',b'x\t1\t3\t1\t2'),('unused_max_stride',f'x\t1\t3\t1\t{maximum}\n'.encode()),('empty_eof',b'x\t0\t64\t0\t0\n')]:
    cases.append((name,row,True))
def expected_queries(row):
    size,offset,bases,width=map(int,row.rstrip(b'\r\n').split(b'\t')[1:])
    content=b'>x\n'+b'ACGT'*15+b'\n'
    result=[]
    for pos in [-1,0,1,2,63,64,2**63-1]:
        for length in [0,1,2,64,maximum]:
            # Independent unbounded integer request range and byte addresses.
            if pos<1 or pos-1+length>size:
                payload=b''
            else:
                payload=bytes(content[offset+(i//bases)*width+i%bases]
                              for i in range(pos-1,pos-1+length)).upper()
            result.append(payload.hex())
    return ';'.join(result)
(root/'bounds-cases.tsv').write_text(''.join(f'{name}\t{row.hex()}\t{int(expected)}\t{expected_queries(row) if expected else "-"}\n' for name,row,expected in cases))
summary={'cases':len(cases),'oracle_accept':sum(e for _,_,e in cases),'oracle_refuse':sum(not e for _,_,e in cases),'actual_fasta_bytes':file_bytes,'seed':80171,'oracle':'Python arbitrary-precision last-address formula and full ASCII uint64 syntax'}
(root/'bounds-input-summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary))
