"""Strict archive/control/science primitives for installed-package CI.

Extracted from the independently accepted receiver. This module never imports
or executes tools from an artifact. A separate service/source receiver must
still establish terminal identity and current-checkout installation provenance.
"""
import csv
import gzip
import hashlib
import json
from pathlib import Path, PureWindowsPath, PurePosixPath
import re
import shutil
import stat
import struct
import sys
import unicodedata
import zipfile

PACKAGE_ROOT = Path(__file__).resolve().parents[3]
BENCH = PACKAGE_ROOT / "benchmarks/2026-11-19/performance-upgrade"
sys.path.insert(0, str(BENCH))
from verify import digest, stream_digest, verify_normalized

MIB = 1 << 20
MAX_MEMBERS = 32768
MAX_DIRECTORY = 16*MIB
EVIDENCE_LIMIT = 512_000_000
OUTSIDE_LIMITS = {name:MIB for name in ("preflight.json","job-start.json","q/deadlines.json",
    "q/UNSEALED.json","q/artifact-status.json","q/job-control-stop.json","q/qualify-stop.json","q/seal-stop.json")}
OUTSIDE_LIMITS.update({"q/"+mode+"."+suffix:8*MIB for mode in ("job-control","qualify","seal")
                       for suffix in ("stdout","stderr")})
OUTSIDE_LIMITS["q/artifact-inventory.json"] = MAX_DIRECTORY
OUTSIDE_LIMITS["q/upload-omissions.json"] = MIB
LINUX_OUTSIDE_LIMITS = {name:MIB for name in ("preflight.json", "job-start.json", "q/deadlines.json",
    "q/linux-stop.json", "q/artifact-status.json", "q/UNSEALED.json", "q/upload-omissions.json")}
LINUX_OUTSIDE_LIMITS["q/artifact-inventory.json"] = MAX_DIRECTORY
LINUX_SUCCESS_OUTSIDE = set(LINUX_OUTSIDE_LIMITS) - {"q/UNSEALED.json", "q/upload-omissions.json"}
SUCCESS_OUTSIDE = set(OUTSIDE_LIMITS)-{"q/upload-omissions.json","q/UNSEALED.json","q/qualify.stdout","q/qualify.stderr",
                                      "q/job-control.stdout","q/job-control.stderr"}


PARSER_GOLD = {
    "failure": (14, "e9bf8672bac2a2dee8500b9b073adde1bfebcfa6223a42833e119b6ec3d33fcc"),
    "width": (4, "37823f11f1e75ee997f458a9fab9deac48bc960bdeb670a0979df4f7ed7b000b"),
    "precedence": (16, "40677dc96d506be51f8ce07a64ecfef1b9a518adf92c55df3ab726fc0876dd31"),
}


def inventory_link(value, native_root):
    """Match the absolute producer inventory path without resolving it locally."""
    if re.match(r"^[A-Za-z]:[\\/]", str(native_root)):
        return canonical_work_path(value) == PureWindowsPath(native_root) / "artifact-inventory.json"
    require(isinstance(value, str) and value.startswith("/") and isinstance(native_root, str)
            and native_root.startswith("/"), "Invalid native Linux inventory path")
    for item in (value, native_root):
        require(not any(part in ("", ".", "..") or "\x00" in part for part in item.split("/")[1:]),
                "Noncanonical native Linux inventory path")
    return PurePosixPath(value) == PurePosixPath(native_root) / "artifact-inventory.json"


def plain_path(path, directory):
    info = path.lstat()
    mode = stat.S_ISDIR if directory else stat.S_ISREG
    if (not mode(info.st_mode) or getattr(info, "st_file_attributes", 0) & stat.FILE_ATTRIBUTE_REPARSE_POINT
            or not directory and info.st_nlink != 1 or path.resolve(strict=True) != path):
        raise ValueError(f"Expected a plain physical path: {path}")
    return info


def verify_suite_fixture(path, spec):
    plain_path(path, False)
    value = digest(path)
    if any(value[key] != spec[key] for key in ("bytes", "sha256")):
        raise ValueError("Installed or archive ALS VCF fixture differs")
    return dict(path=str(path), **value)


def read_table(path, fields, *, max_rows=None):
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.reader(stream, delimiter="\t", strict=True)
        if next(reader, None) != fields:
            raise ValueError(f"Unexpected columns: {path}")
        rows = []
        for values in reader:
            if max_rows is not None and len(rows)>=max_rows:
                raise ValueError("Too many TSV rows: "+str(path))
            if len(values) != len(fields):
                raise ValueError(f"Unexpected row width: {path}")
            rows.append(dict(zip(fields, values)))
        return rows


def verify_parser_gold(path, mode):
    count, expected = PARSER_GOLD[mode]
    with path.open(encoding="utf-8", newline="") as stream:
        contents = list(csv.reader(stream, delimiter="\t", strict=True))
    encoded = json.dumps(contents, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    if len(contents) != count+1 or hashlib.sha256(encoded).hexdigest() != expected:
        raise ValueError(f"Unexpected public parser fields for {mode}")
    return dict(records=count, fields_sha256=expected, file=digest(path))


def verify_lookup_pressure(path):
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) != 6:
        raise ValueError("Incomplete counter-enabled lookup pressure transcript")
    result = []
    for index, fmt in ((0,2),(3,3)):
        warm = re.fullmatch(r"warm negative allocations (\d+) for 1024 queries", lines[index+1])
        counters = re.fullmatch(r"shuffled unique 3000 ordinary probes (\d+) loaded (\d+) fence hits (\d+) "
                    r"key page requests (\d+) reads (\d+) offset page requests (\d+) reads (\d+) rank bytes (\d+)",
                    lines[index+2])
        if lines[index] != f"format {fmt}" or not warm or not counters:
            raise ValueError("Missing lookup pressure counters")
        total, loads, hits, key_requests, key_reads, offset_requests, offset_reads, rank = map(int,counters.groups())
        # The pinned lower bounds for target ranks 401..3400 in 3403 entries
        # take 35,389 comparisons in any order. Each unique hit also visits
        # its own entry and the following nonmatch (both have a successor).
        entry_visits = loads + 2*3000
        if (int(warm[1]) > 64 or total != 35389 or loads <= 0 or loads+hits != total or 2*loads >= total
                or not 0 < rank <= 1024*1024 or min(key_requests,offset_requests) < entry_visits
                or key_reads > key_requests or offset_reads > offset_requests):
            raise ValueError("Invalid lookup counters or warm negative allocations")
        result.append(dict(format=fmt, allocations=int(warm[1]), ordinary_probes=total, loaded=loads,
                           fence_hits=hits, key_requests=key_requests, key_reads=key_reads,
                           offset_requests=offset_requests, offset_reads=offset_reads, rank_bytes=rank))
    return result


def verify_typed_boundaries(path):
    data = path.read_bytes().replace(b"\r\n", b"\n")
    marker = b"typed boundary newline at 16383\n"
    # Frozen final section checks all newline, malformed-offset and split-UTF8
    # queries, budgets and returned candidate fields. Preserve the complete raw log.
    expected = "a940b9548576b4cbc195d289fdead691243577348740a6544589ff9de11693f6"
    if data.count(marker) != 1 or hashlib.sha256(data[data.index(marker):]).hexdigest() != expected:
        raise ValueError("Typed key boundary evidence differs")
    return expected


def canonical_work_path(value):
    # Canonical native work is on the runner's local drive. Check text before
    # PureWindowsPath can discard dot components or redundant separators.
    if not isinstance(value, str) or not re.match(r"^[A-Za-z]:[\\/]", value):
        raise ValueError("Expected a fully qualified local Windows path")
    parts = re.split(r"[\\/]", value[3:])
    for part in parts:
        if (not part or part.endswith((".", " ")) or re.search(r'[<>:"|?*\x00-\x1f]', part)
                or re.fullmatch(r"(?:CON|PRN|AUX|NUL|CONIN\$|CONOUT\$|(?:COM|LPT)[1-9¹²³])",
                                part.split(".",1)[0].rstrip(" "), flags=re.IGNORECASE)):
            raise ValueError("Noncanonical Windows path component")
    return PureWindowsPath(value)


def verify_vcf_locations(row, home, purpose):
    expected_tmp = canonical_work_path(home).parent/"tmp"
    tmp, tempdir, source, bgzf, tabix = (canonical_work_path(row[key]) for key in
                                        ("tmp","tempdir","input","bgzf","tabix"))
    if (tmp != expected_tmp or tempdir.parent != tmp or source.parent.parent != tempdir
            or not tempdir.name.startswith("Rtmp") or purpose not in ("imputation","roundtrip")
            or not re.fullmatch("alsvcf-"+purpose+"-[a-zA-Z0-9]+", source.parent.name) or source.name != "ALSvcf.vcf"
            or bgzf != source.with_name("ALSvcf.vcf.bgz") or tabix != source.with_name("ALSvcf.vcf.bgz.tbi")):
        raise ValueError("Native VCF round-trip files are outside the owned role temp directory")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, "Duplicate JSON key: "+key)
        result[key] = value
    return result


def json_value(data):
    def invalid(value):
        raise ValueError("Invalid JSON number: "+value)
    return json.loads(data.decode("utf-8-sig"),object_pairs_hook=unique_object,parse_constant=invalid)


def document(path):
    require(path.stat().st_size <= MAX_DIRECTORY, "Oversized receipt: "+path.name)
    return json_value(path.read_bytes())


def match(value, expected, label):
    # Python equality equates false with 0 and true with 1.
    require(json.dumps(value,sort_keys=True,allow_nan=False)==json.dumps(expected,sort_keys=True,allow_nan=False),
            "Inconsistent "+label)


def relative_name(name):
    require(isinstance(name,str) and len(name.encode("utf-8"))<=1024 and "\\" not in name
            and len(name.split("/"))<=64, "Invalid ZIP or inventory path")
    canonical_work_path("C:/"+name)
    return name


def directory_bounds(stream, size):
    # ZipFile first loads the whole central directory. Inspect the fixed footer
    # and each bounded central record before allowing that allocation. Actions
    # uploads are single-disk, have no ZIP comment, and may use a ZIP64 footer.
    require(size>=22,"Truncated ZIP footer")
    stream.seek(size-22)
    magic, disk, start, count1,count,length,offset,comment = struct.unpack("<4s4H2LH",stream.read(22))
    require(magic==b"PK\x05\x06" and disk==start==comment==0 and count==count1,"Unsupported ZIP footer")
    end=size-22
    if end>=20:
        stream.seek(end-20)
        locator=stream.read(20)
        if locator.startswith(b"PK\x06\x07"):
            _, disk64, at, disks=struct.unpack("<4sLQL",locator)
            require(disk64==0 and disks==1 and at+56==end-20,"Invalid ZIP64 locator")
            stream.seek(at)
            magic64,recordlen,_,_,disk,start,n1,n,ln,off=struct.unpack("<4sQ2H2L4Q",stream.read(56))
            require(magic64==b"PK\x06\x06" and recordlen==44 and disk==start==0 and n1==n
                    and count in (n,65535) and length in (ln,0xffffffff) and offset in (off,0xffffffff),
                    "Inconsistent ZIP64 footer")
            count,length,offset,end=n,ln,off,at
    require(0<count<=MAX_MEMBERS and length<=MAX_DIRECTORY and offset+length==end,"Unbounded ZIP directory")
    stream.seek(offset)
    for _ in range(count):
        header=stream.read(46)
        require(len(header)==46 and header[:4]==b"PK\x01\x02","Incorrect ZIP member count")
        name,extra,comment=struct.unpack_from("<3H",header,28)
        require(0<name<=1024 and extra<=1024 and comment<=1024,"Oversized ZIP member metadata")
        stream.seek(name+extra+comment,1)
    require(stream.tell()==end,"Incorrect ZIP directory length or count")


    return offset


def local_headers(stream, archive, central_offset):
    # A central directory may point past an unlisted member, or describe a
    # different local filename. Traverse every physical record, including the
    # streaming descriptors emitted when the uploader cannot seek backwards.
    members=sorted(archive.infolist(),key=lambda i:i.header_offset)
    boundaries=[i.header_offset for i in members[1:]]+[central_offset]
    require(members[0].header_offset==0,"ZIP prefix or unselected local member")
    for item,end in zip(members,boundaries):
        stream.seek(item.header_offset)
        header=stream.read(30)
        require(len(header)==30 and header[:4]==b"PK\x03\x04","Missing local ZIP header")
        _,_,flags,method,_,_,crc,compressed,expanded,n,e=struct.unpack("<4s5H3L2H",header)
        require(flags==item.flag_bits and flags & ~0x80e==0 and method==item.compress_type and 0<n<=1024 and e<=1024,
                "Inconsistent or oversized local ZIP metadata")
        require(stream.read(n)==item.orig_filename.encode("utf-8" if flags&0x800 else "cp437"),"Local ZIP path differs")
        extra=stream.read(e)
        zip64=None
        while extra:
            require(len(extra)>=4,"Incomplete local ZIP extra field")
            key,size=struct.unpack_from("<HH",extra)
            require(size<=len(extra)-4,"Invalid local ZIP extra field")
            if key==1:
                require(zip64 is None,"Duplicate local ZIP64 record")
                zip64=extra[4:4+size]
            extra=extra[4+size:]
        sizes=[expanded,compressed]
        for i,value in enumerate(sizes):
            if value==0xffffffff:
                require(zip64 is not None and len(zip64)>=8,"Missing local ZIP64 length")
                sizes[i],=struct.unpack_from("<Q",zip64)
                zip64=zip64[8:]
        descriptor=bool(flags&8)
        for value,want in zip([crc,*sizes],[item.CRC,item.file_size,item.compress_size]):
            require(value==want or descriptor and value==0,"Local and central ZIP sizes or CRC differ")
        after=stream.tell()+item.compress_size
        if descriptor:
            length=end-after
            require(length in (12,16,20,24),"Unexpected ZIP descriptor or hidden local member")
            stream.seek(after)
            data=stream.read(length)
            if length in (16,24):
                require(data[:4]==b"PK\x07\x08","Invalid ZIP descriptor signature")
                data=data[4:]
            require(struct.unpack("<L2L" if len(data)==12 else "<L2Q",data)==(item.CRC,item.compress_size,item.file_size),
                    "ZIP descriptor contradicts central sizes")
        else:
            require(after==end,"Unexpected data or unselected local ZIP member")


def member_names(archive, *, lane="windows"):
    require(lane in {"windows", "linux-gold"}, "Trusted lane schema required")
    limits = OUTSIDE_LIMITS if lane == "windows" else LINUX_OUTSIDE_LIMITS
    files, dirs, spellings = {}, set(), {}
    total=0
    for item in archive.infolist():
        name=relative_name(item.orig_filename.removesuffix("/") if item.is_dir() else item.orig_filename)
        require(name not in dirs and name not in files and item.filename==item.orig_filename,
                "Duplicate or aliased ZIP member: "+name)
        kind=stat.S_IFMT(item.external_attr>>16)
        require(kind in (0,stat.S_IFDIR if item.is_dir() else stat.S_IFREG)
                and not item.external_attr & (0x400|0x08) and not item.flag_bits & 1
                and item.compress_type in (zipfile.ZIP_STORED,zipfile.ZIP_DEFLATED)
                and (not item.external_attr & 0x10 or item.is_dir()), "Special ZIP member: "+name)
        if item.is_dir():
            require(item.file_size==0,"Nonempty ZIP directory")
            dirs.add(name)
        else:
            limit=EVIDENCE_LIMIT if name.startswith("q/evidence/") else limits.get(name,-1)
            require(item.file_size<=limit, "Unexpected or oversized ZIP member: "+name)
            files[name]=item
            total+=item.file_size
        register_spelling(name,spellings)
    require(total<=EVIDENCE_LIMIT+sum(limits.values()),"ZIP logical size exceeds bound")
    return files,dirs,spellings


def register_spelling(name, spellings):
    parts=name.split("/")
    for end in range(1,len(parts)+1):
        prefix="/".join(parts[:end])
        folded=unicodedata.normalize("NFD",prefix).casefold()
        require(spellings.setdefault(folded,prefix)==prefix, "Windows or Unicode path alias: "+prefix)


def unpack_sealed(zip_path, target, size, sha256, native_root, *, lane="windows"):
    require(lane in {"windows", "linux-gold"}, "Trusted lane schema required")
    limits = OUTSIDE_LIMITS if lane == "windows" else LINUX_OUTSIDE_LIMITS
    success = SUCCESS_OUTSIDE if lane == "windows" else LINUX_SUCCESS_OUTSIDE
    """Validate member selection before creating the fresh extraction root."""
    plain_path(zip_path,False)
    require(type(size) is int and zip_path.stat().st_size==size<=EVIDENCE_LIMIT+96*MIB,"ZIP size differs")
    with zip_path.open("rb") as stream:
        require(stream_digest(stream)["sha256"]==sha256,"Service ZIP digest differs")
        central_offset=directory_bounds(stream,size)
        with zipfile.ZipFile(stream) as archive:
            files,dirs,spellings=member_names(archive, lane=lane)
            local_headers(stream,archive,central_offset)
            outside=set(files).intersection(limits)
            require(outside==success,"Unsealed, missing or unexpected success receipt")
            status=json_value(archive.read(files["q/artifact-status.json"]))
            for field,expected in (("qualified",True),("source_sealed",True),("native_qualification",True),
                                   ("receiver_verified",False)):
                match(status[field],expected,"artifact status "+field)
            raw=archive.read(files["q/artifact-inventory.json"])
            require(inventory_link(status["inventory"]["path"],native_root)
                    and hashlib.sha256(raw).hexdigest()==status["inventory"]["sha256"],"Sealer inventory pin differs")
            inventory=json_value(raw)
            sealed=inventory["files"]
            directories=inventory["directories"]
            require(isinstance(sealed,dict) and isinstance(directories,list)
                    and len(directories)+len(sealed)<=MAX_MEMBERS,"Oversized inventory")
            listed_dirs=set(directories)
            require(len(listed_dirs)==len(directories) and not listed_dirs.intersection(sealed),"Duplicate inventory path")
            for name in [*directories,*sealed]:
                relative_name(name)
                register_spelling("q/evidence/"+name,spellings)
                require("/" not in name or name.rsplit("/",1)[0] in listed_dirs,"Missing inventory parent")
            prefix="q/evidence/"
            allowed_dirs={"q","q/evidence"}|{prefix+n for n in listed_dirs}
            require(dirs<=allowed_dirs and set(files)-outside=={prefix+n for n in sealed},"ZIP and sealed path sets differ")
            total=0
            for name,entry in sealed.items():
                require(set(entry)=={"bytes","md5","sha256"} and type(entry["bytes"]) is int
                        and files[prefix+name].file_size==entry["bytes"] and re.fullmatch("[0-9a-f]{32}",entry["md5"])
                        and re.fullmatch("[0-9a-f]{64}",entry["sha256"]),"Invalid sealed size or digest: "+name)
                total+=entry["bytes"]
            require(type(inventory["bytes"]) is int and inventory["bytes"]==total<=EVIDENCE_LIMIT,"Evidence byte total differs")
            plain_path(target.parent,True)
            target.mkdir()
            for name in sorted(allowed_dirs,key=lambda s:s.count("/")):
                (target/name).mkdir()
            for name,item in files.items():
                dest=target/name
                with archive.open(item) as inp,dest.open("xb") as out:
                    shutil.copyfileobj(inp,out,length=MIB)
                if name.startswith(prefix):
                    match(digest(dest),sealed[name[len(prefix):]],"sealed digest "+name)
            return inventory


def check_job_control(row, folder):
    start=row["controller_started_monotonic"]
    require(row["mode"]=="job-control" and row["assigned"] is True and row["drained"] is True
            and row["control_passed"] is True and row["reason"]=="job_control_cutoff"
            and row["active_before"]>=2 and row["active_after"]==0 and row["exit_code"]==99 and start>=0,
            "Native Job control or drain failed")
    ends={name:start+sec for name,sec in (("readiness",5),("compute",6),("drain",10),("outer",12))}
    require(row["fixture_deadlines_monotonic"]==ends,"Job control clocks changed")
    controller,writer=(document(folder/name) for name in ("controller.json","writer.json"))
    require(controller["pid"]==row["pid"] and writer["parent"]==row["pid"] and
            len({controller["parent"],controller["pid"],writer["pid"]})==3 and
            all(type(v) is int and v>0 for v in (row["pid"],controller["parent"],writer["pid"])),
            "Invalid native controller or writer identity")
    require(start<=controller["started_monotonic"]<=writer["started_monotonic"] and
            all(p["natural_exit_monotonic"]>=p["started_monotonic"]+16 for p in (controller,writer)),
            "Invalid natural-exit bound")
    match(controller["attempt"],row["attempt"],"controller attempt")
    match(writer["attempt"],row["attempt"],"writer attempt")
    a,b=row["readiness_observations"]
    trigger=row["deadline_trigger"]
    drain=row["drain_elapsed_seconds"]
    first,last=row["post_drain_observations"]
    require(0<=writer["started_monotonic"]-start<=a["elapsed_seconds"]<b["elapsed_seconds"]<=5
            and trigger["phase"]=="compute" and 6<=trigger["elapsed_seconds"]<=drain<=10
            and drain<=first["elapsed_seconds"]<last["elapsed_seconds"]<=row["elapsed_seconds"]<=12,
            "Job control readiness, cutoff, drain or observation clocks differ")
    final=digest(folder/"progress.bin")
    require(0<a["bytes"]<b["bytes"]<=final["bytes"]<=4096
            and all(r["active_members"]>=2 and r["controller_pid"]==row["pid"] and r["writer_pid"]==writer["pid"]
                    for r in (a,b))
            and all(r["bytes"]==final["bytes"] and r["sha256"]==final["sha256"] for r in (first,last)),
            "Invalid growing or post-drain Job writer evidence")


def verify_job_commit_admission(evidence):
    import windows_job as w
    path=evidence/"job-commit-admission.jsonl"
    plain_path(path,False)
    require(0<path.stat().st_size<=16384,"Missing or oversized native commitment receipts")
    rows=[json.loads(line) for line in path.read_text(encoding="ascii").splitlines()]
    require(len(rows)==2 and [r["phase"] for r in rows]==["before-worker","before-release"],
            "Incomplete or reordered native commitment admission")
    for row in rows:
        expected=w.commitment_record(row["system_commit_limit_bytes"],row["system_commit_used_bytes"],
            row["job_current_before_bytes"],row["job_current_after_bytes"],row["active_members"],row["phase"])
        require(expected["decision"]=="admitted" and type(row["observer_pid"]) is int and row["observer_pid"]>0,
                "Rejected native commitment admission or missing observer")
        match(row,dict(**expected,anchor=row["anchor"],observer_pid=row["observer_pid"]),"native commitment arithmetic and policy")
    require(rows[0]["observer_pid"]==rows[1]["observer_pid"],"Commitment admission observers differ")
    return rows


def tiny_bgzf(path):
    # The reviewed controls contain at most 47,124 decoded bytes. Track each
    # block's physical offset so tabix virtual offsets can be checked.
    import zlib
    require(path.stat().st_size <= 65536, "Oversized VCF control " + path.name)
    raw = path.read_bytes()
    starts, result, pos = {}, bytearray(), 0
    while pos < len(raw):
        header = raw[pos:pos+18]
        require(len(header) == 18 and header[:4] == b"\x1f\x8b\x08\x04"
                and header[10:16] == b"\x06\x00BC\x02\x00", "Invalid control BGZF header")
        length = int.from_bytes(header[16:18], "little") + 1
        require(26 <= length <= len(raw)-pos, "Invalid control BGZF block length")
        decoder = zlib.decompressobj(31)
        chunk = decoder.decompress(raw[pos:pos+length], 65537-len(result))
        require(len(result)+len(chunk) <= 65536 and decoder.eof and not decoder.unused_data,
                "Invalid or oversized decoded control BGZF")
        starts[pos] = (len(result), len(chunk))
        result.extend(chunk)
        pos += length
        if not chunk:
            require(pos == len(raw), "Unexpected empty control BGZF block")
    require(starts and not chunk, "Missing control BGZF terminator")
    starts[len(raw)] = (len(result), 0)  # htslib can report the physical EOF offset.
    return bytes(result), starts


def check_vcf_index(path, blocks, records):
    data, _ = tiny_bgzf(path)
    require(data.startswith(b"TBI\x01"), "Invalid control tabix magic")
    pos = 4
    def take(fmt):
        nonlocal pos
        size = struct.calcsize(fmt)
        require(pos+size <= len(data), "Truncated control tabix index")
        values = struct.unpack_from(fmt, data, pos)
        pos += size
        return values
    def offset(value):
        block, inside = divmod(value, 65536)
        require(block in blocks and inside <= blocks[block][1], "Invalid data offset in control tabix index")
        return blocks[block][0] + inside
    require(take("<8i") == (1, 2, 1, 2, 0, 35, 0, 2) and take("<2s") == (b"1\0",),
            "Unexpected sequence or coordinate system in control tabix index")
    bins, = take("<i")
    require(1 <= bins <= min(37451, len(records)+1), "Unexpected control tabix bin count")
    coverage = {}
    for _ in range(bins):
        bin_id, chunks = take("<Ii")
        require(bin_id not in coverage and 0 <= bin_id <= 37450 and 1 <= chunks <= max(2,len(records)),
                "Invalid control tabix chunks")
        coverage[bin_id] = []
        for i in range(chunks):
            start, end = take("<QQ")
            if bin_id == 37450 and i == 1:
                require(chunks == 2 and (start,end) == (len(records),0), "Wrong tabix mapped count")
            else:
                coverage[bin_id].append((offset(start),offset(end)))
    intervals, = take("<i")
    require(intervals == 1+((records[-1][0]-1) >> 14), "Unexpected control tabix linear index")
    first = [offset(take("<Q")[0]) for _ in range(intervals)]
    for bp,start,end in records:
        x = bp-1
        parents = (0, 1+(x>>26), 9+(x>>23), 73+(x>>20), 585+(x>>17), 4681+(x>>14))
        require(first[x>>14] <= start and any(a<=start and end<=b for parent in parents
                for a,b in coverage.get(parent,())), "Control tabix index does not cover decoded records")
    require(data[pos:] in (b"", bytes(8)), "Unexpected control tabix trailer")


def verify_als_triplets(root, home, fixture, *, recorded, retained_root=None):
    purposes = {"imputation", "roundtrip"}
    require({p.name for p in root.iterdir()} == purposes, "Incomplete ALS fixture purpose roster")
    files = {"input":"ALSvcf.vcf", "bgzf":"ALSvcf.vcf.bgz", "tabix":"ALSvcf.vcf.bgz.tbi",
             "query":"indexed-records.txt"}
    observed = {}
    processes = set()
    for purpose in sorted(purposes):
        area = root/purpose
        wanted = set(files.values()) | {"locations.tsv"}
        if retained_root is not None:
            wanted.add("retained-locations.tsv")
        if recorded:
            wanted.add("files.json")
        require({p.name for p in area.iterdir()} == wanted and
                all(p.is_file() and p.stat().st_size <= 65536 for p in area.iterdir()),
                "Incomplete or oversized ALS fixture evidence")
        row, = read_table(area/"locations.tsv", "purpose tmp tempdir input bgzf tabix".split(), max_rows=1)
        require(row["purpose"] == purpose, "Wrong ALS fixture purpose")
        verify_vcf_locations(row, str(home), purpose)
        if retained_root is not None:
            fields = "purpose retained_directory query_input query_bgzf query_tabix original_absent reconstructed_original_paths query_records empty_range_records".split()
            kept, = read_table(area/"retained-locations.tsv", fields, max_rows=1)
            original = canonical_work_path(row["input"])
            retained = canonical_work_path(str(retained_root))/original.parent.name
            require(kept["purpose"] == purpose and canonical_work_path(kept["retained_directory"]) == retained
                    and kept["original_absent"] == kept["reconstructed_original_paths"] == "TRUE"
                    and kept["query_records"] == "415" and kept["empty_range_records"] == "0"
                    and canonical_work_path(kept["query_input"]) == retained/"ALSvcf.vcf"
                    and canonical_work_path(kept["query_bgzf"]) == retained/"ALSvcf.vcf.bgz"
                    and canonical_work_path(kept["query_tabix"]) == retained/"ALSvcf.vcf.bgz.tbi",
                    "Retained ALS query paths/mapping differ from byte-verified helper copies")
        processes.add(canonical_work_path(row["tempdir"]))
        verify_suite_fixture(area/files["input"], fixture)
        data, blocks = tiny_bgzf(area/files["bgzf"])
        require(data == (area/files["input"]).read_bytes(), "Decoded ALS BGZF differs from frozen input")
        records, offsets, start = [], [], 0
        for line in data.splitlines(keepends=True):
            if not line.startswith(b"#"):
                chrom, bp, _ = line.split(b"\t",2)
                require(chrom == b"1", "Unexpected chromosome in ALS fixture")
                records.append(line)
                offsets.append((int(bp),start,start+len(line)))
            start += len(line)
        require(len(records) == 415 and (area/files["query"]).read_bytes() == b"".join(records),
                "Indexed ALS observations differ from frozen records")
        check_vcf_index(area/files["tabix"], blocks, offsets)
        observed[purpose] = {key:digest(area/name) for key,name in files.items()}
        if recorded:
            match(document(area/"files.json"), observed[purpose], "ALS file pins for "+purpose)
    require(len(processes) == 1, "ALS fixture copies came from different suite processes")
    return observed


def verify_vcf_error_provenance(err, library, commit, temp_root, native_output):
    """Native Windows provenance independent of the receiver's operating system."""
    library = canonical_work_path(str(library))
    temp_root = canonical_work_path(str(temp_root))
    native_output = canonical_work_path(str(native_output))
    tmp = canonical_work_path(err["temp"])
    require(err["commit"] == commit and canonical_work_path(err["library"]) == library
            and tmp.parent == temp_root and tmp.name.startswith("Rtmp")
            and err["class"] == "simpleError,error,condition" and err["owned_stages"] == "0",
            "Wrong control error provenance or cleanup")
    dest = native_output/"blocked"
    messages = [err[k] for k in ("direct", "wrapped")]
    inputs = []
    for msg in messages:
        hits = set(re.findall(r"\n  file: ([^\n\r]+)", msg))
        require(len(hits) == 1 and "opening 'dest': " in msg
                and "\n  dest: " in msg, "Missing native bgzip error")
        inputs.append(hits.pop())
    raw = canonical_work_path(inputs[1])
    require(canonical_work_path(inputs[0]) == dest/"ss.vcf"
            and raw.parent.parent == tmp and raw.parent.name.startswith("mss-indexed-vcf-") and raw.name == "sorted.vcf"
            and messages[0].replace("\n  file: "+inputs[0],"\n  file: [input]") ==
                messages[1].replace("\n  file: "+inputs[1],"\n  file: [input]")
            and ("\n  dest: "+(dest/"ss.vcf.bgz").as_posix()).lower() in messages[0].replace("\\","/").lower(),
            "Changed native bgzip error or control paths")


def verify_vcf_control(out, library, commit, temp_root, *, native_output):
    sentinel = b"caller-owned uncompressed VCF sibling\n"
    names = {"checks.tsv", "error.tsv", "blocked"}
    core = "1:101:A:G,1:202:C:T"
    expected = []
    for path in out.rglob("*"):
        if path.is_file():
            require(path.stat().st_size <= 65536, "Oversized VCF control evidence")
    for mode in ("normal", "opengwas"):
        for suffix in ("vcf", "vcf.gz", "vcf.bgz"):
            for indexed in (False, True):
                label = f"{mode}-{suffix}-{int(indexed)}"
                names.add(label)
                folder = out/label
                output = "ss.vcf.bgz" if indexed else "ss."+suffix
                path = folder/output
                actual = {p.name for p in folder.iterdir()}
                require(actual == ({output, "ss.vcf", output+".tbi"} if indexed else {output}),
                        "Unexpected control VCF output roster")
                require(path.stat().st_size <= 65536, "Oversized control VCF")
                if output.endswith((".gz", ".bgz")):
                    data, blocks = tiny_bgzf(path)
                else:
                    data = path.read_bytes()
                header, records, ranges, start = [], [], [], 0
                for line in data.splitlines(keepends=True):
                    text = line.decode("utf-8").rstrip("\r\n")
                    if text.startswith("#"):
                        header.append(text)
                    else:
                        records.append(text.split("\t"))
                        ranges.append((start, start+len(line)))
                    start += len(line)
                require(header and header[-1].split("\t") == [
                    "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","GWAS"],
                    "Unexpected control VCF sample header")
                require(len(records) == 2 and all(len(row) == 10 for row in records),
                        "Unexpected control VCF record count")
                for row, bp, ref, alt, rs, p, beta, se, lp in zip(records, ("101","202"), "AC", "GT",
                        ("rs101","rs202"), (.01,.001), (.25,.5), (.125,.25), (2,3)):
                    require([row[k] for k in (0,1,3,4)] == ["1",bp,ref,alt], "Wrong control VCF alleles")
                    keys, values = row[8].split(":"), row[9].split(":")
                    require(len(keys) == len(set(keys)) == len(values), "Ambiguous control VCF sample fields")
                    fields = dict(zip(keys,values))
                    if mode == "opengwas":
                        require("RSID="+rs in row[7].split(";") and float(fields["ES"]) == beta
                                and float(fields["LP"]) == lp, "Wrong OpenGWAS control VCF data")
                    else:
                        require(fields["SNP"] == rs and float(fields["BETA"]) == beta
                                and abs(float(fields["P"]) - p) <= p*1e-6, "Wrong ordinary control VCF data")
                    require(float(fields["SE"]) == se, "Wrong control VCF standard error")
                if indexed:
                    require((folder/"ss.vcf").read_bytes() == sentinel, "Changed indexed VCF sibling")
                    check_vcf_index(folder/(output+".tbi"), blocks, [(int(row[1]),*span) for row,span in zip(records,ranges)])
                expected.append(dict(case=label,output=output,decoded=core,query=core+",empty" if indexed else "",
                                     owned_stages="0"))
    match(read_table(out/"checks.tsv", list(expected[0]), max_rows=12), expected, "hosted VCF matrix observations")
    require({p.name for p in out.iterdir()} == names, "Unexpected control evidence roster")
    blocked = out/"blocked"
    require({p.name for p in blocked.iterdir()} == {"ss.vcf","ss.vcf.bgz"}
            and {p.name for p in (blocked/"ss.vcf.bgz").iterdir()} == {"marker"}
            and (blocked/"ss.vcf").read_bytes() == (blocked/"ss.vcf.bgz/marker").read_bytes() == sentinel,
            "Changed blocked-destination VCF sentinel")
    err, = read_table(out/"error.tsv","commit library temp direct wrapped class owned_stages".split(),max_rows=1)
    verify_vcf_error_provenance(err, library, commit, temp_root, native_output)


def read_r_object(path,version):
    """Read the bounded XDR subset emitted by the frozen native R helpers."""
    with gzip.open(path,"rb") as stream:
        data=stream.read(MIB+1)
    require(len(data)<=MIB,"Oversized deferral result")
    pos=0
    symbols=[]
    nodes=0
    def take(n):
        nonlocal pos
        require(0<=n<=len(data)-pos,"Truncated R serialization")
        start=pos;pos+=n
        return data[start:pos]
    def integer():
        return struct.unpack(">i",take(4))[0]
    def vector_length():
        n=integer()
        require(0<=n<=65536,"Invalid R vector length")
        return n
    require(take(2)==b"X\n" and integer()==version and integer()==0x40503
            and integer()=={2:0x20300,3:0x30500}[version],"Unexpected R serialization header")
    encoding=take(vector_length()).decode("ascii") if version==3 else "utf-8"
    def node(depth=0):
        nonlocal nodes
        nodes+=1
        require(depth<128 and nodes<65536,"Excessive R object nesting or count")
        flags=integer();kind=flags&255
        def child():return node(depth+1)
        if kind==254:return None
        if kind==255:
            ref=flags>>8 or integer()
            require(1<=ref<=len(symbols),"Invalid R symbol reference")
            return symbols[ref-1]
        if kind==1:
            value=child()
            require(isinstance(value,str),"Invalid R symbol name")
            symbols.append(value)
            return value
        if kind==2:
            require(flags&1024 and not flags&512,"Expected a tagged attribute pairlist")
            key=child();value=child();tail=child() or {}
            require(isinstance(key,str) and key not in tail,"Repeated R attribute")
            return {key:value,**tail}
        if kind==9:
            n=integer()
            if n==-1:return None
            raw=take(n)
            return raw.decode("utf-8" if flags&(8<<12) else "latin-1" if flags&(4<<12) else encoding)
        if kind==238:
            info=child();state=child();attrs=child()
            return ["altrep",info,state,attrs]
        require(kind in (10,13,14,16,19),"Unsupported R result type "+str(kind))
        n=vector_length()
        if kind in (10,13):values=list(struct.unpack(">"+str(n)+"i",take(4*n)))
        elif kind==14:values=[take(8).hex() for _ in range(n)]
        else:values=[child() for _ in range(n)]
        attrs=child() if flags&512 else {}
        return [kind,values,attrs]
    result=node()
    require(pos==len(data),"Trailing bytes in R serialization")
    return result


def verify_reference_admission(path, capacity):
    """Replay the bounded available-RAM wait, then bind its final admission."""
    from datetime import datetime, timedelta
    import math
    import windows_job as w
    plain_path(path,False)
    require(path.stat().st_size<=65536,"Oversized six-reference admission log")
    with path.open("rb") as stream:
        data=stream.read(65537)
    require(len(data)<=65536 and data.endswith(b"\n"),"Incomplete six-reference admission log")
    lines=data.splitlines()
    require(0<len(lines)<=25,"Invalid six-reference admission sample count")
    fields={"utc","growth_bytes","elapsed_seconds","total_ram_bytes","available_ram_bytes","volumes","decision"}
    previous=-1
    for index,line in enumerate(lines):
        row=json_value(line)
        require(isinstance(row,dict) and row.get("decision") in ("wait","admitted"),"Invalid admission decision")
        require(set(row)==fields|({"reason"} if row["decision"]=="wait" else set()),"Invalid admission sample fields")
        elapsed=row["elapsed_seconds"]
        require(type(elapsed) in (int,float) and math.isfinite(elapsed) and 0<=elapsed<120 and elapsed>previous,
                "Invalid admission elapsed time")
        previous=elapsed
        require(isinstance(row["utc"],str),"Invalid admission UTC timestamp")
        stamp=datetime.fromisoformat(row["utc"])
        require(stamp.utcoffset()==timedelta(0),"Invalid admission UTC timestamp")
        match(row["growth_bytes"],capacity["budget_bytes"],"six-reference growth")
        volumes=row["volumes"]
        require(isinstance(volumes,dict) and set(volumes)=={"tmp","dependencies","dependency-archives"},
                "Invalid admission volume roster")
        for key,volume in volumes.items():
            require(isinstance(volume,dict) and set(volume)=={"path","free_bytes"}
                    and type(volume["free_bytes"]) is int and volume["free_bytes"]>=0,
                    "Invalid admission volume sample")
            match(volume["path"],capacity["volumes"][key]["path"],"six-reference volume path")
        free=min(volume["free_bytes"] for volume in volumes.values())
        try:
            expected=w.capacity_record(row["total_ram_bytes"],row["available_ram_bytes"],free,
                                       row["growth_bytes"],"six-reference-minimum")
        except w.AvailableRAMShortfall as error:
            require(row["decision"]=="wait" and row["reason"]==str(error),"Invalid available-RAM wait")
        else:
            require(row["decision"]=="admitted" and index==len(lines)-1,"Admission must be the final sample")
            expected["volumes"]=volumes
            match(expected,capacity,"six-reference final admission")
    require(row["decision"]=="admitted","Missing final six-reference admission")


def named_r_values(node,kind):
    require(isinstance(node,list) and len(node)==3 and node[0]==kind and set(node[2])=={"names"},
            "Expected named R values")
    names=node[2]["names"]
    require(names[0]==16 and names[2]=={} and len(names[1])==len(set(names[1]))==len(node[1]),
            "Invalid R object names")
    return dict(zip(names[1],node[1]))


def stable_edition_manifest(path,reference,files,assets):
    manifest=named_r_values(read_r_object(path,2),19)
    created=manifest.pop("created_utc")
    require(created[0]==16 and len(created[1])==1 and isinstance(created[1][0],str)
            and created[1][0] and "path" not in manifest,"Invalid native edition manifest")
    require(manifest["core_version"]==[16,[f'mss-preservation-{reference["format"]}'],{}] and
            canonical_work_path(manifest["reference_path"][1][0])==canonical_work_path(reference["directory"]),
            "Native edition manifest selected another reference")
    match(named_r_values(manifest["checksums"],16),
          {n:r["md5"] for n,r in files.items() if n not in ("checksums.tsv","manifest.rds")},"native manifest edition checksums")
    expected={n:assets[str(PureWindowsPath(reference["directory"])/n)]["md5"] for n in reference["files"] if n!="manifest.rds"}
    expected["fasta"]=assets[reference["fasta"]["path"]]["md5"]
    match(named_r_values(manifest["reference_checksums"],16),expected,"native manifest reference checksums")
    # R identical() compares named list elements in order. Attribute mappings
    # remain order-independent, as in identical() with attrib.as.set=TRUE.
    return list(manifest.items())


def edition_decision(edition,case,reference,files,assets,local):
    from decimal import Decimal
    from verify import decoded_digest,tsv_pairs,NORMALIZED_COLUMNS
    source=assets[case["input"]["path"]]
    decoded=decoded_digest(local(case["input"]["path"]),source)
    match(files["source-artifact"],source,"source artifact")
    match(files["source.tsv"],decoded,"decoded source")
    match(tsv_pairs(edition/"checksums.tsv",("filename","md5")),
          {n:e["md5"] for n,e in files.items() if n not in ("manifest.rds","checksums.tsv")},"edition checksums")
    opt=case["settings"];columns=opt["columns"]
    match(tsv_pairs(edition/"columns.tsv",("field","source_column")),columns,"edition column mapping")
    wanted=dict(core_version=f'mss-preservation-{reference["format"]}',source_rows=str(case["rows"]),normalized_rows=str(case["rows"]),
                input_checksum=source["md5"],checksum_algorithm="MD5",source_file="source.tsv",projection="not_attempted",
                acquisition_receipt="verified" if opt["verify_acquisition_receipt"] else "not_supplied",
                source_strand=("column:"+columns["source_strand"] if "source_strand" in columns else "not_declared"),
                reference_build="SYNTHETIC",reference_version="literal-scientific-v1")
    wanted.update({k:str(opt[k]) for k in ("frequency_type","effect_type","delimiter","batch_bytes","max_record_bytes")})
    metadata=tsv_pairs(edition/"metadata.tsv",("key","value"))
    require(canonical_work_path(metadata["reference_path"])==canonical_work_path(reference["directory"]) and
            all(metadata.get(k)==v for k,v in wanted.items()),"Edition metadata differs from call")
    refchecks=tsv_pairs(edition/"reference_checksums.tsv",("name","md5"))
    expected={n:assets[str(PureWindowsPath(reference["directory"])/n)]["md5"] for n in reference["files"] if n!="manifest.rds"}
    expected["fasta"]=assets[reference["fasta"]["path"]]["md5"]
    match(refchecks,expected,"edition reference checksums")
    result=verify_normalized(edition/"normalized.tsv",case["rows"],local(case["oracle"]["path"]),encoding="utf-8",
                             max_record_bytes=opt["max_record_bytes"],require_full_oracle=True)
    audit=tsv_pairs(edition/"audit.tsv",("status","count"))
    require({k:Decimal(v) for k,v in audit.items()}==result["status_counts"],"Edition audit differs from normalized rows")
    return dict(result,source_copy_verified=True,edition_count=1,normalized_columns=list(NORMALIZED_COLUMNS),
                source_artifact=source,decoded_source=decoded)


def verify_installed_suite(receipt, contract):
    """Verify new-source named coverage; never import a historical pass receipt."""
    from collections import Counter
    files = contract["source_files"]
    require(contract["schema"] == "mss-installed-suite-v1" and contract["legacy_direct_parity_required"] is True,
            "Unknown installed-suite contract")
    root = receipt / "tests"
    plain_path(root, True)
    actual = {path.relative_to(root).as_posix(): digest(path)
              for path in root.rglob("*") if path.is_file()}
    require(set(actual) == set(files), "Installed/copied test membership differs from current checkout")
    for name, pin in files.items():
        plain_path(root / name, False)
        require(actual[name]["bytes"] == pin["bytes"] and actual[name]["sha256"] == pin["sha256"],
                "Installed test bytes differ from current checkout: " + name)
    columns = ["file", "test", "nb", "passed", "skipped", "failed", "error", "warning"]
    cases = read_table(receipt / "tests.tsv", columns)
    require(cases, "Empty installed suite")
    observed = Counter()
    warnings = Counter()
    for row in cases:
        require(row["file"] in files and row["file"].startswith("test-") and row["file"].endswith(".R")
                and row["test"] and row["skipped"] == row["error"] == "FALSE" and row["failed"] == "0",
                "Unexpected/failed/skipped installed case")
        require(all(re.fullmatch("[0-9]+", row[name]) for name in ("nb", "passed", "failed", "warning")),
                "Invalid installed assertion count")
        passed, warning = int(row["passed"]), int(row["warning"])
        require(passed > 0 and int(row["nb"]) == passed + warning, "Empty or inconsistent installed case")
        observed[(row["file"], row["test"], passed)] += 1
        warnings[(row["file"], row["test"])] += warning
    expected = Counter({(file, name, count): copies for file, name, count, copies in contract["case_counts"]})
    require(observed == expected, "New current-checkout case/assertion profile differs from reviewed source intent")
    scripts = {name for name in files if name.startswith("test-") and name.endswith(".R")}
    require({row["file"] for row in cases} == scripts - set(contract["comment_only_scripts"]),
            "An installed test script did not contribute its reviewed cases")
    conditions = read_table(receipt / "conditions.tsv", ["file", "test", "type", "message"])
    detailed_warnings = Counter()
    suite_warnings = 0
    for row in conditions:
        if row["type"] == "suite_warning" and row["file"] == row["test"] == "":
            suite_warnings += 1
        else:
            require(row["type"] == "expectation_warning", "Unexpected installed error/skip/condition")
            detailed_warnings[(row["file"], row["test"])] += 1
    require(detailed_warnings == +warnings and sum(warnings.values()) == contract["expected_test_warnings"]
            and suite_warnings == contract["expected_suite_warnings"], "Installed warning profile differs")
    summary = dict(candidate_suite_pass="TRUE", tests=str(len(cases)), passed=str(sum(int(row["passed"]) for row in cases)),
                   skipped="0", failed="0", errors="0", test_warnings=str(sum(warnings.values())),
                   suite_warnings=str(suite_warnings), source_contract_verification_required="TRUE",
                   legacy_direct_parity_required="TRUE")
    actual_summary, = read_table(receipt / "qualification.tsv", list(summary))
    match(actual_summary, summary, "current installed-suite summary")
    runtime, = read_table(receipt / "runtime.tsv", ["package_path", "package_version", "dll", "dll_md5", "R_version",
                                                  "R_home", "home", "development"])
    require(runtime["development"] == "FALSE", "Developer namespace cannot satisfy installed-package CI")
    return {"installed_suite_verified": True, "tests": len(cases), "passed": int(summary["passed"]),
            "legacy_direct_parity_still_required": True, "runtime": runtime}


def check_probe_executable(path):
    plain_path(path,False)
    with path.open("rb") as stream:
        header=stream.read(64)
        require(len(header)==64 and header[:2]==b"MZ","Missing native probe executable")
        offset=struct.unpack_from("<I",header,60)[0]
        require(64<=offset<=path.stat().st_size-26,"Invalid Windows executable header offset")
        stream.seek(offset)
        h=stream.read(26)
        machine,sections=struct.unpack_from("<HH",h,4)
        optional,flags,magic=struct.unpack_from("<HHH",h,20)
        require(h[:4]==bytes((80,69,0,0)) and machine==0x8664 and 1<=sections<=96 and flags&2 and not flags&0x2000
                and optional>=112 and magic==0x20b and offset+24+optional+40*sections<=path.stat().st_size,
                "Compiled probe lacks a native x64 PE executable header")

def deferral_result(path,mode):
    result=read_r_object(path,3)
    gold={"startup":"437278d9809ed6d55771a2682271aa1888719047bd26569a5185781b22a6f410",
          "bp":"0226c9a4be4103233c2d5df45026c8f0581e7dcad757a2fe8cdde886da96b8af",
          "parallel-vcf":"4360587b6ca27530440d27d5d36b532276ed9a106abdb11c3f0578574a76f3af"}
    warnings=[]
    def collect(value):
        if not isinstance(value,list):return
        if len(value)==3 and value[0]==19 and value[2].get("names")==[16,["value","warnings"],{}]:
            values=value[1]
            require(len(values)==2 and values[1][0]==16 and values[1][2]=={}
                    and all(isinstance(v,str) and v for v in values[1][1]),"Invalid R warning vector")
            warnings.append(values[1][1])
            values[1][1]=[]
        for child_value in value:collect(child_value)
    collect(result)
    identity=hashlib.sha256(json.dumps(result,sort_keys=True,ensure_ascii=False,separators=(",",":")).encode()).hexdigest()
    require(mode in gold and identity==gold[mode] and len(warnings)==(1 if mode=="parallel-vcf" else 3),
            "Mode-specific R values, classes or warning containers differ")
    return identity,warnings

def normalization_phases(path):
    # Both selected Git headers are uninstrumented. The archived driver emits
    # exactly these case, row, result, error and cap lines. Its four fallback
    # counters remain zero because only the separate experiment injected hooks.
    require(path.stat().st_size<=MAX_DIRECTORY,"Oversized normalization transcript")
    records=path.read_bytes().splitlines()
    identity=hashlib.sha256(b"\n".join(records)+b"\n").hexdigest()
    require(len(records)==190 and identity=="1004ffc4fe5e52e807b2f41b51b1e1923194fcb55501b0b2e7b7d5ff7d5dc64b",
            "Native normalization rows, values or error observations differ")
