"""Finite pinned public source acquisition/replay, adapted from accepted producer.

No live version resolution, credentials or package-manager mutation. The caller
supplies the existing owned execution controller and its original clock.
"""
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import shutil
import tarfile
import urllib.error
import urllib.parse
import urllib.request

from evidence import plain_path, require
from verify import digest, write_json

RESERVE = 8_000_000_000
SIX_REFERENCE_DISK_FLOOR = 65_337_374_704


def url(value):
    parsed = urllib.parse.urlsplit(value)
    require(parsed.scheme == "https" and parsed.hostname and not parsed.username and not parsed.password
            and not parsed.query and not parsed.fragment and not any(char.isspace() for char in value),
            "Public credential-free source HTTPS URL required")
    return value


def plan(closure, lane):
    require(lane in {"windows", "linux-gold"}, "Unknown dependency lane")
    platform = "windows" if lane == "windows" else "linux"
    selected = {name for name, row in closure["packages"].items() if row["kind"] == "code"}
    require(len(selected) == 148, "Require the exact148 code package closure")
    selected.update(closure["lanes"][lane]["reference_packages"])
    order = closure["installation_order"][platform]
    require(len(set(order)) == len(order) == 154 and set(order) == set(closure["packages"]),
            "Frozen complete topological package order differs")
    rows = []
    for name in order:
        if name not in selected:
            continue
        row = closure["packages"][name]
        artifact = row["artifacts"][platform]
        require(re.fullmatch(r"[A-Za-z][A-Za-z0-9.]*", name)
                and artifact["archive"] == name + "_" + row["version"] + ".tar.gz"
                and re.fullmatch("[0-9a-f]{64}", artifact["sha256"])
                and str(artifact["bytes"]).isdecimal() and int(artifact["bytes"]) > 0,
                "Malformed pinned public archive")
        source = url(artifact["source"])
        allowed = ("https://cran.r-project.org/src/contrib/", "https://cloud.r-project.org/src/contrib/",
                   *("https://bioconductor.org/packages/3.22/" + section + "/src/contrib/"
                     for section in ("bioc", "data/annotation", "data/experiment", "workflows", "books")))
        require(any(source == prefix + artifact["archive"] for prefix in allowed), "Archive outside official frozen repositories")
        request = source.replace("https://bioconductor.org/packages/3.22/",
                                 "https://bioconductor.posit.co/packages/3.22/", 1)
        rows.append(dict(package=name, version=row["version"], kind=row["kind"], **artifact,
                         request_url=request, unpacked_bytes=closure["unpacked_bytes"][name]))
    require(len(rows) == (154 if lane == "windows" else 150), "Incomplete lane closure")
    return rows


def archive_info(path, row, clock):
    total = 0
    descriptions = []
    seen = set()
    with tarfile.open(path, "r:gz") as archive:
        for member in archive:
            clock.remaining("compute")
            name = PurePosixPath(member.name)
            require(name.parts and name.parts[0] == row["package"] and not name.is_absolute()
                    and ".." not in name.parts and "\\" not in member.name
                    and member.name not in seen and len(seen) < 100000
                    and (member.isdir() or member.isfile() or member.issym() or member.islnk()),
                    "Unexpected pinned source archive object")
            seen.add(member.name)
            if member.issym() or member.islnk():
                target = PurePosixPath(member.linkname)
                require(not target.is_absolute() and ".." not in target.parts,
                        "Indirect source archive object escapes package")
            total += member.size
            require(total <= row["unpacked_bytes"], "Source expansion exceeds fixed actual archive census")
            if member.name == row["package"] + "/DESCRIPTION":
                require(member.isfile() and member.size <= 1048576, "Invalid source DESCRIPTION")
                descriptions.append(archive.extractfile(member).read().decode("utf-8"))
    require(total == row["unpacked_bytes"] and len(descriptions) == 1
            and all(re.search("(?m)^" + field + r":\s*" + re.escape(row[key]) + r"\s*$", descriptions[0])
                    for field, key in (("Package", "package"), ("Version", "version"))),
            "Pinned archive DESCRIPTION or complete expansion changed")
    return total


def open_source(row, clock, attempts):
    """Only404/410 may try the identical version under official CRAN Archive."""
    candidates = [row["request_url"]]
    source = urllib.parse.urlsplit(row["source"])
    if (row["kind"] == "code" and source.hostname in {"cran.r-project.org", "cloud.r-project.org"}
            and source.path == "/src/contrib/" + row["archive"]):
        candidates.append("https://" + source.hostname + "/src/contrib/Archive/" + row["package"] + "/" + row["archive"])
    for index, candidate in enumerate(candidates):
        clock.remaining("compute")
        request = urllib.request.Request(candidate, headers={"User-Agent": "MungeSumstats-installed-CI/1"})
        try:
            response = urllib.request.urlopen(request, timeout=min(60, clock.remaining("compute")))
        except urllib.error.HTTPError as error:
            attempts.write(json.dumps(dict(package=row["package"], source=row["source"],
                attempted_url=candidate, http_status=error.code, outcome="HTTP_error")) + "\n")
            attempts.flush()
            if error.code not in {404, 410} or index + 1 == len(candidates):
                raise
            error.close()
        except BaseException as error:
            attempts.write(json.dumps(dict(package=row["package"], source=row["source"],
                attempted_url=candidate, http_status=None, outcome=type(error).__name__)) + "\n")
            attempts.flush()
            raise
        else:
            final = url(response.url)
            attempts.write(json.dumps(dict(package=row["package"], source=row["source"],
                attempted_url=candidate, final_url=final, http_status=response.status,
                outcome="opened")) + "\n")
            attempts.flush()
            return response, candidate
    raise RuntimeError("All exact frozen source locations were unavailable")


def acquire(execution, closure, lane, cache, planning):
    rows = plan(closure, lane)
    require(not cache.exists() and not planning.exists(), "Fresh dependency source/metadata roots")
    cache.mkdir()
    planning.mkdir()
    archive_bytes = sum(int(row["bytes"]) for row in rows)
    expanded = [row["unpacked_bytes"] for row in rows]
    # The accepted source-install bound funds all retained sources, accumulated
    # library, one serial unpacked build and generated objects. Never assume the
    # advertised runner disk specification is an admission observation.
    growth = archive_bytes + 2 * sum(expanded) + max(expanded) + 2_000_000_000
    if lane == "windows":
        growth = max(growth, SIX_REFERENCE_DISK_FLOOR - RESERVE)
    capacity = execution.capacity("six-reference-minimum" if lane == "windows" else "two-reference-minimum", growth)
    write_json(planning / "acquisition-capacity.json", dict(lane=lane, archive_bytes=archive_bytes,
        total_expansion_bytes=sum(expanded), additional_growth_bytes=growth, admission=capacity))
    with (planning / "verified-get.jsonl").open("x", encoding="utf-8") as log, \
            (planning / "source-attempts.jsonl").open("x", encoding="utf-8") as attempts:
        for row in rows:
            path = cache / row["archive"]
            observed = None
            try:
                execution.clock.remaining("compute")
                response, actual_request_url = open_source(row, execution.clock, attempts)
                with response:
                    final_url = url(response.url)
                    require(response.status == 200, "Exact pinned source requires HTTP200")
                    size = int(row["bytes"])
                    count = 0
                    with path.open("xb") as stream:
                        while True:
                            execution.clock.remaining("compute")
                            require(shutil.disk_usage(cache).free >= RESERVE, "Dependency acquisition disk reserve reached")
                            block = response.read(min(1048576, size + 1 - count))
                            if not block:
                                break
                            stream.write(block)
                            count += len(block)
                            require(count <= size, "Source exceeds frozen byte bound")
                    observed = digest(path)
                    require(count == size and observed["sha256"] == row["sha256"], "Frozen source bytes/SHA256 differ")
                    archive_info(path, row, execution.clock)
                    log.write(json.dumps(dict(row, verified=True, actual_request_url=actual_request_url, final_url=final_url, observed=observed)) + "\n")
                    log.flush()
            except BaseException as error:
                write_json(planning / "acquisition-failure.json", dict(package=row["package"], source=row["source"],
                    request_url=row["request_url"], error_type=type(error).__name__, error=str(error), observed=observed))
                raise
    return rows


def install(execution, rows, cache, library, r):
    require(not library.exists(), "Fresh dependency installation library required")
    library.mkdir()
    execution.environment["R_LIBS_USER"] = str(library)
    for row in rows:
        source = cache / row["archive"]
        require(digest(source)["sha256"] == row["sha256"], "Dependency source changed before ordinary install")
        execution.run("dependency-" + row["package"], [r, "CMD", "INSTALL", "--preclean",
                      "--library=" + str(library), source], timeout=1800,
                      provenance=dict(package=row["package"], version=row["version"], archive=digest(source)))


def reference_payloads(execution, rows, cache, library):
    result = {}
    for row in rows:
        if row["kind"] != "reference":
            continue
        source = cache / row["archive"]
        expected = {}
        with tarfile.open(source, "r:gz") as archive:
            for member in archive:
                execution.clock.remaining("compute")
                prefix = row["package"] + "/inst/extdata/"
                if member.name.startswith(prefix) and member.isfile():
                    relative = member.name[len(row["package"] + "/inst/"):]
                    require(relative not in expected, "Repeated source reference payload")
                    h = hashlib.sha256()
                    total = 0
                    with archive.extractfile(member) as stream:
                        while block := stream.read(1048576):
                            execution.clock.remaining("compute")
                            h.update(block)
                            total += len(block)
                    require(total == member.size, "Truncated source reference payload")
                    expected[relative] = {"bytes": total, "sha256": h.hexdigest()}
        installed = {}
        payload = library / row["package"] / "extdata"
        plain_path(payload, True)
        for path in payload.rglob("*"):
            execution.clock.remaining("compute")
            if path.is_dir():
                plain_path(path, True)
                continue
            plain_path(path, False)
            identity = digest(path)
            installed[path.relative_to(library / row["package"]).as_posix()] = {
                "bytes": identity["bytes"], "sha256": identity["sha256"]}
        require(expected and expected == installed, "Complete reference extdata source/installed payload mismatch")
        result[row["package"]] = {"version": row["version"],
            "archive": {"bytes": int(row["bytes"]), "sha256": row["sha256"]},
            "source_payload": expected, "installed_payload": installed,
            "installed_payload_matches_source": True}
    return result
