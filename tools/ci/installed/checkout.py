"""Freeze/copy the actual triggering checkout before ordinary R build/install.

No historic package archive, network, credentials or package execution.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import stat
import subprocess


def require(value, message):
    if not value:
        raise ValueError(message)


def git(repo, *args):
    return subprocess.check_output(["git", "-C", str(repo), *args], timeout=60)


def plain(path, directory=False):
    info = path.lstat()
    require(path.resolve(strict=True) == path and not path.is_symlink()
            and not getattr(info, "st_file_attributes", 0) & 0x400
            and (stat.S_ISDIR(info.st_mode) if directory else stat.S_ISREG(info.st_mode) and info.st_nlink == 1),
            "Indirect/nonordinary source input: " + str(path))
    return info


def read_source(path, blob, *, relative_path=None, diagnostics=None):
    before = plain(path)
    require(before.st_size <= 64 * 1024**2, "Oversized tracked source file")
    with path.open("rb") as stream:
        held = os.fstat(stream.fileno())
        raw = stream.read(64 * 1024**2 + 1)
        after = os.fstat(stream.fileno())
    final = path.lstat()
    samples = dict(path_before=before, fd_held=held, fd_after=after, path_after=final)
    common = ("st_dev", "st_ino", "st_mode", "st_nlink", "st_size", "st_mtime_ns")
    # CPython3.12 Windows lstat ctime is CreationTime, but fstat ctime is
    # ChangeTime. Compare creation across APIs and each ctime within its domain.
    timestamp = "st_birthtime_ns" if os.name == "nt" else "st_ctime_ns"
    available = all(type(getattr(item, timestamp, None)) is int for item in samples.values())
    identity = lambda item: tuple(getattr(item, field, None) for field in common + (timestamp,))
    stable = available and identity(before) == identity(held) == identity(after) == identity(final)
    stable = stable and before.st_ctime_ns == final.st_ctime_ns and held.st_ctime_ns == after.st_ctime_ns
    complete = len(raw) == before.st_size
    actual_blob = hashlib.sha1(b"blob " + str(len(raw)).encode("ascii") + b"\0" + raw).hexdigest()
    blob_matches = actual_blob == blob
    cross_domain = os.name == "nt" and (before.st_ctime_ns != held.st_ctime_ns or
                                         final.st_ctime_ns != after.st_ctime_ns)
    diagnostic = None
    if not (stable and complete and blob_matches) or cross_domain:
        fields = common + ("st_ctime_ns", "st_birthtime_ns", "st_file_attributes", "st_reparse_tag")
        snapshots = {name: {field: getattr(item, field, None) for field in fields}
                     for name, item in samples.items()}
        import sys
        diagnostic = dict(schema="mss-source-identity-observation-v1",
            status="windows_cross_domain_only" if stable and complete and blob_matches else "rejected_source",
            tracked_path=(relative_path if relative_path is not None else path.name)[:1024],
            cross_api_timestamp=timestamp, common_identity_stable=stable, read_complete=complete,
            path_ctime_stable=before.st_ctime_ns == final.st_ctime_ns,
            fd_ctime_stable=held.st_ctime_ns == after.st_ctime_ns,
            snapshots=snapshots, read_bytes=len(raw), expected_git_blob=blob,
            actual_git_blob=actual_blob, raw_sha256=hashlib.sha256(raw).hexdigest(),
            python_version=sys.version[:512], python_executable=sys.executable[:1024])
        encoded = json.dumps(diagnostic, sort_keys=True)
        require(len(encoded.encode("utf-8")) <= 16384, "Oversized source identity diagnostic")
        if diagnostics is not None:
            name = "windows-stat-domains.json" if diagnostic["status"] == "windows_cross_domain_only" else "source-identity-failure.json"
            try:
                with (diagnostics / name).open("x", encoding="utf-8") as stream:
                    stream.write(encoded + "\n")
            except FileExistsError:
                pass  # Retain the first bounded observation in this fresh freeze.
    detail = "" if diagnostic is None else ": " + json.dumps(diagnostic, sort_keys=True)
    require(stable and complete, "Source changed while frozen" + detail)
    require(blob_matches, "Checkout bytes differ from Git blob; disable autocrlf/filter transformations" + detail)
    return raw


def freeze_checkout(checkout, output):
    checkout = Path(checkout).resolve(strict=True)
    output = Path(output).absolute()
    plain(checkout, True)
    require(not output.exists() and not output.is_symlink(), "Fresh source receipt directory required")
    require(git(checkout, "diff", "--name-only", "HEAD") == b"", "Modified tracked checkout")
    require(git(checkout, "diff", "--cached", "--name-only") == b"", "Modified Git index")
    commit = git(checkout, "rev-parse", "HEAD").decode("ascii").strip()
    tree = git(checkout, "rev-parse", "HEAD^{tree}").decode("ascii").strip()
    require(re.fullmatch("[0-9a-f]{40}", commit) and re.fullmatch("[0-9a-f]{40}", tree), "GitHub SHA1 source identity required")
    listing = git(checkout, "ls-tree", "-r", "-z", "--full-tree", "HEAD")
    rows = []
    for record in listing.split(b"\0"):
        if not record:
            continue
        metadata, raw_name = record.split(b"\t", 1)
        mode, kind, blob = metadata.decode("ascii").split()
        name = raw_name.decode("utf-8")
        path = PurePosixPath(name)
        require(kind == "blob" and mode in {"100644", "100755"} and not path.is_absolute()
                and ".." not in path.parts and "\\" not in name and len(raw_name) <= 1024,
                "Unsupported tracked source entry")
        rows.append((name, mode, blob))
    require(rows and len(rows) <= 50000 and len({name for name, _, _ in rows}) == len(rows), "Tracked source roster bound")
    output.mkdir(mode=0o700)
    source = output / "source"
    source.mkdir(mode=0o700)
    files = {}
    total = 0
    for name, mode, blob in rows:
        original = checkout / name
        for parent in original.parents:
            if parent == checkout:
                break
            plain(parent, True)
        raw = read_source(original, blob, relative_path=name, diagnostics=output)
        total += len(raw)
        require(total <= 512 * 1024**2, "Tracked source byte bound")
        target = source / name
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open("xb") as stream:
            stream.write(raw)
        target.chmod(0o700 if mode == "100755" else 0o600)
        files[name] = {"git_blob": blob, "git_mode": mode, "bytes": len(raw),
                       "sha256": hashlib.sha256(raw).hexdigest()}
    require(git(checkout, "rev-parse", "HEAD").decode("ascii").strip() == commit
            and git(checkout, "diff", "--name-only", "HEAD") == b"", "Checkout changed during freeze")
    receipt = {"schema": "mss-checkout-v1", "commit": commit, "tree": tree, "files": files,
               "bytes": total, "source": str(source), "normal_current_checkout": True}
    with (output / "checkout.json").open("x", encoding="utf-8") as stream:
        json.dump(receipt, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return receipt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkout", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    receipt = freeze_checkout(arguments.checkout, arguments.output)
    print(json.dumps({key: receipt[key] for key in ("schema", "commit", "tree", "bytes", "normal_current_checkout")}, sort_keys=True))
