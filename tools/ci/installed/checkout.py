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


def read_source(path, blob):
    before = plain(path)
    require(before.st_size <= 64 * 1024**2, "Oversized tracked source file")
    with path.open("rb") as stream:
        held = os.fstat(stream.fileno())
        raw = stream.read(64 * 1024**2 + 1)
        after = os.fstat(stream.fileno())
    identity = lambda item: (item.st_dev, item.st_ino, item.st_mode, item.st_nlink,
                             item.st_size, item.st_mtime_ns, item.st_ctime_ns)
    require(identity(before) == identity(held) == identity(after) == identity(path.lstat())
            and len(raw) == before.st_size, "Source changed while frozen")
    require(hashlib.sha1(b"blob " + str(len(raw)).encode("ascii") + b"\0" + raw).hexdigest() == blob,
            "Checkout bytes differ from Git blob; disable autocrlf/filter transformations")
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
        raw = read_source(original, blob)
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
