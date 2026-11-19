"""Inventory and verify the quiescent evidence roots selected for rsync."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import stat
import sys


ROOTS = ("admission.txt", "archives", "guard.jsonl", "inputs", "results")


def paths(root, *, receiving):
    if receiving and {p.name for p in root.iterdir()} != set(ROOTS):
        raise ValueError("Destination top-level entries differ from retained roots")
    entries = {}
    def visit(path):
        mode = path.lstat().st_mode
        name = path.relative_to(root).as_posix()
        if stat.S_ISDIR(mode):
            entries[name] = "directory"
            for child in sorted(path.iterdir()):
                visit(child)
        elif stat.S_ISREG(mode):
            entries[name] = "file"
        else:
            raise ValueError(f"Unexpected retained object type: {name}")
    for name in ROOTS:
        visit(root / name)
    return entries


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as source:
        before = os.fstat(source.fileno())
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
        after = os.fstat(source.fileno())
        if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
            raise ValueError(f"Retained file changed while reading: {path}")
    return digest.hexdigest()


def seal(root, manifest):
    entries = paths(root, receiving=False)
    hashes = {name: sha256(root / name) for name, kind in entries.items() if kind == "file"}
    # Read failures occur before the manifest is opened for publication.
    with manifest.open("x", encoding="utf-8") as f:
        json.dump(dict(schema=1, roots=ROOTS, entries=entries, sha256=hashes), f,
                  indent=2, ensure_ascii=True, sort_keys=True)
        f.write("\n")


def verify(root, manifest):
    with manifest.open(encoding="utf-8") as f:
        expected = json.load(f)
    if expected["schema"] != 1 or expected["roots"] != list(ROOTS):
        raise ValueError("Unknown retained inventory")
    actual = paths(root, receiving=True)
    if actual != expected["entries"]:
        raise ValueError("Destination relative paths or file types differ from source")
    files = {name for name, kind in actual.items() if kind == "file"}
    if files != set(expected["sha256"]):
        raise ValueError("Incomplete retained digest inventory")
    for name in sorted(files):
        if sha256(root / name) != expected["sha256"][name]:
            raise ValueError(f"Destination SHA256 mismatch: {name}")
    print(f"Verified {len(files)} files and {len(actual) - len(files)} directories")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_subparsers(dest="action", required=True)
    actions.add_parser("roots", help="write the canonical NUL-separated rsync root list")
    for name in ("seal", "verify"):
        command = actions.add_parser(name)
        command.add_argument("root", type=Path)
        command.add_argument("manifest", type=Path, help="separate file outside retained roots")
    args = parser.parse_args()
    if args.action == "roots":
        sys.stdout.buffer.write(b"".join(os.fsencode(name) + b"\0" for name in ROOTS))
    elif args.action == "seal":
        seal(args.root, args.manifest)
    else:
        verify(args.root, args.manifest)


if __name__ == "__main__":
    main()
