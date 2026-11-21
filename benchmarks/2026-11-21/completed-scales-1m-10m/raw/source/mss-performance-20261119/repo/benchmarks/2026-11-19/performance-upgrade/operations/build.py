"""Build both normal installations offline from verified uncompressed git archives."""
import argparse
import csv
import hashlib
import os
from pathlib import Path
import platform
import re
import shlex
import shutil
import subprocess
import sys
import tarfile

BENCH = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BENCH))
from run import BASELINE, package_tree
from verify import digest, write_json


def pin(path):
    path = Path(path).resolve(strict=True)
    return {"path": str(path), "sha256": digest(path)["sha256"]}


def checked(argv):
    return subprocess.check_output(argv, text=True).strip()


def inventory_library(path):
    # Distribution libraries may link files into /usr/share. Bind both link
    # spelling and dereferenced content; the package under comparison instead
    # uses the runner's strict package_tree definition throughout.
    identities = {}
    for package in sorted(path.iterdir()):
        if not (package / "DESCRIPTION").is_file():
            continue
        identity = hashlib.sha256()
        def visit(directory, ancestors):
            resolved = directory.resolve(strict=True)
            if resolved in ancestors:
                raise ValueError(f"Cyclic dependency directory: {directory}")
            for entry in sorted(directory.iterdir()):
                identity.update(entry.relative_to(path).as_posix().encode() + b"\0")
                if entry.is_symlink():
                    identity.update(b"link\0" + os.fsencode(os.readlink(entry)) + b"\0")
                if entry.is_dir():
                    identity.update(b"dir\0")
                    visit(entry, ancestors | {resolved})
                elif entry.is_file():
                    identity.update(digest(entry)["sha256"].encode() + b"\0")
                else:
                    raise ValueError(f"Unsupported dependency entry: {entry}")
        visit(package, set())
        identities[package.name] = identity.hexdigest()
    return identities


def archive_source(archive, expected_hash, commit, target):
    actual = pin(archive)
    if actual["sha256"] != expected_hash:
        raise ValueError("Source archive SHA256 differs from transfer receipt")
    with archive.open("rb") as inp:
        embedded = subprocess.check_output(["git", "get-tar-commit-id"], stdin=inp, text=True).strip()
    if embedded != commit:
        raise ValueError("Archive commit differs from the approved full Git ID")
    target.mkdir()
    with tarfile.open(archive, "r:") as source:
        # The operator supplies plain git archive output without --prefix.
        for member in source:
            if (Path(member.name).is_absolute() or ".." in Path(member.name).parts
                    or not (member.isdir() or member.isfile())):
                raise ValueError(f"Unexpected archived path or object: {member.name}")
            source.extract(member, target, filter="data")
    if not (target / "src/Makevars").is_file():
        raise ValueError("Expected MungeSumstats package sources at archive root")
    return actual


def runtime_roles(root, rscript, engine):
    output = root / "runtime.tsv"
    subprocess.run([rscript, "--vanilla", BENCH / "worker.R", "runtime", output], check=True)
    with output.open(encoding="utf-8", newline="") as stream:
        runtime = {row["key"]: row["value"] for row in csv.DictReader(stream, delimiter="\t")}
    if Path(runtime["R_engine"]).resolve(strict=True) != engine:
        raise ValueError("Pinned launcher and observed R engine differ")
    roles = {"r_shared": "R_shared", "blas": "BLAS_library", "lapack": "LAPACK_library"}
    return [dict(role=role, **pin(runtime[field])) for role, field in roles.items()
            if runtime[field] != "unavailable"]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("root", type=Path, help="new output directory")
    ap.add_argument("--candidate", required=True)
    ap.add_argument("--archives", type=Path, required=True)
    for role in ("baseline", "candidate"):
        ap.add_argument(f"--{role}-sha256", required=True)
    ap.add_argument("--dependencies", required=True, type=Path,
                    help="shared R package library in the pinned Linux image")
    ap.add_argument("--image-digest", required=True, help="immutable provisioned image sha256:...")
    args = ap.parse_args()
    if platform.system() != "Linux" or not re.fullmatch("sha256:[0-9a-f]{64}", args.image_digest):
        raise ValueError("Require Linux and an immutable pre-provisioned image digest")
    if not re.fullmatch("[0-9a-f]{40}", args.candidate) or args.candidate == BASELINE:
        raise ValueError("Supply a distinct immutable candidate commit")
    root = args.root.resolve()
    root.mkdir()
    libs = args.dependencies.resolve(strict=True)
    if (libs / "MungeSumstats").exists() or not (libs / "testthat" / "DESCRIPTION").is_file():
        raise ValueError("Need shared checked dependency library including testthat, without MSS")
    # Run from a clean shell (README). No personal R startup / Makevars.
    os.environ.update(R_LIBS_USER=str(libs), R_LIBS_SITE=str(libs), R_MAKEVARS_USER="/dev/null",
                      R_ENVIRON_USER="/dev/null", R_PROFILE_USER="/dev/null", MAKEFLAGS="-j1",
                      LANG="C.UTF-8", LC_ALL="C.UTF-8", TZ="UTC")
    r = Path(shutil.which("R")).resolve(strict=True)
    r_home = Path(checked([r, "RHOME"])).resolve(strict=True)
    rscript = (r_home / "bin/Rscript").resolve(strict=True)
    engine = (r_home / "bin/exec/R").resolve(strict=True)
    def expr(code):
        return checked([rscript, "--vanilla", "-e", code])
    # Preserve the complete effective configuration without requiring optional
    # make variables from a different R release to exist.
    config = checked([r, "CMD", "config", "--all"])
    (root / "R-CMD-config.txt").write_text(config + "\n")
    flags = {k: checked([r, "CMD", "config", k])
             for k in ("CC", "CFLAGS", "CPPFLAGS", "CXX17", "CXX17FLAGS", "LDFLAGS")}
    compiler = Path(shutil.which(shlex.split(flags["CXX17"])[0])).resolve(strict=True)
    common = {"compiler": checked([compiler, "--version"]), "compiler_file": pin(compiler),
              "compiler_flags": flags, "R_version": expr("cat(R.version.string)"),
              "R_platform": expr("cat(R.version$platform)"), "image_digest": args.image_digest,
              "makeconf": pin(r_home / "etc/Makeconf"), "all_flags": pin(root / "R-CMD-config.txt"),
              "dependency_trees": inventory_library(libs),
              "base_library_trees": inventory_library(r_home / "library")}
    packages = expr("write.table(utils::installed.packages()[,c('Package','Version','LibPath')],"
                    "file=stdout(),sep='\\t',row.names=FALSE)")
    # The tabular inventory also records system package versions.
    (root / "dpkg.tsv").write_text(checked(["dpkg-query", "-W", "-f=${binary:Package}\\t${Version}\\n"])+"\n")
    (root / "r-packages.txt").write_text(packages+"\n")
    roles = {}
    for name, commit in [("baseline", BASELINE), ("candidate", args.candidate)]:
        directory = root / name
        directory.mkdir()
        source = directory / "source"
        archive = archive_source(args.archives / f"{commit}.tar",
                                 getattr(args, name + "_sha256"), commit, source)
        library = directory / "library"
        library.mkdir()
        command = [str(r), "CMD", "INSTALL", "--preclean", f"--library={library}", str(source)]
        log = directory / "install.log"
        with log.open("xb") as transcript:
            subprocess.run(command, stdout=transcript, stderr=subprocess.STDOUT, check=True)
        package = library / "MungeSumstats"
        dll = pin(package / "libs/MungeSumstats.so")
        tree = package_tree(package)
        write_json(directory / "installed-files.json", tree)
        build = directory / "build.json"
        write_json(build, {"source_commit": commit, "source_archive": archive, "install_argv": command,
                           "install_log": pin(log), "source_makevars": pin(source / "src/Makevars"),
                           "toolchain": common, "package_tree_sha256": tree["tree_sha256"],
                           "dll_sha256": dll["sha256"]})
        roles[name] = dict(commit=commit, library=str(library), dll=dll, build_receipt=pin(build),
                           package_tree_sha256=tree["tree_sha256"])
    if inventory_library(libs) != common["dependency_trees"]:
        raise ValueError("Dependency tree changed during builds")
    linked = runtime_roles(root, rscript, engine)
    write_json(root / "installations.json", {**roles, "dependency_libraries": [str(libs)],
               "r_executable": pin(rscript), "r_engine": pin(engine), "runtime_files": linked})


if __name__ == "__main__":
    main()
