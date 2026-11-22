"""Normal current-checkout R build/install and exact included-code receipt."""
import hashlib
import json
from pathlib import Path, PurePosixPath
import shutil
import tarfile

from evidence import plain_path, require
from verify import digest, write_json


def description(path):
    data = {}
    key = None
    for line in path.read_text(encoding="utf-8").splitlines():
        if line[:1].isspace():
            require(key is not None, "Unbound DESCRIPTION continuation")
            data[key] += " " + line.strip()
        elif line:
            key, separator, value = line.partition(":")
            require(separator and key not in data, "Malformed source DESCRIPTION")
            data[key] = value.strip()
    require(data["Package"] == "MungeSumstats" and data["Version"], "Expected actual MSS checkout DESCRIPTION")
    return data


def pin(path):
    plain_path(path, False)
    return dict(path=str(path), **digest(path))


def included_source(archive, checkout, package):
    observed = {}
    total = 0
    with tarfile.open(archive, "r:gz") as stream:
        for member in stream:
            name = PurePosixPath(member.name)
            require(name.parts and name.parts[0] == package and not name.is_absolute()
                    and ".." not in name.parts and "\\" not in member.name
                    and (member.isdir() or member.isfile()), "Unexpected current R build archive member")
            if member.isdir():
                continue
            relative = name.relative_to(package).as_posix()
            require(relative not in observed and member.size <= 64 * 1024**2
                    and len(observed) < 50000, "Bounded unique source build members required")
            with stream.extractfile(member) as incoming:
                raw = incoming.read(64 * 1024**2 + 1)
            total += len(raw)
            require(len(raw) == member.size and total <= 512 * 1024**2, "Current build archive byte bound")
            observed[relative] = {"bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()}
    wanted = {name: {"bytes": item["bytes"], "sha256": item["sha256"]}
              for name, item in checkout["files"].items()
              if name.startswith(("R/", "src/", "tests/testthat/")) or name in {"NAMESPACE", "tests/testthat.R"}}
    actual_code = {name for name in observed if name.startswith(("R/", "src/", "tests/testthat/"))
                   or name in {"NAMESPACE", "tests/testthat.R"}}
    require(wanted and actual_code == set(wanted) and all(observed.get(name) == value for name, value in wanted.items()),
            "R CMD build omitted/changed current native/R/installed-test source")
    return dict(included_code_and_test_files=wanted, archive_files=observed, complete_included_code_matches_checkout=True)


def build_install(execution, checkout, dependencies, r_command, rscript, ci_source, evidence):
    frozen = Path(checkout["source"])
    build = evidence / "build"
    build.mkdir()
    work = execution.root / "package-build"
    work.mkdir()
    # R's normal build/install may write derived files. The frozen source copy
    # remains immutable; compiler input comes from this exact verified copy.
    source = work / "source"
    shutil.copytree(frozen, source)
    for name, item in checkout["files"].items():
        actual = digest(source / name)
        require(actual["bytes"] == item["bytes"] and actual["sha256"] == item["sha256"],
                "Current build work copy changed")
    doc = description(source / "DESCRIPTION")
    build_log = execution.run("current-R-CMD-build", [r_command, "CMD", "build", "--no-build-vignettes",
        "--no-manual", source], cwd=work, provenance=dict(checkout_commit=checkout["commit"], checkout_tree=checkout["tree"]))
    archive = work / (doc["Package"] + "_" + doc["Version"] + ".tar.gz")
    require({path.name for path in work.glob("*.tar.gz")} == {archive.name},
            "Normal R build did not produce the one expected current archive")
    binding = included_source(archive, checkout, doc["Package"])
    retained = build / archive.name
    shutil.copyfile(archive, retained)
    require(digest(retained) == digest(archive), "Generated current R archive retention changed")
    library = execution.root / "package-library"
    library.mkdir()
    install_argv = [r_command, "CMD", "INSTALL", "--install-tests", "--preclean",
                    "--library=" + str(library), retained]
    install_log = execution.run("current-R-CMD-INSTALL", install_argv,
        provenance=dict(checkout_commit=checkout["commit"], generated_current_archive=pin(retained)))
    receipt = build / "installation.json"
    execution.run("current-installed-identity", [rscript, "--vanilla", ci_source / "install-receipt.R",
        library, dependencies, evidence / "checkout/checkout.json", ci_source / "closure.json",
        json.loads((execution.root / "request.json").read_bytes())["lane"], receipt],
        provenance=dict(checkout_commit=checkout["commit"], driver=pin(ci_source / "install-receipt.R")))
    actual = json.loads(receipt.read_bytes())
    require(actual["schema"] == "mss-current-install-v1" and actual["normal_current_checkout"] is True
            and actual["checkout_commit"] == checkout["commit"] and actual["checkout_tree"] == checkout["tree"]
            and actual["source_files"] == checkout["files"] and actual["development"] is False,
            "Actual current installation/source identity differs")
    write_json(build / "build.json", dict(schema="mss-normal-build-v1", checkout_commit=checkout["commit"],
        checkout_tree=checkout["tree"], generated_current_archive=pin(retained), source_binding=binding,
        build_stdout=pin(build_log), build_stderr=pin(build_log.with_suffix(".stderr")),
        build_execution=pin(build_log.with_suffix(".json")), install_stdout=pin(install_log),
        install_stderr=pin(install_log.with_suffix(".stderr")), install_execution=pin(install_log.with_suffix(".json")),
        install_argv=[str(value) for value in install_argv], installed_identity=pin(receipt),
        historical_archive_used=False, tests_installed_by_R_CMD_INSTALL=True))
    return actual
