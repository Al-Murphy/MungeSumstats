"""Bounded Linux source-sealed upload; no live scientific tree enumeration."""
import hashlib
from pathlib import Path
import shutil
import tempfile

import windows_job as w
from evidence import document, require
from verify import digest, write_json

LINUX_LIMITS = {name:1 << 20 for name in ("preflight.json", "job-start.json", "q/deadlines.json",
    "q/linux-stop.json", "q/artifact-status.json", "q/UNSEALED.json", "q/upload-omissions.json")}
LINUX_LIMITS["q/artifact-inventory.json"] = 16 << 20
SUCCESS = set(LINUX_LIMITS) - {"q/UNSEALED.json", "q/upload-omissions.json"}


def prepare_linux_snapshot(attempt, area):
    attempt, area = attempt.absolute(), area.absolute()
    require(area != attempt and attempt not in area.parents, "Snapshot outside scientific attempt required")
    with w.physical_directory(area):
        private = Path(tempfile.mkdtemp(prefix="mss-linux-upload-", dir=area))
    payload = private / "payload"
    payload.mkdir()
    pending = private / "pending"
    pending.mkdir()
    omitted = {}
    sequence = 0
    def take(name, limit, live=False):
        nonlocal sequence
        sequence += 1
        target, incoming, source = payload / name, pending / str(sequence), attempt / name
        try:
            with w.physical_directory(source.parent) as fd:
                _, truncated = w.transfer_regular(source, fd, limit, incoming, live)
            target.parent.mkdir(parents=True, exist_ok=True)
            incoming.rename(target)
            if truncated:
                omitted[name] = "bounded observed prefix retained"
        except FileNotFoundError:
            omitted[name] = "absent"
        except (OSError, ValueError, RuntimeError) as error:
            omitted[name] = "unreadable or nonphysical: " + type(error).__name__
    for name, limit in LINUX_LIMITS.items():
        if name != "q/upload-omissions.json":
            take(name, limit)
    try:
        status = document(payload / "q/artifact-status.json")
        stop = document(payload / "q/linux-stop.json")
        require(omitted.get("q/UNSEALED.json") == "absent" and status["source_sealed"] is True
                and status["receiver_verified"] is False and status["lane"] == "linux-gold"
                and stop["owned_tree_quiescent"] is True, "Missing Linux source seal/closure")
        require(not any(name in omitted for name in SUCCESS), "Incomplete source-sealed Linux outside receipts")
        inventory = payload / "q/artifact-inventory.json"
        require(digest(inventory)["sha256"] == status["inventory"]["sha256"], "Source inventory identity differs")
        expected = document(inventory)
        require(0 <= expected["bytes"] <= w.ARTIFACT_LIMIT and shutil.disk_usage(area).free >= w.RESERVE + expected["bytes"],
                "Bounded Linux snapshot capacity")
        complete = private / "complete-evidence"
        complete.mkdir()
        actual = w.evidence_inventory(attempt / "q/evidence", destination=complete)
        require(actual == expected, "Linux immutable upload copy differs from independent source seal")
        complete.rename(payload / "q/evidence")
        return payload, True
    except (KeyError, TypeError, ValueError, OSError, RuntimeError) as error:
        omitted["q/evidence"] = "full tree withheld: " + type(error).__name__
    # No recursive read follows absent/unconfirmed closure. Existing source-only
    # failure and candidate metadata may contribute bounded observed prefixes.
    for name in ("failure.json", "qualification.json", "qualification-candidate.json"):
        take("q/evidence/" + name, 1 << 20, True)
    (payload / "q").mkdir(exist_ok=True)
    write_json(payload / "q/upload-omissions.json", dict(full_tree_uploaded=False, omitted=omitted))
    return payload, False
