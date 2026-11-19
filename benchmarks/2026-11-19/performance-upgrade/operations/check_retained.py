"""Small export controls using synthetic files only; retain all control artifacts."""
from pathlib import Path
import shutil
import subprocess
import sys
from unittest.mock import patch

import retained


def reject(action, contains):
    try:
        action()
    except (OSError, ValueError) as e:
        if contains not in str(e):
            raise
    else:
        raise AssertionError(f"Expected rejection: {contains}")


def main(root):
    root.mkdir()
    source = root / "source"
    source.mkdir()
    for name in ("inputs", "archives", "results", "empty-home", "tmp"):
        (source / name).mkdir()
    for name in ("admission.txt", "guard.jsonl", "inputs/panel", "results/space and\nnewline", "tmp/transient"):
        (source / name).write_text("tiny fixture\n")
    (source / "archives/empty-directory").mkdir()
    manifest = root / "inventory.json"
    retained.seal(source, manifest)
    selected = root / "selected"
    selected.mkdir()
    roots = root / "roots.nul"
    with roots.open("xb") as output:
        subprocess.run([sys.executable, "-B", retained.__file__, "roots"], stdout=output, check=True)
    subprocess.run(["rsync", "-a", "--recursive", "--partial", "--from0", f"--files-from={roots}",
                    str(source) + "/", str(selected) + "/"], check=True)
    retained.verify(selected, manifest)
    for where in ("extra", "results/extra"):
        broken = shutil.copytree(selected, root / ("extra-" + where.replace("/", "-")))
        (broken / where).write_text("unexpected")
        reject(lambda: retained.verify(broken, manifest), "differ")
    missing = shutil.copytree(selected, root / "missing", ignore=shutil.ignore_patterns("panel"))
    reject(lambda: retained.verify(missing, manifest), "differ")
    changed = shutil.copytree(selected, root / "changed")
    (changed / "inputs/panel").write_text("different")
    reject(lambda: retained.verify(changed, manifest), "SHA256 mismatch")
    # Inject EIO at the file open boundary because chmod(000) is ineffective
    # under privileged test accounts.
    real_open = Path.open
    def bad_open(path, *args, **kwargs):
        if path.name == "panel":
            raise OSError("injected EIO reading panel")
        return real_open(path, *args, **kwargs)
    rejected = root / "incomplete.json"
    with patch.object(Path, "open", bad_open):
        reject(lambda: retained.seal(source, rejected), "injected EIO")
        reject(lambda: retained.verify(selected, manifest), "injected EIO")
    if rejected.exists():
        raise AssertionError("Source read error published an inventory")
    # A non-regular member is rejected before an attempted read.
    (source / "inputs/dangling").symlink_to("absent")
    reject(lambda: retained.seal(source, rejected), "object type")
    if rejected.exists():
        raise AssertionError("Unsupported object published an inventory")
    print("Retained export controls passed")


if __name__ == "__main__":
    main(Path(sys.argv[1]))
