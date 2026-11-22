"""Validate frozen output and query order; report signed candidate-minus-baseline differences."""
import json
import math
import statistics
import sys
from pathlib import Path

SOURCE = Path("/workspace/mss-key-page-8k-experiment-20261122")
ORDERS = ("matched", "shuffled")
ROLES = ("baseline", "candidate")
SCHEDULE = json.loads((SOURCE / "expected/order.json").read_bytes())
GOLD = SOURCE / "expected"


def read_timing(prefix, mode, order):
    receipt = json.loads(prefix.with_suffix(".timing.json").read_bytes())
    elapsed = receipt.pop("harmonise_ns")
    batch = receipt.pop("batch_rank", None)
    if prefix.name.endswith("-candidate"):
        required = {"attempts", "admitted", "arena_peak_bytes", "gather_estimate_bytes", "gather_max_keys", "gather_max_key_bytes"}
        if not isinstance(batch, dict) or set(batch) != required or any(type(v) is not int or v < 0 for v in batch.values()):
            raise ValueError("Missing candidate batch-rank diagnostics")
        if batch["attempts"] != receipt["batches"] or batch["admitted"] > batch["attempts"] or batch["arena_peak_bytes"] > 32*1024**2 or batch["gather_max_keys"] > 200000 or batch["gather_max_key_bytes"] > 224:
            raise ValueError("Invalid batch-rank bounds")
        if mode == "panel" and batch["admitted"] != 10:
            raise ValueError("Panel prefetch did not admit all ten batches")
    elif batch is not None:
        raise ValueError("Baseline must be unmodified ordinary row processing")
    if receipt != SCHEDULE[mode][order] or type(elapsed) is not int:
        raise ValueError(f"Incomplete batch or query order at {prefix}")
    if elapsed < 0 or (mode == "empty") != (elapsed == 0):
        raise ValueError(f"Invalid native duration at {prefix}")
    return elapsed


def controls(root):
    checked = []
    for role in ROLES:
        if (root / f"lookup-parity-{role}.txt").read_bytes() != (GOLD / "lookup-parity.txt").read_bytes():
            raise ValueError(f"Lookup parity failed for {role}")
        for order in ORDERS:
            for fmt in (2, 3):
                for mode in ("fixture", "short", "empty"):
                    prefix = root / f"{mode}-format{fmt}-{order}-{role}"
                    read_timing(prefix, mode, order)
                    if mode == "fixture":
                        value = prefix.with_suffix(".json").read_bytes()
                        expected = (GOLD / f"fixture-format{fmt}.json").read_bytes()
                    else:
                        value = prefix.with_suffix(".rows.tsv").read_bytes()
                        expected = (GOLD / "short.rows.tsv").read_bytes() if mode == "short" else b""
                    if value != expected:
                        raise ValueError(f"Output mismatch at {prefix}")
                    checked.append(prefix.name)
    return {"passed": checked, "lookup_parity": list(ROLES)}


def panel(root, order, pair, role):
    prefix = root / f"panel-{order}-{pair}-{role}"
    elapsed = read_timing(prefix, "panel", order) / 1e9
    if prefix.with_suffix(".json").read_bytes() != (GOLD / "panel.json").read_bytes():
        raise ValueError(f"Panel output mismatch at {prefix}")
    receipt = json.loads(prefix.with_suffix(".process.json").read_bytes())
    expected_command = [str(root / role), "panel", "/workspace/mss-task4/inputs/reference/matched-index",
                        "/workspace/mss-task4/inputs/reference/GRCh37.ensembl104.chr22.fa", order]
    if receipt["exit_code"] != 0 or receipt["command"] != expected_command or receipt["pid"] <= 0:
        raise ValueError(f"Incomplete wait4 receipt at {prefix}")
    metrics = {k: receipt[k] for k in ("wall_seconds", "user_seconds", "system_seconds", "peak_rss_bytes")}
    metrics.update(native_seconds=elapsed, cpu_seconds=metrics["user_seconds"] + metrics["system_seconds"])
    if any(not math.isfinite(v) or v < 0 for v in metrics.values()):
        raise ValueError(f"Invalid native process metrics at {prefix}")
    if metrics["peak_rss_bytes"] <= 0 or metrics["wall_seconds"] < elapsed or metrics["cpu_seconds"] <= 0:
        raise ValueError(f"Inconsistent native process metrics at {prefix}")
    return metrics


def medians(pairs, field):
    return {metric: statistics.median(p[field][metric] for p in pairs)
            for metric in pairs[0][field] if all(p[field][metric] is not None for p in pairs)}


def paired(root):
    if json.loads((root / "controls.json").read_bytes()) != controls(root):
        raise ValueError("Missing or stale control receipts")
    output = {}
    for order in ORDERS:
        pairs = []
        for number in (1, 2, 3):
            roles = ROLES[::-1] if number == 2 else ROLES
            pair = {"pair": number, "execution_order": roles}
            for role in roles:
                pair[role] = panel(root, order, number, role)
            baseline, candidate = pair["baseline"], pair["candidate"]
            delta = {k: candidate[k] - baseline[k] for k in baseline}
            pair["signed_delta"] = delta
            pair["signed_percent"] = {k: 100 * d / baseline[k] if baseline[k] else None for k, d in delta.items()}
            pairs.append(pair)
        output[order] = {"pairs": pairs, "median_signed_delta": medians(pairs, "signed_delta"),
                         "median_signed_percent": medians(pairs, "signed_percent")}
    return output


if __name__ == "__main__":
    root = Path(sys.argv[1]).resolve()
    result = {"controls": controls, "paired": paired}[sys.argv[2]](root)
    print(json.dumps(result, indent=2, allow_nan=False))
