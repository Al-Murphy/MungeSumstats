# Completed-scales verifier fix review

**Spec and quality pass for guarded execution.** The one-line change at `verify-completed-scales.py:164` resolves the sole Important finding. The collector now uses the authenticated helper’s `ROLE_LOGS` for each of the twelve original workers. This retains stdout and stderr coverage and adds each `r-heap.rds` to the closing hash check at lines 279–280 and to `RESULT.json.original_small_controls` at line 308. The existing bounded control reader rejects a missing or changed log. I found no new issue in the change.

I verified source SHA256 `fb61736c23dcde4bfe16fedfa450c72a67d8803273536022193d4d39004e2848` and fix diff SHA256 `12335cb92cc6640992b436e9c9283ba315e9461ce488509f1963137a1631dfcf`, inspected the frozen log list and read the appended author report. This verdict covers source only. Guarded execution, both independent oracles and observer closure remain pending.
