# Native batch-rank benchmark — prepared only

Baseline is the selected 16 KiB header5881233f... from R-only commit25c26fe2...;
candidate is the unchanged prototype header9f71fa92.... This replaces the
unlaunched8 KiB counter diagnostic as Root's priority; that source/history is
retained and was never run. No Rcpp/R API integration or package gain is claimed.

Use the existing original native packet only for helpers, input pins, mixed
fixture, short/empty controls and frozen13-field/order/native-parity outputs.
The original packet manifest remains bd681ee3.... Both new roles use16KiB×512
key pages, ordinary C++17/O2/g flags and the same disk/reference/1M input.
Only the candidate driver adds MSS_BATCH_RANK_BENCH to call the prototype API.
The product core is copied byte-for-byte and is not edited by this benchmark.

For each nonempty batch, the candidate clock begins BEFORE gathering initial
C/I keys, estimation, prefetch sorting/search, ordinary original row-order
harmonisation, directory clear and destruction of the caller vector/strings.
Input parsing, seeded permutation/inverse setup, result-vector allocation and
scientific output stay outside that native clock in both roles. Baseline has
no artificial gather work; its preprocessing removes the conditional block.
Both processes retain one Reference across ten100,000-row batches; only the
candidate's batch directory is cleared. Existing caches remain ordinary.

At most200,000 initial keys of length<=224 are retained per100,000-row batch.
The receipt records max key count/length and a caller estimate:
vector.capacity()*sizeof(string)+sum(string.capacity()+1). This estimate can
double count SSO bytes and excludes allocator overhead and transient constructed
keys before the length check. The frozen panel has short chr22/SNP keys; this
is not a generic caller heap bound. It is not an allocation
cap. The core's32MiB incremental arena budget is separate from this vector,
existing page/rank/pivot caches and other caller/input/result buffers. There is
no32MiB whole-method claim. Whole-child wait4 RSS captures the entire native
process; the accepted outer monitor bounds the whole owned tree independently.

The fixed recipe builds both drivers and both untimed parity binaries. Existing
format2/3 mixed fixtures, matched/shuffled short/empty controls and historical
native parity must pass before timing. Then three fresh serial AB/BA/AB pairs
run in EACH matched/shuffled order (12 workers). All scientific rows stream
through the unchanged digest checker and are not spooled as1M files. Only tiny
short-fixture rows are retained for existing controls. Each of12 panel digests
must exactly match the full frozen13-field ordered digest, including both
binary64 values and status/ordinal/ID associations. Seed/order/count receipts
must match. All ten candidate panel prefetch calls must admit; otherwise the
benchmark fails rather than timing an inactive speculation path.

recipe.json pins the actual RED, original failed GREEN, allocation diagnostic
and corrected full GREEN2 receipts. Original failures are preserved. The fixed
1GiB ROOT allowance for intervening sequential compiles is RESERVED/UNMEASURED,
not an observed highwater and not summed once per phase. Resume original
69b589 combined budget/baseline/MAX facts first; add the single reservation
conservatively within original80GB accounting. Fresh workspace measurement
counts retained sources/binaries/crash/fixture files. The runtime adapter must
still pass fresh provider/canonicalguard/writer/mount/capacity admission.

Prospective END07:30 UTC remains inside original compute07:37:49.716 and
provider09:37:49.716. The adapter keeps600s whole phase,5s kill,60s return,
3GiB groupRSS/2GiB AS/128MiB files/256MiB new growth and original186GiB workspace,
8GB root-free and cumulative80GB caps. Root must fully review source and release
one concrete phase before any Pod operation. No resource change, deletion,
automatic retry, Mac compiler/native/R work or cloud call occurred in preparation.
