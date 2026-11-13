# Boosting controls and GBM3 backend

Boosting uses
[`gbm3::gbmt()`](https://gbm-developers.github.io/gbm3/reference/gbmt.html)
for count, continuous, and survival outcomes. Existing public argument
names are translated internally.

## Details

Additional arguments passed through `...` can include `shrinkage`
(default `0.1`), `n.minobsinnode` (`10`), `bag.fraction` (`0.5`),
`train.fraction` (`1`), `var.monotone`, `mFeatures`, `keep.data`,
`class.stratify.cv`, and `fold.id`. Set `n.cores` to select the number
of GBM3 threads, or supply `par.details = gbm3::gbmParallel(...)`; do
not supply both. With neither argument, the `gbm.parallel` option is
respected and the fallback is one thread.

GBM3's learning-rate and training-sample defaults are overridden to
retain precmed's previous GBM settings. Fits use the `GBMFit` class
rather than `gbm`; use the generic
[`predict()`](https://rdrr.io/r/stats/predict.html) method for
prediction. Predictions exclude exposure offsets, preserving the
previous rate-scale behavior. Results may differ from GBM even with the
same settings and seed.
