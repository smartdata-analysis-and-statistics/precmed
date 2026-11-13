# Generate K-fold Indices for Cross-Validation

This function generates indices for K-fold cross-validation based on the
total sample size \`N\` and the number of folds \`Kfold\`. If \`reverse
= TRUE\`, the remainder indices will be assigned in reverse order.

## Usage

``` r
generate_kfold_indices(N, Kfold, reverse = FALSE)
```

## Arguments

- N:

  Integer. Total sample size (number of observations).

- Kfold:

  Integer. The number of folds to split the data into.

- reverse:

  Logical. Whether to reverse the remainder indices when \`N\` is not
  divisible by \`Kfold\`. Defaults to \`FALSE\`.

## Value

A vector of length \`N\` containing the fold assignments (from 1 to
\`Kfold\`).

## Author

Thomas Debray
