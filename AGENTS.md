# Repository guidance

- This repository contains the `precmed` R package. Follow existing package conventions and keep documentation and tests aligned with changes.
- Maintain `NEWS.md` for user-visible changes since the latest CRAN release. Review the development entries against that release, retain relevant unreleased changes, and distinguish major changes, minor changes, and bug fixes in separate sections. Verify the latest CRAN release when establishing the baseline; avoid rewriting historical release notes unnecessarily.
- Preserve existing behavior where practical. Document changes to function arguments, defaults, returned objects, or results in NEWS and relevant documentation.
- For changes to statistical methods or C backends, test against the R implementation or established results where applicable, allowing for reasonable floating-point differences.
- Add or update tests for behavior changes, including relevant edge cases such as missing treatment arms, small samples, and survival horizons. Use fixed seeds in stochastic tests and examples, and avoid relying on exact random results across platforms.
- Update roxygen documentation and regenerate help files rather than editing generated help files by hand. Avoid committing compiled artifacts such as `.o` and `.so` files.
- Run tests appropriate to the change and, when practical, check the package with `R CMD check`.
