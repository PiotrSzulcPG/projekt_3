# projekt3_test
## Installation

- Clone this repository
- `pip install ./scikit_build_example`

## Test call

```python
import scikit_build_example

scikit_build_example.add(1, 2)
```

## Deleted files from original scikit_build_example

* `docs/`: Documentation
* `.pre-commit-config.yaml`: Configuration for the fantastic static-check runner [pre-commit][].
* `noxfile.py`: Configuration for the [nox][] task runner, which helps make setup easier for contributors.

## Files considered for deletion

* `.github`: configuration for [Dependabot][] and [GitHub Actions][]
* `conda.recipe`: Example recipe. Normally you should submit projects to conda-forge instead of building them yourself, but this is useful for testing the example.

## License

pybind11 is provided under a BSD-style license that can be found in the LICENSE
file. By using, distributing, or contributing to this project, you agree to the
terms and conditions of this license.

[cibuildwheel]: https://cibuildwheel.readthedocs.io
[scientific-python development guide]: https://learn.scientific-python.org/development
[dependabot]: https://docs.github.com/en/code-security/dependabot
[github actions]: https://docs.github.com/en/actions
[pre-commit]: https://pre-commit.com
[nox]: https://nox.thea.codes
[pybind11]: https://pybind11.readthedocs.io
[scikit-build-core]: https://scikit-build-core.readthedocs.io
[scikit-build (classic)]: https://scikit-build.readthedocs.io
