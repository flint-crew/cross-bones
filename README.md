# Cross-BONES 🏴‍☠️

[![Actions Status][actions-badge]][actions-link]
[![Documentation Status][rtd-badge]][rtd-link]

[![PyPI version][pypi-version]][pypi-link]
[![PyPI platforms][pypi-platforms]][pypi-link]

<!-- SPHINX-START -->

<!-- prettier-ignore-start -->
[actions-badge]:            https://github.com/flint-crew/cross-bones/workflows/CI/badge.svg
[actions-link]:             https://github.com/flint-crew/cross-bones/actions
[pypi-link]:                https://pypi.org/project/cross-bones/
[pypi-platforms]:           https://img.shields.io/pypi/pyversions/cross-bones
[pypi-version]:             https://img.shields.io/pypi/v/cross-bones
[rtd-badge]:                https://readthedocs.org/projects/cross-bones/badge/?version=latest
[rtd-link]:                 https://cross-bones.readthedocs.io/en/latest/?badge=latest

<!-- prettier-ignore-end -->

🏴‍☠️ Cross-match By Offsetting Neighbouring Extracted Sources 🏴‍☠️

Attempt to align catalogues of radio sources onto a self-consistent grid. It
implements a fairly simple iterative procedure that aims to reduce separations
between sources in common between pairs of catalogues.

## Installation

All dependencies can be installed via `pip`. We strongly recommend using [uv](https://docs.astral.sh/uv/).


### Install from GitHub (latest)

```bash
pip install git+https://github.com/flint-crew/cross-bones
```

### Install from PyPI release (stable)

```bash
pip install cross-bones
```

## Usage and documentation

### Command-line

Full documentation is hosted on [ReadTheDocs](https://cross-bones.readthedocs.io/en/latest/).

## Contributing

Contributions are welcome! Please open an issue to discuss ahead of opening a pull request.

### Dev tools

The dev tooling can be installed by using the `dev` option flags e.g.

```bash
git clone https://github.com/flint-crew/cross-bones.git
cd cross-bones
pip install -e .[dev]
pre-commit install
```
