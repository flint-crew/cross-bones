# Cross-BONES 🏴‍☠️

[![Actions Status][actions-badge]][actions-link]
[![Documentation Status][rtd-badge]][rtd-link]

[![PyPI version][pypi-version]][pypi-link]
[![PyPI platforms][pypi-platforms]][pypi-link]

[![GitHub Discussion][github-discussions-badge]][github-discussions-link]

<!-- SPHINX-START -->

<!-- prettier-ignore-start -->
[actions-badge]:            https://github.com/flint-crew/cross-bones/workflows/CI/badge.svg
[actions-link]:             https://github.com/flint-crew/cross-bones/actions
[github-discussions-badge]: https://img.shields.io/static/v1?label=Discussions&message=Ask&color=blue&logo=github
[github-discussions-link]:  https://github.com/flint-crew/cross-bones/discussions
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

## Usage

### Command-line

```bash
$ cross_bones -h
# usage: cross_bones [-h] [-o OUTPUT_PREFIX] [--passes PASSES] paths [paths ...]

# Looking at per-beam shifts

# positional arguments:
#   paths                 The beam wise catalogues to examine

# options:
#   -h, --help            show this help message and exit
#   -o OUTPUT_PREFIX, --output-prefix OUTPUT_PREFIX
#                         The prefix to base outputs onto
#   --passes PASSES       Number of passes over the data should the iterative method attempt
```
