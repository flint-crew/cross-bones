from __future__ import annotations

import importlib.metadata

import cross_bones as m


def test_version():
    assert importlib.metadata.version("cross_bones") == m.__version__
