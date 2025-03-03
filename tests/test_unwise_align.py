"""Testing functions used throughout the external (e.g. unwise) align grid search"""

from __future__ import annotations

from pathlib import Path

import pytest

from cross_bones.unwise_align_catalogues import guess_sbid_and_field_racs


def test_guess_sbid_and_field_racs():
    """Attempt to get the SBID and field name out from a catalogue path"""

    example = "SB1234.RACS_1234-345.round3.comp_fits"
    example_path = Path(example)
    for ex in (example, example_path):
        sbid, field = guess_sbid_and_field_racs(catalogue_path=ex)
        assert sbid == 1234
        assert field == "RACS_1234-345"

    with pytest.raises(RuntimeError):
        _ = guess_sbid_and_field_racs(catalogue_path="RACS_1234-345.round3.comp_fits")

    with pytest.raises(RuntimeError):
        _ = guess_sbid_and_field_racs(catalogue_path="SB1234aa.round3.comp_fits")
