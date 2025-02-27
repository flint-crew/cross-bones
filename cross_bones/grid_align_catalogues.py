"""Implements a brute force grid search of catalogue A
to catalogue B
"""

from __future__ import annotations

from argparse import ArgumentParser
from dataclasses import dataclass
from pathlib import Path

from astropy.coordinates import SkyCoord
from astropy.table import Table

from cross_bones.catalogue import Catalogue, load_catalogue
from cross_bones.logging import logger


@dataclass
class FixedCatalogue:
    """Represent a catalogue that is fixed and used a reference for a grid search"""

    table: Table
    """The table loaded"""
    identifier: Path | str
    """Original path to the loaded catalogue, or string to download from vizier"""
    center: SkyCoord
    """Rough beam center derived from coordinates of componetns in catalogue"""
    sky_coords: SkyCoord | None = None
    """The positions from the table loaded"""


def load_fixed_catalogue(fixed_catalogue_1: Path | str) -> FixedCatalogue:
    """Load the catalogue that will be fixed and used as a reference
    to match to

    Args:
        fixed_catalogue_1 (Path | str): The catalogue to load. If a ``Path`` is it assumed to be on disk. If a ``str`` it is assumed to be an ``astroquery.vizier`` table.

    Raises:
        NotImplemented: Raised if ``fixed_catalogue`` is a string

    Returns:
        FixedCatalogue: The loaded catalogue and associated meta-data
    """
    if isinstance(fixed_catalogue_1, str):
        msg = f"{fixed_catalogue_1=} is a string, which means downloading from vizier. Not implemented yet."

        raise NotImplementedError(msg)

    if (
        isinstance(fixed_catalogue_1, Path)  # type: ignore[redundant-expr]
        and fixed_catalogue_1.exists()
        and fixed_catalogue_1.is_file()
    ):
        fixed_catalogue: Catalogue = load_catalogue(catalogue_path=fixed_catalogue_1)

        return FixedCatalogue(
            table=fixed_catalogue.table,
            identifier=fixed_catalogue.path,
            center=fixed_catalogue.center,
            sky_coords=fixed_catalogue.sky_coords,
        )

    msg = "Something went wrong"
    raise NotImplementedError(msg)


def grid_match_catalogues(
    catalogue_1: Path,
    fixed_catalogue_1: Path | str,
) -> None:
    logger.info(f"{catalogue_1=} will be matched to {fixed_catalogue_1=}")
    shift_catalogue: Catalogue = load_catalogue(catalogue_path=catalogue_1)
    fixed_catalogue: FixedCatalogue = load_fixed_catalogue(
        fixed_catalogue_1=fixed_catalogue_1
    )
    logger.debug(shift_catalogue)
    logger.debug(fixed_catalogue)


def get_parser() -> ArgumentParser:
    """Construct and return an argument parser"""
    parser = ArgumentParser(description="Align catalogues through a grid search")

    parser.add_argument(
        "catalogue_1", type=Path, help="The catalogue that will be matched"
    )

    parser.add_argument(
        "fixed_catalogue",
        type=str,
        help="The reference catalogue that the catalogue will be matched to",
    )

    return parser


def cli() -> None:
    parser = get_parser()

    args = parser.parse_args()

    grid_match_catalogues(
        catalogue_1=args.catalogue_1, fixed_catalogue_1=args.fixed_catalogue
    )


if __name__ == "__main__":
    cli()
