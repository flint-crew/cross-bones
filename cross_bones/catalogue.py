from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table
from numpy.typing import NDArray

from cross_bones.logging import logger

Paths = tuple[Path, ...]


@dataclass
class Offset:
    """Contains offsets in the RA and Dec directions in arcsec"""

    ra: float = 0.0
    """Offset in RA direction"""
    dec: float = 0.0
    """Offset in Dec direction"""


@dataclass
class Catalogue:
    """Represent a per-beam ASKAP component catalogue"""

    table: Table
    """The table loaded"""
    path: Path
    """Original path to the loaded catalogue"""
    center: SkyCoord
    """Rough beam center derived from coordinates of componetns in catalogue"""
    sky_coords: SkyCoord | None = None
    """The positions from the table loaded"""
    fixed: bool = False
    """Indicates whether beam has been fixed into a place"""
    offset: Offset = field(default_factory=Offset)
    """Per beam offsets, if known, in arcsec"""
    idx: int | None = None
    """Some optional identifier"""

    def __repr__(self) -> str:
        return f"Catalogue(idx={self.idx}, table={len(self.table)} sources, path={self.path}, fixed={self.fixed})"


Catalogues = list[Catalogue]


def make_sky_coords(table: Table | Catalogue) -> SkyCoord:
    """Create the sky-coordinates from a cataloguue table

    Args:
        table (Table | Catalogue): Loaded table or catalogue

    Returns:
        SkyCoord: Sky-positions loaded
    """
    if isinstance(table, Catalogue) and table.sky_coords:
        return table.sky_coords

    table = table.table if isinstance(table, Catalogue) else table
    return SkyCoord(table["ra"], table["dec"], unit=(u.deg, u.deg))


def estimate_skycoord_centre(
    sky_positions: SkyCoord, final_frame: str = "fk5"
) -> SkyCoord:
    """Estimate the central position of a set of positions by taking the
    mean of sky-coordinates in their XYZ geocentric frame. Quick approach
    not intended for accuracy.

    Args:
        sky_positions (SkyCoord): A set of sky positions to get the rough center of
        final_frame (str, optional): The final frame to convert the mean position to. Defaults to "fk5".

    Returns:
        SkyCoord: The rough center position
    """

    xyz_positions = sky_positions.cartesian.xyz
    xyz_mean_position = np.median(xyz_positions, axis=1)

    return SkyCoord(*xyz_mean_position, representation_type="cartesian").transform_to(
        final_frame
    )


def filter_table(table: Table) -> NDArray[bool]:
    """Filter radio components out of an aegean radio catalogue
    based on their distance to neighbouring components and compactness.

    Args:
        table (Table): Aegean radio component catalogue

    Returns:
        np.ndarray: Boolean array of components to keep.
    """
    sky_coord = SkyCoord(table["ra"], table["dec"], unit=(u.deg, u.deg))

    isolation_mask = sky_coord.match_to_catalog_sky(sky_coord, nthneighbor=2)[1] > (
        0.01 * u.deg
    )

    ratio = table["int_flux"] / table["peak_flux"]
    ratio_mask = (ratio > 0.8) & (ratio < 1.2)

    return isolation_mask & ratio_mask


def load_catalogue(catalogue_path: Path, idx: int | None = None) -> Catalogue:
    """Load a beam catalogue astropy table

    Args:
        catalogue_path (Path): Path to load catalogue from
        idx (int | None, optional): Some optional identifier added to Catalogue. Defaults to None.

    Returns:
        Catalogue: Loaded catalogue
    """
    logger.info(f"Loading {catalogue_path}")
    table = Table.read(catalogue_path)

    table_mask = filter_table(table=table)
    sub_table = table[table_mask]

    sky_coords = make_sky_coords(table=table)

    center = estimate_skycoord_centre(
        SkyCoord(table["ra"], table["dec"], unit=(u.deg, u.deg))
    )

    return Catalogue(
        table=sub_table,
        path=catalogue_path,
        sky_coords=sky_coords,
        center=center,
        idx=idx,
    )


def load_catalogues(catalogue_paths: Paths) -> Catalogues:
    """Load in all of the catalgues"""
    return [
        load_catalogue(catalogue_path=catalogue_path, idx=idx)
        for idx, catalogue_path in enumerate(catalogue_paths)
    ]


def save_catalogue_shift_positions(
    catalogues: Catalogues, output_path: Path | None = None
) -> Path:
    from pandas import DataFrame

    output_path = output_path if output_path else Path("shifts.csv")

    shifts_df = DataFrame(
        [
            {
                "path": catalogue.path,
                "d_ra": catalogue.offset.ra,
                "d_dec": catalogue.offset.dec,
            }
            for catalogue in catalogues
        ]
    )

    logger.info(f"Writing {output_path}")
    shifts_df.to_csv(output_path, index=False)

    return output_path
