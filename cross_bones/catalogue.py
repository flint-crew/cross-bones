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
class CatalogueKeys:
    """Contains catalogue column names/keys. Defaults are aegean-based."""

    ra: str = "ra"
    """RA column key."""
    dec: str = "dec"
    """DEC column key."""
    int_flux: str = "int_flux"
    """Integrated flux density column key."""
    peak_flux: str = "peak_flux"
    """Peak flux density column key."""
    local_rms: str = "local_rms"
    """Local rms column key."""


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
    table_keys: CatalogueKeys = field(default_factory=CatalogueKeys)
    """Column names/keys with defaults from aegean component lists."""

    def __repr__(self) -> str:
        return f"Catalogue(idx={self.idx}, table={len(self.table)} sources, path={self.path}, fixed={self.fixed})"


Catalogues = list[Catalogue]


# TODO: Expand this to include other best guesses and units
def _guess_table_keys(table: Table | None = None) -> dict[str, str]:
    """Get table keys on a best effort basis"""

    known_ra_dec = (
        ("ra", "dec"),
        ("RAJ2000", "DEJ2000"),
        ("RAJ2000", "DECJ2000"),
        ("RA", "DEC"),
    )

    if table is not None:
        logger.warn("Stub function is ignoring the table input")
        for ra, dec in known_ra_dec:
            if ra in table.colnames and dec in table.colnames:
                # for the moment this is best effort. Need something more robust
                return {
                    "ra": ra,
                    "dec": dec,
                    "int_flux": "int_flux",
                    "peak_flux": "peak_flux",
                    "local_rms": "local_rms",
                }

    return {
        "ra": "ra",
        "dec": "dec",
        "int_flux": "int_flux",
        "peak_flux": "peak_flux",
        "local_rms": "local_rms",
    }


def make_sky_coords(
    table: Table | Catalogue, ra_key: str = "ra", dec_key: str = "dec"
) -> SkyCoord:
    """Create the sky-coordinates from a cataloguue table

    Args:
        table (Table | Catalogue): Loaded table or catalogue
        ra_key (str): RA column key name
        dec_key (str): DEC column key name

    Returns:
        SkyCoord: Sky-positions loaded
    """
    if isinstance(table, Catalogue) and table.sky_coords:
        return table.sky_coords

    table = table.table if isinstance(table, Catalogue) else table
    return SkyCoord(table[ra_key], table[dec_key], unit=(u.deg, u.deg))


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


def filter_table(
    table: Table,
    min_snr: float = 10.0,
    min_iso: float = 36.0,
    table_keys: dict[str, str] | None = None,
) -> NDArray[bool]:
    """Filter radio components out of a radio catalogue
    based on their distance to neighbouring components, compactness,
    and optionally minimum SNR.

    Args:
        table (Table): Aegean radio component catalogue
        min_snr (float): Minimum SNR of component. Default 10.
        min_iso (float): Minimum separation from a neighbour in arcsec. Default 30.
        table_keys (dict[str, str] | None, optional): Dictionary of column keys. If None is is assumed to be aegean. Defaults to None.

    Returns:
        np.ndarray: Boolean array of components to keep.
    """

    if table_keys is None:
        table_keys = _guess_table_keys()

    sky_coord = SkyCoord(
        ra=table[table_keys["ra"]], dec=table[table_keys["dec"]], unit=(u.deg, u.deg)
    )

    isolation_mask = sky_coord.match_to_catalog_sky(sky_coord, nthneighbor=2)[1] > (
        min_iso * u.arcsec
    )

    ratio = table[table_keys["int_flux"]] / table[table_keys["peak_flux"]]
    ratio_mask = (ratio > 0.8) & (ratio < 1.2)

    snr = table[table_keys["peak_flux"]] / table[table_keys["local_rms"]]
    snr_mask = snr > min_snr

    return isolation_mask & ratio_mask & snr_mask


def load_catalogue(
    catalogue_path: Path,
    idx: int | None = None,
    min_snr: float = 10.0,
    min_iso: float = 36.0,
    table_keys: dict[str, str] | None = None,
) -> Catalogue:
    """Load a beam catalogue astropy table

    Args:
        catalogue_path (Path): Path to load catalogue from
        idx (int | None, optional): Some optional identifier added to Catalogue. Defaults to None.
        min_snr (float): Minimum SNR of component. Default 10.
        min_iso (float): Minimum separation from a neighbour in arcsec. Default 30.
        table_keys (dict[str, str] | None, optional): Dictionary of column keys. If None is is assumed to be aegean. Defaults to None.

    Returns:
        Catalogue: Loaded catalogue
    """
    if table_keys is None:
        table_keys = _guess_table_keys()

    logger.info(f"Loading {catalogue_path}")
    table = Table.read(catalogue_path)

    table_mask = filter_table(
        table=table, table_keys=table_keys, min_snr=min_snr, min_iso=min_iso
    )
    sub_table = table[table_mask]

    sky_coords = make_sky_coords(
        table=table, ra_key=table_keys["ra"], dec_key=table_keys["dec"]
    )

    center = estimate_skycoord_centre(
        SkyCoord(
            ra=table[table_keys["ra"]],
            dec=table[table_keys["dec"]],
            unit=(u.deg, u.deg),
        )
    )

    catalogue_keys = CatalogueKeys(
        ra=table_keys["ra"],
        dec=table_keys["dec"],
        int_flux=table_keys["int_flux"],
        peak_flux=table_keys["peak_flux"],
        local_rms=table_keys["local_rms"],
    )

    return Catalogue(
        table=sub_table,
        path=catalogue_path,
        sky_coords=sky_coords,
        center=center,
        idx=idx,
        table_keys=catalogue_keys,
    )


def load_catalogues(
    catalogue_paths: Paths,
    min_snr: float = 10.0,
    min_iso: float = 36.0,
    table_keys: dict[str, str] | None = None,
) -> Catalogues:
    """Load in all of the catalgues"""

    return [
        load_catalogue(
            catalogue_path=catalogue_path,
            idx=idx,
            table_keys=table_keys,
            min_snr=min_snr,
            min_iso=min_iso,
        )
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
