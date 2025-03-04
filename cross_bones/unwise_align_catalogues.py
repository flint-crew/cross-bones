from __future__ import annotations

from argparse import ArgumentParser
from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.coordinates import (
    SkyCoord,
    concatenate,
    match_coordinates_sky,
)
from astropy.table import Table, unique, vstack
from astroquery import vizier
from numpy.typing import NDArray
from typing_extensions import TypeAlias

from cross_bones.catalogue import (
    Catalogue,
    Catalogues,
    TableKeys,
    load_catalogues,
)
from cross_bones.logging import logger
from cross_bones.matching import (
    OffsetGridSpace,
    find_minimum_offset_space,
)
from cross_bones.plotting import plot_offset_grid_space, plot_offsets_in_field

Paths = tuple[Path, ...]
MatchMatrix: TypeAlias = NDArray[int]


def guess_sbid_and_field_racs(
    catalogue_path: str | Path,
) -> tuple[int, str]:
    """Attempt to extract the SBID and field name from a file path.
    The filenames assumes some deliminted name scheme, with '.' as
    the field marker.

    The field name is the second item, and it taken as it. The SBID
    is taken as the first field, and requires a 'SB' prefix.

    Args:
        catalogue_path (str | Path): The file path to extract

    Raises:
        RuntimeError: Raised when a path does not follow expectations

    Returns:
        tuple[int, str]: And SBID and field name
    """
    name_components = str(Path(catalogue_path).name).split(".")
    logger.info(f"These are the split name components: {name_components}")
    field_name = name_components[1]

    if name_components[0][:2] != "SB":
        msg = f"SBID is not found in {name_components=}"
        raise RuntimeError(msg)

    sbid_str = name_components[0][2:]
    try:
        sbid = int(sbid_str)
    except ValueError as err:
        msg = f"Can not convert {sbid_str=} to an int"
        raise RuntimeError(msg) from err

    return sbid, field_name


def _get_output_table_path(
    output_dir: Path | str, output_name: str, name_prefix: str | None = None
) -> Path:
    """Make the output table download path and perform some basic checks"""
    base_output_dir: Path = Path(output_dir)

    if not base_output_dir.exists():
        logger.info(f"Creating {base_output_dir=}")
        base_output_dir.mkdir(parents=True, exist_ok=True)
    else:
        assert base_output_dir.is_dir(), (
            f"{base_output_dir} already exists and is not a directory"
        )

    # The name prefix may be the name of the survery being downloaded etc
    out_name = [output_name]
    if name_prefix:
        out_name.insert(0, name_prefix)

    return Path(f"{output_dir!s}/{'_'.join(out_name)}.fits")


def _download_vizer_id_to_table(
    sky_coord_center: SkyCoord, vizier_id: str, radius_deg: float, max_retries: int = 3
) -> Table:
    """Internal method to download a Vizer catalogue around a region

    Args:
        sky_coord_center (SkyCoord): The position to center the region on
        vizier_id (str): The Vizer ID of the catalogue to download
        radius_deg (float): Radius to obtain data through
        max_retries (int, optional): The maximum number of attempts to download a table before a RunTimeError is raised

    Raises:
        RuntimeError: Raised when multiple attempts to download the table has failed

    Returns:
        Table: The downloaded table
    """
    attempt = 1
    while attempt <= max_retries:
        try:
            logger.info(f"{attempt=} to download {vizier_id=}")
            vizier.Vizier.ROW_LIMIT = -1
            vizier_tables = vizier.Vizier(
                columns=["RAJ2000", "DEJ2000"], row_limit=-1, timeout=720
            ).query_region(
                sky_coord_center, catalog=vizier_id, width=radius_deg * u.deg
            )
            break

        except (ConnectionAbortedError, ConnectionError, ValueError) as e:
            logger.warning(f"Caught {e=}")
            logger.warning("Retry VizieR in 120s")
            from time import sleep

            sleep(120)
            attempt += 1
            continue
    else:
        msg = f"Too many attempts to download vizier table. {max_retries=}"
        raise RuntimeError(msg)

    assert len(vizier_tables) == 1, "More tables downloaded than expected"

    # We need to clean the table up so it can be written nicely
    vizier_table = vizier_tables[0]
    for icol, _ in enumerate(vizier_table.itercols()):
        vizier_table.columns[icol].description = ""
        vizier_table.meta["description"] = ""

    return vizier_table


def download_vizier_catalogue(
    field_name: str,
    beam_skycoords: list[SkyCoord],
    radius_deg: float = 5.0,
    unwise_table_location: str | Path = "./",
    vizier_id: str = "II/363/unwise",
) -> Table:
    """Download tables from vizier given a catalogue ID and positions to query a region around.
    The tables will be concatenated tpgether, removed of duplicates and returned.

    Args:
        field_name (str): The field name of the region being downloaded.
        beam_skycoords (list[SkyCoord]): A collection of coordinates defining a region
        radius_deg (float, optional): The size of the region to download. Defaults to 5.0.
        unwise_table_location (str | Path, optional): Location of the path of the output table. Defaults to "./".
        vizier_id (str, optional): The vizier catalogue ID to download. Defaults to "II/363/unwise".

    Returns:
        Table: _description_
    """
    download_table_path = _get_output_table_path(
        output_dir=unwise_table_location, output_name=field_name
    )

    if download_table_path.exists() and download_table_path.is_file():
        logger.info(f"Found existing {download_table_path=}, loading.")
        field_cat = Table.read(download_table_path)
        logger.debug(f"{len(field_cat)} rows loaded")

        return field_cat

    beam_catalogues: list[Table] = []
    for beam_idx, beam_skycoord in enumerate(beam_skycoords):
        logger.debug(
            f"Downloading {vizier_id} table for field {field_name} and beam {beam_idx}"
        )

        beam_cat_path = _get_output_table_path(
            output_dir=unwise_table_location, output_name=f"{field_name}_{beam_idx:02d}"
        )

        if beam_cat_path.exists():
            logger.info(f"Found cached file, reading {beam_cat_path}")
            beam_cat = Table.read(beam_cat_path)
        else:
            beam_cat = _download_vizer_id_to_table(
                sky_coord_center=beam_skycoord,
                vizier_id=vizier_id,
                radius_deg=radius_deg,
            )
            beam_cat.write(beam_cat_path, overwrite=False)

        logger.info(f"{len(beam_cat)} rows downloaded")
        logger.info("Merging into field catalogue")
        beam_catalogues.append(beam_cat)
        logger.info(f"Unlinking {beam_cat_path}")
        beam_cat_path.unlink()

    field_catalogue: Table = unique(vstack(beam_catalogues)) if field_cat else beam_cat
    field_catalogue.write(download_table_path, format="fits", overwrite=True)

    return field_catalogue


def add_offset_to_coords_skyframeoffset(
    sky_coords: SkyCoord, offset: tuple[float, float]
) -> SkyCoord:
    """Add angular offsets to a ``SkyCoord`` object.

    Args:
        sky_coords (SkyCoord): The input positions that will be shifted
        offset (tuple[float, float]): The delta RA and Dec units

    Returns:
        SkyCoord: The shifted units
    """
    d_ra = -offset[0] * u.arcsec
    d_dec = -offset[1] * u.arcsec

    return sky_coords.spherical_offsets_by(d_ra, d_dec)


def get_offset_space(
    catalogue: Catalogue,
    unwise_table: Table,
    window: tuple[float, float, float, float, float],
    beam: int | None = None,
) -> OffsetGridSpace:
    # TODO create skycoord object earlier, and pass between beams
    unwise_sky = SkyCoord(
        ra=unwise_table["RAJ2000"], dec=unwise_table["DEJ2000"], unit=(u.deg, u.deg)
    )

    shifted_table = catalogue.table.copy()
    cata_sky = SkyCoord(
        ra=shifted_table[catalogue.table_keys.ra],
        dec=shifted_table[catalogue.table_keys.dec],
        unit=(u.deg, u.deg),
    )

    # The range of RA and Dec searched
    ra_extent = np.abs(window[1] - window[0])
    dec_extent = np.abs(window[3] - window[2])
    # The number of bins in RA and Dec
    ra_bins = int(ra_extent / window[4])
    dec_bins = int(dec_extent / window[4])

    ras = np.linspace(window[0], window[1], ra_bins)
    decs = np.linspace(window[2], window[3], dec_bins)

    coords: list[SkyCoord] = []

    n_sources = len(cata_sky)
    n_delta = n_sources * len(decs) * len(ras)

    broadcast_d_ra = np.zeros(n_delta)
    broadcast_d_dec = np.zeros(n_delta)

    for d_idx, dec in enumerate(decs):
        for r_idx, ra in enumerate(ras):
            logger.debug(f"{d_idx=} {r_idx=}")
            i = len(coords)
            j = i + 1
            broadcast_d_ra[i * n_sources : j * n_sources] = ra
            broadcast_d_dec[i * n_sources : j * n_sources] = dec

            coords.append(cata_sky)
    collected_coords = concatenate(coords)

    shifted_sky = add_offset_to_coords_skyframeoffset(
        collected_coords, (broadcast_d_ra, broadcast_d_dec)
    )

    matches = match_coordinates_sky(
        shifted_sky, unwise_sky, nthneighbor=1, storekdtree=True
    )

    accumulated_seps = np.zeros((dec_bins, ra_bins))

    results = {}
    pair_decra = []

    for d_idx, dec in enumerate(decs):
        for r_idx, ra in enumerate(ras):
            start_idx = d_idx * len(ras) + r_idx
            end_idx = start_idx + 1
            k = (d_idx, r_idx)
            v = (
                np.sum(matches[1][start_idx * n_sources : end_idx * n_sources].value)
                / n_sources
            )

            seps = (dec, ra, v)
            pair_decra.append((dec, ra))
            results[k] = seps
            accumulated_seps[k] = v

    array_decra = np.array(pair_decra)

    return OffsetGridSpace(
        dec_offsets=array_decra[:, 0],
        ra_offsets=array_decra[:, 1],
        beam=beam if beam else 0,
        n_sources=n_sources,
        seps=accumulated_seps,
    )


def unwise_shifts(
    catalogue_paths: Paths,
    table_keys: TableKeys,
    output_prefix: str | None = None,
    sbid: int | None = None,
    field_name: str | None = None,
    beam_table: str = "closepack36_beams.fits",
    unwise_table_location: Path = Path("./"),
    min_snr: float = 10.0,
    min_iso: float = 36.0,
) -> Path:
    if output_prefix:
        output_parent = Path(output_prefix).parent
        output_parent.mkdir(exist_ok=True, parents=True)

    field_beams = Table.read(beam_table)

    # assuming things are named correctly.
    _catalogue_paths = [Path(path) for path in catalogue_paths]
    _catalogue_paths.sort()
    catalogue_paths = tuple(_catalogue_paths)

    # if SBID and field name not provide, guess from first input catalogue:
    if sbid is None or field_name is None:
        sbid, field_name = guess_sbid_and_field_racs(catalogue_path=catalogue_paths[0])

    beam_inf = field_beams[np.where(field_beams["FIELD_NAME"] == field_name)[0]]
    beam_skycoords = SkyCoord(
        ra=beam_inf["RA_DEG"] * u.deg, dec=beam_inf["DEC_DEG"] * u.deg, frame="fk5"
    )

    logger.info(f"Will be processing {len(catalogue_paths)} catalogues")
    catalogues: Catalogues = load_catalogues(
        catalogue_paths=catalogue_paths,
        table_keys=table_keys,
        min_iso=min_iso,
        min_snr=min_snr,
    )

    unwise_field_cat = download_vizier_catalogue(
        field_name=field_name,
        beam_skycoords=beam_skycoords,
        unwise_table_location=unwise_table_location,
    )

    window_incs = [10.0, 1.0, 0.1]
    windows = [(wi, wi, wi, wi, wi / 10.0) for wi in window_incs]

    # TODO add other stats??
    ra_offsets = np.full((len(catalogues),), np.nan)
    dec_offsets = np.full((len(catalogues),), np.nan)

    all_offset_results = []

    for beam in range(36):
        logger.debug(f"Working on beam {beam}")

        min_ra, min_dec = 0.0, 0.0

        for i in range(len(windows)):
            window = (
                min_ra - windows[i][0],
                min_ra + windows[i][1],
                min_dec - windows[i][2],
                min_dec + windows[i][3],
                windows[i][4],
            )

            logger.debug(f"Working on window: {window}")

            offset_results: OffsetGridSpace = get_offset_space(
                catalogue=catalogues[beam],
                unwise_table=unwise_field_cat,
                window=window,
                beam=beam,
            )

            min_ra, min_dec, min_sep = find_minimum_offset_space(offset_results)

        # per window?
        plot_offset_grid_space(
            f"SB{sbid}.{field_name}_beam{beam:02d}.offset_grid.png",
            offset_results,
            window=window,
        )

        ra_offsets[beam] = min_ra
        dec_offsets[beam] = min_dec

        all_offset_results.append(offset_results)

    plot_offsets_in_field(
        offset_results=all_offset_results,
        fname=f"SB{sbid}.{field_name}_offset_grid.png",
    )

    shift_table = Table(
        [catalogue_paths, ra_offsets, dec_offsets], names=["path", "d_ra", "d_dec"]
    )

    if output_prefix is None:
        outname = "SB{sbid}.{field_name}-unwise-shifts.csv"
    else:
        outname = output_prefix + "-unwise-shifts.csv"

    output_path = Path(outname)
    shift_table.write(output_path, format="ascii.csv", overwrite=True)

    return output_path


def get_parser() -> ArgumentParser:
    parser = ArgumentParser(description="Looking at per-beam shifts")

    parser.add_argument(
        "paths", nargs="+", type=Path, help="The beam wise catalogues to examine"
    )

    parser.add_argument(
        "-s",
        "--sbid",
        type=int,
        default=None,
        help="ASKAP SBID, if none provided it will be guessed from the filenames.",
    )

    parser.add_argument(
        "-f",
        "--field-name",
        type=str,
        default=None,
        help="ASKAP field name, if none provided it will assumed a RACS observation and guessed from the filename.",
    )

    parser.add_argument(
        "--beam-table",
        default="closepack36_beams.fits",
        type=str,
        help="Table of beam positions for a given ASKAP footprint. Default 'closepack36_beams.fits",
    )

    parser.add_argument(
        "--unwise-table-location",
        default=Path("./"),
        type=Path,
        help="Directory of the unWISE tables. Default './'",
    )

    parser.add_argument(
        "-o",
        "--output-prefix",
        type=str,
        default=None,
        help="The prefix to base outputs onto",
    )
    parser.add_argument(
        "--coord-keys",
        nargs=2,
        default=["ra", "dec"],
        type=str,
        help="Column names/keys in tables for (ra, dec). [Default ('ra', 'dec')]",
    )
    parser.add_argument(
        "--flux-keys",
        nargs=2,
        default=["int_flux", "peak_flux"],
        type=str,
        help="Column names/keys in tables for integrated and peak flux density. [Default ('int_flux', 'peak_flux')]",
    )
    parser.add_argument(
        "--rms-key",
        default="local_rms",
        type=str,
        help="Local rms column name/key in tables. Default 'local_rms'",
    )

    parser.add_argument(
        "--snr-min", default=5.0, type=float, help="Minimum SNR of sources. Default 5"
    )

    parser.add_argument(
        "--iso-min",
        default=36.0,
        type=float,
        help="Minimum separation between close neighbours in arcsec. Default 36",
    )

    return parser


def cli() -> None:
    parser = get_parser()

    args = parser.parse_args()

    table_keys = TableKeys(
        ra=args.coord_keys[0],
        dec=args.coord_keys[1],
        int_flux=args.flux_keys[0],
        peak_flux=args.flux_keys[1],
        local_rms=args.rms_key,
    )

    unwise_shifts(
        catalogue_paths=args.paths,
        table_keys=table_keys,
        output_prefix=args.output_prefix,
        sbid=args.sbid,
        field_name=args.field_name,
        beam_table=args.beam_table,
        unwise_table_location=args.unwise_table_location,
        min_snr=args.snr_min,
        min_iso=args.iso_min,
    )


if __name__ == "__main__":
    cli()
