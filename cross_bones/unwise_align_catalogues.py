#! /usr/bin/env python

from __future__ import annotations

import os

from argparse import ArgumentParser
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import (
    SkyCoord,
    search_around_sky,
    concatenate,
    match_coordinates_sky
)
from astropy.table import Table, unique, vstack
from astroquery import vizier
from time import time

from numpy.typing import NDArray
from typing_extensions import TypeAlias

from cross_bones.catalogue import (
    Catalogue,
    Catalogues,
    Offset,
    load_catalogues,
    make_sky_coords,
    save_catalogue_shift_positions,
)
from cross_bones.logging import logger
from cross_bones.matching import (
    Match, 
    calculate_matches,
    OffsetGridSpace,
    find_minimum_offset_space
)
from cross_bones.plotting import plot_offset_grid_space, plot_offsets_in_field

Paths = tuple[Path, ...]
MatchMatrix: TypeAlias = NDArray[int]


def guess_sbid_and_field_racs(
    catalogue_path: str,
) -> tuple[int, str]:
    """Attempt to guess SBID and field name from a catalogue name.
    
    TODO: generalise for other ASKAP surveys"""

    catalogue_path = os.path.basename(catalogue_path)
    field_name_pos = catalogue_path.find("RACS_")
    if field_name_pos == -1:
        raise RuntimeError("No field name found in {}".format(catalogue_path))
    field_name = catalogue_path[field_name_pos:field_name_pos+12]

    sbid_pos = catalogue_path.find("SB")
    if sbid_pos == -1:
        raise RuntimeError("No SB found in {}".format(catalogue_path))
    sbid = int(catalogue_path[sbid_pos+2:].split(".")[0])

    return sbid, field_name


def download_unwise(
    field_name: str, 
    beam_sc, 
    radius_deg: float = 5.,
    unwise_table_location: str = "./",
    unwise_vizier: str = "II/363/unwise") -> Table:

    unwise_table = f"{unwise_table_location}/unwise_{field_name}.fits"
    
    field_cat = None

    if os.path.exists(unwise_table):

        field_cat = Table.read(unwise_table)
        logger.debug(f"{len(field_cat)} rows loaded")
        
        return field_cat
    
    for beam in range(36):
        
        logger.debug(f"Downloading {unwise_vizier} table for field {field_name} and beam {beam}")
        
        unwise_beam_table = f"{unwise_table_location}/unwise_{field_name}_beam{beam:02d}.fits"
        
        if os.path.exists(unwise_beam_table):
            
            # try:
            #     logger.warning("Removing astroquery cache")
            #     os.system("rm -fr {}/.astropy/cache/astroquery/Vizier".format(os.environ["HOME"]))
            # except Exception:
            #     # check what except might occur
            #     pass
            while True:
                try:

                    vizier.Vizier.ROW_LIMIT = -1
                    unwise_tables = vizier.Vizier(
                        columns=["RAJ2000", "DEJ2000"], 
                        row_limit=-1, 
                        timeout=720
                    ).query_region(
                        beam_sc[beam], 
                        catalog=unwise_vizier, 
                        width=radius_deg*u.deg
                    )

                except:
                    # possible remove
                    logger.warning("Retry VizieR in 120s")
                    time.sleep(120)
                    continue
        
                break

            assert len(unwise_tables) == 1
            
            beam_cat = unwise_tables[0]
            
            for icol, column in enumerate(beam_cat.itercols()):
                beam_cat.columns[icol].description = ""
                beam_cat.meta["description"] = ""

            beam_cat.write(unwise_beam_table, overwrite=False)

        else:
            logger.info(f"Found cached file, reading {unwise_beam_table}")
            beam_cat = Table.read(unwise_beam_table)

        logger.info(f"{len(beam_cat)} rows downloaded")
        logger.info("Merging into field catalogue")
        if beam == 0:
            field_cat = beam_cat
        else:
            field_cat = unique(vstack([field_cat, beam_cat]))

    field_cat.write(unwise_table, format="fits", overwrite=True)
    logger.info("Cleaning cached beam catalogues")
    os.system(f"rm -fr {unwise_table_location}/unwise_{field_name}_beam*.fits")
    
    return field_cat



def add_offset_to_coords_skyframeoffset(sky_coords, offset):
    d_ra = -offset[0]*u.arcsec
    d_dec = -offset[1]*u.arcsec
    
    new_coords = sky_coords.spherical_offsets_by(d_ra, d_dec)
    return new_coords



def get_offset_space(
    catalogue: Catalogue, 
    unwise_table: Table, 
    window: tuple[float], 
    plot: bool = True, 
    radius_deg: float = 5.,
    beam : int | None = None) -> OffsetGridSpace:

    # TODO create skycoord object earlier, and pass between beams
    unwise_sky = SkyCoord(
        ra=unwise_table["RAJ2000"], 
        dec=unwise_table["DEJ2000"], 
        unit=(u.deg, u.deg)
    )
    
    shifted_table = catalogue.table.copy()
    cata_sky = SkyCoord(
        ra=shifted_table[catalogue.table_keys.ra], 
        dec=shifted_table[catalogue.table_keys.dec], 
        unit=(u.deg, u.deg)
    )

    # The range of RA and Dec searched
    ra_extent = np.abs(window[1] - window[0])
    dec_extent = np.abs(window[3] - window[2])
    # The number of bins in RA and Dec
    ra_bins = int(ra_extent / window[4])
    dec_bins = int(dec_extent / window[4])

    ras = np.linspace(window[0], window[1], ra_bins)
    decs = np.linspace(window[2], window[3], dec_bins)

    coords = []
    
    n_sources = len(cata_sky)
    n_delta = n_sources * len(decs) * len(ras)
    
    broadcast_d_ra = np.zeros(n_delta)
    broadcast_d_dec = np.zeros(n_delta)
    
    for d_idx, dec in enumerate(decs):
        for r_idx, ra in enumerate(ras):
            i = len(coords)
            j = i + 1
            broadcast_d_ra[i * n_sources:j * n_sources] = ra
            broadcast_d_dec[i * n_sources:j * n_sources] = dec
            
            coords.append(cata_sky)
    collected_coords = concatenate(coords)
    
    shifted_sky = add_offset_to_coords_skyframeoffset(collected_coords, (broadcast_d_ra, broadcast_d_dec))
    
    matches = match_coordinates_sky(shifted_sky, unwise_sky, 
        nthneighbor=1, 
        storekdtree=True
    )
    
    accumulated_seps = np.zeros((dec_bins, ra_bins))
    
    results = {}
    minimum_key = None
    minimum_seps = None
    pair_decra = []
    
    for d_idx, dec in enumerate(decs):
        for r_idx, ra in enumerate(ras):
            start_idx = (d_idx*len(ras) + r_idx)
            end_idx = start_idx + 1
            k = (d_idx, r_idx)
            v = np.sum(matches[1][start_idx*n_sources:end_idx*n_sources].value) / n_sources
    
            seps = (dec, ra, v)
            pair_decra.append((dec, ra))
            results[k] = seps
            accumulated_seps[k] = v

    array_decra = np.array(pair_decra)
    
    return OffsetGridSpace(
        dec_offsets=array_decra[:, 0], 
        ra_offsets=array_decra[:, 1], 
        beam=beam, 
        n_sources=n_sources, 
        seps=accumulated_seps
    )


def unwise_shifts(
    catalogue_paths: list[str],
    output_prefix: str | None = None,
    sbid: int | None = None,
    field_name: str | None = None,
    beam_table: str = "closepack36_beams.fits",
    field_table: str = "closepack_fields.fits",
    unwise_table_location: str = "./",
    min_snr: float = 10.,
    min_iso: float = 36.,
    do_plot: bool = True,
    table_keys: dict = {
        "ra": "ra",
        "dec": "dec",
        "int_flux": "int_flux",
        "peak_flux": "peak_flux",
        "local_rms": "local_rms"
    }
) -> Catalogues:

    if output_prefix:
        output_parent = Path(output_prefix).parent
        output_parent.mkdir(exist_ok=True, parents=True)


    field_beams = Table.read(beam_table)
    field_centres = Table.read(field_table)

    # assuming things are named correctly.
    catalogue_paths.sort()

    # if SBID and field name not provide, guess from first input catalogue:
    sbid, field_name = guess_sbid_and_field_racs(
        catalogue_path=catalogue_paths[0]
    )

    beam_inf = field_beams[np.where(field_beams["FIELD_NAME"] == field_name)[0]]
    beam_sc = SkyCoord(
        ra=beam_inf["RA_DEG"]*u.deg,
        dec=beam_inf["DEC_DEG"]*u.deg,
        frame="fk5"
    )

    logger.info(f"Will be processing {len(catalogue_paths)} catalogues")
    catalogues: Catalogues = load_catalogues(
        catalogue_paths=catalogue_paths,
        table_keys=table_keys,
        min_iso=min_iso,
        min_snr=min_snr
    )

    unwise_field_cat = download_unwise(
        field_name=field_name, 
        beam_sc=beam_sc,
        unwise_table_location=unwise_table_location
    )

    
    window_incs = [10., 1., 0.1]
    windows = [(wi, wi, wi, wi, wi/10.) for wi in window_incs]


    # TODO add other stats??
    ra_offsets = np.full((len(catalogues),), np.nan)
    dec_offsets = np.full((len(catalogues),), np.nan)

    all_offset_results = []

    for beam in range(36):
        
        logger.debug("Working on beam {}".format(beam))

        min_ra, min_dec = 0., 0.

        for i in range(len(windows)):

            window = (
                min_ra-windows[i][0],
                min_ra+windows[i][1],
                min_dec-windows[i][2],
                min_dec+windows[i][3],
                windows[i][4]
            )

            logger.debug("Working on window: {}".format(window))

            offset_results = get_offset_space(
                catalogue=catalogues[beam], 
                unwise_table=unwise_field_cat, 
                window=window, 
                plot=do_plot, 
                radius_deg=5.,
                beam=beam
            )

            all_offset_results.append(offset_results)

            min_ra, min_dec, min_sep = find_minimum_offset_space(offset_results)

        # per window?
        plot_offset_grid_space(f"SB{sbid}.{field_name}_beam{beam:02d}.offset_grid.png", 
            offset_results, 
            window=window
        )

        ra_offsets[beam] = min_ra
        dec_offsets[beam] = min_dec
    
    plot_offsets_in_field(
        offset_results=all_offset_results, 
        fname=f"SB{sbid}.{field_name}_offset_grid.png"
    )


    shift_table = Table([
        catalogue_paths, ra_offsets, dec_offsets
    ], names=[
        "path", "d_ra", "d_dec"
    ])

    if output_prefix is None:
        outname = "SB{sbid}.{field_name}-unwise-shifts.csv"
    else:
        outname = output_prefix + "-unwise-shifts.csv"

    shift_table.write(outname, format="ascii.csv", overwrite=True)
            


def get_parser() -> ArgumentParser:
    parser = ArgumentParser(description="Looking at per-beam shifts")

    parser.add_argument(
        "paths", nargs="+", type=Path, help="The beam wise catalogues to examine"
    )

    parser.add_argument(
        "-s", "--sbid",
        type=int,
        default=None,
        help="ASKAP SBID, if none provided it will be guessed from the filenames."
    )

    parser.add_argument(
        "-f", "--field_name",
        type=str,
        default=None,
        help="ASKAP field name, if none provided it will assumed a RACS observation and guessed from the filename."
    )

    parser.add_argument(
        "--field_table",
        default="closepack36_fields.fits",
        type=str,
        help="Table of field names for a given ASKAP survey. Default 'closepack36_fields.fits'"
    )

    parser.add_argument(
        "--beam_table",
        default="closepack36_beams.fits",
        type=str,
        help="Table of beam positions for a given ASKAP footprint. Default 'closepack36_beams.fits"
    )

    parser.add_argument(
        "--unwise_table_location",
        default="./",
        type=str,
        help="Directory of the unWISE tables. Default './'"
    )

    parser.add_argument(
        "-o", "--output-prefix", 
        type=str,
        default=None, 
        help="The prefix to base outputs onto"
    )
    # parser.add_argument(
    #     "--passes",
    #     type=int,
    #     default=1,
    #     help="Number of passes over the data should the iterative method attempt",
    # )
    # parser.add_argument(
    #     "--all-plots",
    #     action="store_true",
    #     help="If provided all plots will be produced. Otherwise a minimumal set will be",
    # )
    # parser.add_argument(
    #     "--report-statistics-throughout",
    #     action="store_true",
    #     help="Collect and report statistics each iteration",
    # )
    parser.add_argument(
        "--coord_keys",
        nargs=2,
        default=["ra", "dec"],
        type=str,
        help="Column names/keys in tables for (ra, dec). [Default ('ra', 'dec')]"
    )
    parser.add_argument(
        "--flux_keys",
        nargs=2,
        default=["int_flux", "peak_flux"],
        type=str,
        help="Column names/keys in tables for integrated and peak flux density. [Default ('int_flux', 'peak_flux')]"
    )
    parser.add_argument(
        "--rms_key",
        default="local_rms",
        type=str,
        help="Local rms column name/key in tables. Default 'local_rms'"
    )

    parser.add_argument(
        "--snr_min",
        default=5.,
        type=float,
        help="Minimum SNR of sources. Default 5"
    )

    parser.add_argument(
        "--iso_min",
        default=36.,
        type=float,
        help="Minimum separation between close neighbours in arcsec. Default 36"
    )

    return parser


def cli() -> None:
    parser = get_parser()

    args = parser.parse_args()

    unwise_shifts(
        catalogue_paths=args.paths,
        output_prefix=args.output_prefix,
        sbid=args.sbid,
        field_name=args.field_name,
        beam_table=args.beam_table,
        field_table=args.field_table,
        unwise_table_location=args.unwise_table_location,
        min_snr=args.snr_min,
        min_iso=args.iso_min,
        table_keys={
            "ra": args.coord_keys[0],
            "dec": args.coord_keys[1],
            "int_flux": args.flux_keys[0],
            "peak_flux": args.flux_keys[1],
            "local_rms": args.rms_key
        }
    )


if __name__ == "__main__":
    cli()