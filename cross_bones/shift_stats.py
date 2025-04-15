from __future__ import annotations

import warnings
from argparse import ArgumentParser
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.coordinates import Angle
from astropy.table import Table

from cross_bones.catalogue import (
    Catalogues,
    TableKeys,
    guess_sbid_and_field_racs,
    load_catalogues,
)
from cross_bones.logger import logger
from cross_bones.matching import calculate_matches

Paths = tuple[Path, ...]


def compare_askap_beams(
    comparison_type: str,
    catalogue_list: Paths,
    table_keys: TableKeys,
    dxs: list[float],
    dys: list[float],
    max_sep_from_beam: float | None = None,
    min_snr: float = 20.0,
    match_radius: float = 15.0,
    output_prefix: str | None = None,
) -> dict[str, object]:  # TODO Move output results to dataclass thing
    """Compare catalogues with shifts from ASKAP beams."""

    sbid, field_name = guess_sbid_and_field_racs(catalogue_list[0])

    all_dx = np.array([], dtype=float)
    all_dy = np.array([], dtype=float)

    catalogues: Catalogues = load_catalogues(
        catalogue_paths=catalogue_list,
        table_keys=table_keys,
        min_snr=min_snr,
    )

    mean_ra_offsets = np.full((len(catalogues) - 1,), np.nan)
    mean_dec_offsets = np.full((len(catalogues) - 1,), np.nan)
    median_ra_offsets = np.full((len(catalogues) - 1,), np.nan)
    median_dec_offsets = np.full((len(catalogues) - 1,), np.nan)
    std_ra_offsets = np.full((len(catalogues) - 1,), np.nan)
    std_dec_offsets = np.full((len(catalogues) - 1,), np.nan)
    n_offsets = np.full((len(catalogues) - 1,), 0.0)
    beams = []

    for b1 in range(len(catalogues) - 1):
        # Use b1 as the reference catalogue
        ref = catalogues[b1]
        beams.append(b1)

        assert ref.sky_coords is not None
        ### Note that this requires astropy>=6.1.0 (hence np 2.0)
        ### otherwise the offsets cannots be calculated against the skycoord list
        ### fudge: make a separate list
        ref_shift = ref.sky_coords.spherical_offsets_by(
            np.full((len(ref.table),), dxs[b1]) * u.arcsec,
            np.full((len(ref.table),), dys[b1]) * u.arcsec,
        )

        if len(ref.table) == 0:
            logger.warning(f"No sources in {catalogue_list[b1]}")
            continue

        ref_ra_offsets = np.array([])
        ref_dec_offsets = np.array([])

        # Compare reference with other neighbouring beams
        for beam in range(b1 + 1, len(catalogues)):
            comp_cat = catalogues[beam]
            # only compare nearest neighbours
            # approx centre
            if comp_cat.center.separation(ref.center).deg > 1.0:
                continue

            assert comp_cat.sky_coords is not None

            comp_shift = comp_cat.sky_coords.spherical_offsets_by(
                np.full((len(comp_cat.table),), dxs[beam]) * u.arcsec,
                np.full((len(comp_cat.table),), dys[beam]) * u.arcsec,
            )

            if max_sep_from_beam is not None:
                seps = comp_cat.center.separation(comp_shift).deg
                comp_shift = comp_shift[np.where(seps < max_sep_from_beam)[0]]

            if len(comp_cat.table) == 0:
                logger.warning(f"No sources in {catalogue_list[beam]}")
                continue

            match = calculate_matches(
                catalogue_1=ref_shift,
                catalogue_2=comp_shift,
                sep_limit_arcsecond=match_radius,
            )

            if match.n == 0:
                # no matches
                logger.warning(
                    f"No matches between {catalogue_list[b1]} and {catalogue_list[beam]}"
                )
                continue

            all_dx = np.append(all_dx, match.err_ra)
            all_dy = np.append(all_dy, match.err_dec)

            ref_ra_offsets = np.append(ref_ra_offsets, match.err_ra)
            ref_dec_offsets = np.append(ref_dec_offsets, match.err_dec)

        assert np.isfinite(ref_ra_offsets).all()
        assert np.isfinite(ref_dec_offsets).all()

        mean_ra_offsets[b1] = np.mean(ref_ra_offsets)
        mean_dec_offsets[b1] = np.mean(ref_dec_offsets)
        median_ra_offsets[b1] = np.median(ref_ra_offsets)
        median_dec_offsets[b1] = np.median(ref_dec_offsets)
        std_ra_offsets[b1] = np.std(ref_ra_offsets)
        std_dec_offsets[b1] = np.std(ref_dec_offsets)
        n_offsets[b1] = len(ref_ra_offsets)

    if output_prefix is not None:
        output_table = Table(
            [
                beams,
                n_offsets,
                mean_ra_offsets,
                median_ra_offsets,
                std_ra_offsets,
                mean_dec_offsets,
                median_dec_offsets,
                std_dec_offsets,
            ],
            names=[
                "beam",
                "n_matches",
                "mean_ra",
                "median_ra",
                "std_ra",
                "mean_dec",
                "median_dec",
                "std_dec",
            ],
        )

        output_table.write(
            output_prefix + "-" + comparison_type + "-beam-stats.csv",
            format="ascii.csv",
            overwrite=True,
        )

    mean_ra = np.nanmean(mean_ra_offsets)
    mean_dec = np.nanmean(mean_dec_offsets)
    median_ra = np.nanmedian(median_ra_offsets)
    median_dec = np.nanmedian(median_dec_offsets)
    std_ra = np.nanstd(mean_ra_offsets)
    std_dec = np.nanstd(mean_dec_offsets)

    loginfo = f"{comparison_type}: {field_name} {sbid} ({len(all_dx)}) {mean_ra:.3f} +/- {std_ra:.3f} {mean_dec:.3f} +/- {std_dec:.3f}"

    logger.info(loginfo)

    return {
        "field_name": field_name,
        "comparison": comparison_type,
        "sbid": sbid,
        "n_matches": len(all_dx),
        "mean_dra": mean_ra,
        "median_dra": median_ra,
        "std_dra": std_ra,
        "mean_ddec": mean_dec,
        "median_ddec": median_dec,
        "std_ddec": std_dec,
    }


def compare_original_to_fitted(
    catalogue_list: list[str],
    offset_file: str,
    table_keys: TableKeys,
    dxs_key: str = "d_ra",
    dys_key: str = "d_dec",
    max_sep_from_beam: float | None = None,
    min_snr: float = 20.0,
    output_prefix: str | None = None,
) -> tuple[dict[str, object], dict[str, object]]:
    """Compare original and fitted offsets by comparison to neighbouring beams.

    Args:
        catalogue_paths  (list):
    """

    warnings.filterwarnings("ignore")

    _catalogue_paths = [Path(path) for path in catalogue_list]
    _catalogue_paths.sort()
    catalogue_paths = tuple(_catalogue_paths)

    shift = Table.read(offset_file, format="csv")

    dxs = Angle(shift[dxs_key].value - shift[dxs_key].value, u.arcsec)
    dys = Angle(shift[dys_key].value - shift[dys_key].value, u.arcsec)

    r_original = compare_askap_beams(
        comparison_type="original",
        catalogue_list=catalogue_paths,
        dxs=dxs,
        dys=dys,
        max_sep_from_beam=max_sep_from_beam,
        table_keys=table_keys,
        min_snr=min_snr,
        output_prefix=output_prefix,
    )

    dxs = -Angle(shift[dxs_key].value, u.arcsec)
    dys = -Angle(shift[dys_key].value, u.arcsec)

    r_fitted = compare_askap_beams(
        comparison_type="fitted",
        catalogue_list=catalogue_paths,
        dxs=dxs,
        dys=dys,
        max_sep_from_beam=max_sep_from_beam,
        table_keys=table_keys,
        min_snr=min_snr,
        output_prefix=output_prefix,
    )

    if output_prefix is not None:
        out_table = Table(
            [
                [r_original["comparison"], r_fitted["comparison"]],
                [r_original["field_name"], r_fitted["field_name"]],
                [r_original["sbid"], r_fitted["sbid"]],
                [r_original["n_matches"], r_fitted["n_matches"]],
                [r_original["mean_dra"], r_fitted["mean_dra"]],
                [r_original["median_dra"], r_fitted["median_dra"]],
                [r_original["std_dra"], r_fitted["std_dra"]],
                [r_original["mean_ddec"], r_fitted["mean_ddec"]],
                [r_original["median_ddec"], r_fitted["median_ddec"]],
                [r_original["std_ddec"], r_fitted["std_ddec"]],
            ],
            names=[
                "comparison",
                "field_name",
                "sbid",
                "n_matches",
                "mean_ra",
                "median_ra",
                "std_ra",
                "mean_dec",
                "median_dec",
                "std_dec",
            ],
        )

        out_table.write(
            output_prefix + "-fit-stats.csv", format="ascii.csv", overwrite=True
        )

    return (r_original, r_fitted)


def get_parser(parent_parser: bool = False) -> ArgumentParser:
    parser = ArgumentParser(
        add_help=not parent_parser,
    )
    parser.add_argument(
        "offset_file",
        type=str,
        help="Output CSV file from cross_bones with per-catalogue offsets.",
    )

    parser.add_argument(
        "catalogues", nargs="*", type=str, help="List of catalogues to inspect."
    )

    parser.add_argument(
        "-r",
        "--max_sep_from_beam",
        type=float,
        default=None,
        help="Maximum distance from an image centre to consider sources. Default None",
    )

    parser.add_argument(
        "-o",
        "--output-prefix",
        type=str,
        default=None,
        help="The prefix to base outputs onto",
    )

    parser.add_argument(
        "--dra_key",
        default="d_ra",
        type=str,
        help="Column name in offset file for RA offsets. Default 'd_ra'.",
    )

    parser.add_argument(
        "--ddec_key",
        default="d_dec",
        type=str,
        help="Column name in offset file for RA offsets. Default 'd_dec'.",
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

    compare_original_to_fitted(
        catalogue_list=args.catalogues,
        offset_file=args.offset_file,
        max_sep_from_beam=args.max_sep_from_beam,
        dxs_key=args.dra_key,
        dys_key=args.ddec_key,
        table_keys=table_keys,
        min_snr=args.snr_min,
        output_prefix=args.output_prefix,
    )


if __name__ == "__main__":
    cli()
