from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord, search_around_sky
from numpy.typing import NDArray

from cross_bones.catalogue import Catalogue, make_sky_coords


@dataclass
class Match:
    """Components around matching Catalogue 1 to Catalogue 2"""

    sky_pos_1: SkyCoord
    """Sky positions from catalogue 1"""
    sky_pos_2: SkyCoord
    """Sky positions from catalogue 2"""
    matches: tuple[NDArray[Any], NDArray[Any], Any, Any]
    """The indices of the matches as returned by ``search_around_sky``"""
    match_1: SkyCoord
    """The sky-coordinate of a match in catalogue 1"""
    match_2: SkyCoord
    """The sky-coordinate of a match in catalogue 2"""
    n: int
    """Number of matches"""
    offset_mean: tuple[float, float]
    """Mean of the offset in arcseconds in the RA and Declination directions"""
    offset_std: tuple[float, float]
    """Std of the offset in arcseconds in the RA and Declination directions"""
    err_ra: NDArray[float]
    """Difference in RA coordinates between matches"""
    err_dec: NDArray[float]
    """Different in Dec coordinates between matches"""


@dataclass
class OffsetGridSpace:
    """Results of the offset surface fitting"""

    ra_offsets: float
    dec_offsets: float
    beam: int
    seps: NDArray[float]
    n_sources: int


def calculate_matches(
    catalogue_1: Catalogue, catalogue_2: Catalogue, sep_limit_arcsecond: float = 9
) -> Match:
    """Match a pair of catalogues to identify the sources in common.

    Args:
        catalogue_1 (Catalogue): The first loaded catalogue
        catalogue_2 (Catalogue): The second loaded catalogue
        sep_limit_arcsecond (float, optional): The separation limit for a match, in arcseconds. Defaults to None.

    Returns:
        Match: The result of the matching
    """

    sky_pos_1 = make_sky_coords(catalogue_1)
    sky_pos_2 = make_sky_coords(catalogue_2)

    matches = search_around_sky(
        sky_pos_1, sky_pos_2, seplimit=sep_limit_arcsecond * u.arcsec
    )
    match_1 = sky_pos_1[matches[0]]
    match_2 = sky_pos_2[matches[1]]

    # Extract the offsets of positions as angular offsets of the sphere
    deltas = match_1.spherical_offsets_to(match_2)
    err_ra = deltas[0].to(u.arcsec).value
    err_dec = deltas[1].to(u.arcsec).value

    mean_ra, mean_dec = np.median(err_ra), np.mean(err_dec)
    std_ra, std_dec = np.std(err_ra), np.std(err_dec)

    return Match(
        sky_pos_1=sky_pos_1,
        sky_pos_2=sky_pos_2,
        matches=matches,
        match_1=match_1,
        match_2=match_2,
        n=len(match_2),
        offset_mean=(float(mean_ra), float(mean_dec)),
        offset_std=(float(std_ra), float(std_dec)),
        err_ra=err_ra,
        err_dec=err_dec,
    )


def find_minimum_offset_space(
    offset_space: OffsetGridSpace,
) -> tuple[float, float, float]:
    """Search the input offset grid space to find the minimum position

    Args:
        offset_space (OffsetGridSpace): Results from the brute force grid search

    Returns:
        tuple[float, float, float]: The minimum minimum delta RA and delta Dec and the corresponding separation
    """
    minimum_sep = None
    minimum_ra = None
    minimum_dec = None
    for dec, ra, sep in zip(  # type: ignore[call-overload]
        offset_space.dec_offsets, offset_space.ra_offsets, offset_space.seps.flatten()
    ):
        if minimum_sep is None or minimum_sep > sep:
            minimum_sep = sep
            minimum_ra = ra
            minimum_dec = dec

    assert minimum_dec is not None
    assert minimum_ra is not None
    assert minimum_sep is not None

    return minimum_ra, minimum_dec, minimum_sep
