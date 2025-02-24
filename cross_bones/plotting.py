from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from cross_bones.catalogue import Catalogue, Catalogues
from cross_bones.matching import calculate_matches


def plot_astrometric_offsets(
    catalogue_1: Catalogue, catalogue_2: Catalogue, ax: plt.axes | None = None
) -> plt.axes:
    """Plot the relative astrometry between two catalogues

    Args:
        catalogue_1 (Catalogue): The first catalogue
        catalogue_2 (Catalogue): The second catalogue
        ax (plt.axes | None, optional): Where to plot. If None it is created on a new Figure. Defaults to None.

    Returns:
        plt.axes: The axes where drawing occurred
    """

    src_matches = calculate_matches(catalogue_1, catalogue_2)
    mean_ra, mean_dec = src_matches.offset_mean
    std_ra, std_dec = src_matches.offset_std

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(3, 3))

    ax.scatter(src_matches.err_ra, src_matches.err_dec, color="red", s=1)
    ax.errorbar(mean_ra, mean_dec, xerr=std_ra, yerr=std_dec)
    ax.set(
        xlabel="Delta RA (arcsec)",
        ylabel="Delta Dec (arcsec)",
        title=f"Cata. 1: {catalogue_1.idx}, Cata. 2: {catalogue_2.idx} {len(src_matches.err_ra)} srcs",
        xlim=[-5, 5],
        ylim=[-5, 5],
    )
    ax.grid()
    ax.axvline(0, color="black", ls=":")
    ax.axhline(0, color="black", ls=":")

    return ax


def plot_beam_locations(
    catalogues: Catalogues,
    catalogue_1: Catalogue | None = None,
    catalogue_2: Catalogue | None = None,
    ax: plt.axes | None = None,
) -> plt.axes:
    """Plot the rough centre of the catalogues, and optionally present a pair of catalogues

    Args:
        catalogues (Catalogues): Collection of catalogues to plot
        catalogue_1 (Catalogue | None, optional): The first catalogue of a pair. Defaults to None.
        catalogue_2 (Catalogue | None, optional): The second catalogue of a pair. Defaults to None.
        ax (plt.axes | None, optional): An axes object to plot onto. If None, it will be created. Defaults to None.

    Returns:
        plt.axes: The axes object that was plotted to
    """
    ras = np.array([c.center.ra.deg for c in catalogues])
    decs = np.array([c.center.dec.deg for c in catalogues])
    fixed = np.array([c.fixed for c in catalogues])

    if ax is None:
        fig, ax = plt.subplots(1, 1)
    ax.scatter(ras[fixed], decs[fixed], color="red", marker="o")
    ax.scatter(ras[~fixed], decs[~fixed], color="green", marker="^")

    if catalogue_1 and catalogue_2:
        catalogue_1_pos = catalogue_1.center
        catalogue_2_pos = catalogue_2.center

        ax.scatter(
            catalogue_1_pos.ra.deg, catalogue_1_pos.dec.deg, color="black", marker="x"
        )
        ax.scatter(
            catalogue_2_pos.ra.deg, catalogue_2_pos.dec.deg, color="black", marker="x"
        )
        ax.plot(
            (catalogue_1_pos.ra.deg, catalogue_2_pos.ra.deg),
            (catalogue_1_pos.dec.deg, catalogue_2_pos.dec.deg),
            ls="--",
            c="black",
        )

    return ax
