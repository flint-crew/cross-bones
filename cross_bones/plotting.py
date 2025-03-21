"""Some basic plotting procedures used throughout ``cross_bones``"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from cross_bones.catalogue import Catalogue, Catalogues
from cross_bones.matching import (
    OffsetGridSpace,
    calculate_matches,
    find_minimum_offset_space,
)


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


def plot_offset_grid_space(
    fname: str | Path,
    offset_grid_space: OffsetGridSpace,
    window: tuple[float, float, float, float, float],
) -> None:
    """Plot the offset surface of a match

    Args:
        fname (str | Path): The output file name to write
        offset_grid_space (OffsetGridSpace): The constructed offset surface
        window (tuple[float, float, float, float, float]): Details on the extent of the surface (min and max)
    """
    min_ra, min_dec, min_sep = find_minimum_offset_space(offset_grid_space)

    fig, ax = plt.subplots(1, 1)

    cim = ax.imshow(
        offset_grid_space.seps,
        extent=(window[0], window[1], window[2], window[3]),
        origin="lower",
    )

    ax.grid()
    ax.axhline(min_dec, ls="--", color="white")
    ax.axvline(min_ra, ls="--", color="white")

    ax.set(
        xlabel="Delta RA (arcsec)",
        ylabel="Delta Dec (arcsec)",
        title=f"Beam {offset_grid_space.beam} ({min_ra:.4f} {min_dec:.4f}) arcsec {offset_grid_space.n_sources} beam srcs",
    )
    _ = fig.colorbar(cim, label="Summed offsets (Degrees)")

    fig.tight_layout()
    plt.savefig(fname, dpi=150)
    plt.close()


def plot_offsets_in_field(
    offset_results: list[OffsetGridSpace],
    windows: list[tuple[float, float, float, float, float]],
    fname: str | Path,
) -> None:
    """Create a figure of each of the ``OffsetGridSpace`` where
    each surface is its own panel

    Args:
        offset_results (list[OffsetGridSpace]): The characterised surfaces to plot
        fname (str | Path): The output file name
    """
    from math import ceil

    num_columns = ceil(len(offset_results) ** 0.5)
    num_rows = ceil(len(offset_results) / num_columns)

    assert num_columns * num_rows >= len(offset_results), (
        f"The grid {num_columns=} {num_rows=} is not large enough for {len(offset_results)} results"
    )

    fig, axes = plt.subplots(num_columns, num_rows, figsize=(10, 10))

    for offset_result, ax, window in zip(offset_results, axes.flatten()[::-1], windows):
        minimum_point = find_minimum_offset_space(offset_space=offset_result)

        min_dec = minimum_point[1]
        min_ra = minimum_point[0]

        _ = ax.imshow(
            offset_result.seps,
            extent=(window[0], window[1], window[2], window[3]),
            origin="lower",
        )

        ax.set(title=f'({min_ra:.2f}, {min_dec:.2f})"')

        ax.grid()
        ax.axhline(min_dec, ls="--", color="white")
        ax.axvline(min_ra, ls="--", color="white")

    for ax in axes.flatten():
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_ticks([])
        ax.grid()

    fig.tight_layout()
    fig.savefig(fname, dpi=150)
    plt.close(fig)
