from __future__ import annotations

import astropy.units as u
from astropy.coordinates import (
    SkyCoord,
)
from astropy.io import fits
from astropy.table import Table

from cross_bones.logging import logger


def apply_shift_to_image(
    image: str, dra: float, ddec: float, outname_suffix: str = ".shift", slice: int = 0
) -> None:
    """Apply RA, dec shifts to images by shifting reference pixels.


    Args:
        image (str): Input image name/path.
        dra (float): RA shift in degrees.
        ddec (float): DEC shift in degrees.
        outname_suffix (str): Suffix to append to output filename.
            Input image will be overwritten if this is set to an empty string ('').
            Default '.shift'.
        slice (int): HUD list slice. Default 0.


    """

    # TODO sanity check FITS/fits extensions...
    outname = image.replace(".fits", "") + outname_suffix + ".fits"
    logger.debug(f"Setting output FITS name: {outname}")

    with fits.open(image) as hdu:
        if "CRVAL1" in hdu[slice].header:
            keys = ("CRVAL1", "CRVAL2")
        elif "CD1_1" in hdu[slice].header:
            keys = ("CD1_1", "CD2_2")
        else:
            e = "Cannot determine reference coordinates keys in image header."
            raise RuntimeError(e)

        ref_coords = SkyCoord(
            ra=hdu[slice].header[keys[0]] * u.deg,
            dec=hdu[slice].header[keys[1]] * u.deg,
        )

        shifted_coords = ref_coords.spherical_offsets_by(
            dra * u.arcsec, ddec * u.arcsec
        )

        hdu[slice].header[keys[0]] = shifted_coords.ra.value
        hdu[slice].header[keys[1]] = shifted_coords.dec.value

        hdu.writeto(outname, overwrite=True)


def shift_image_collections(
    images: list[str],
    offset_file: str,
    outname_suffix: str = ".shift",
    dxs_key: str = "d_ra",
    dys_key: str = "d_dec",
) -> None:
    """
    Shift a collection of images (e.g. from ASKAP) based on an offset catalogue.

    Args:
        images (list[str]): List of images to shift. Must match order of offsets in `offset_file`.
        offset_file (str): File of offsets - must be in the same order as `images`.
        outname_suffix (str): Suffix to append to output filename.
            Input image will be overwritten if this is set to an empty string ('').
            Default '.shift'.
        dxs_key (str): RA shift key in `offsets`. Default 'd_ra'.
        dys_key (str): DEC shift key in `offsets`. Default 'd_dec'.

    """

    offset_table = Table.read(offset_file, format="csv")
    if len(offset_table) != len(images):
        e = f"Number of images {len(images)} must equal number of entries in the offset file {len(offset_table)}."
        raise RuntimeError(e)

    for i in range(len(offset_table)):
        dra = offset_table[i][dxs_key]
        ddec = offset_table[i][dys_key]

        logger.info(f'Offsetting {images[i]} with {dra}", {ddec}"')

        apply_shift_to_image(
            image=images[i], dra=dra, ddec=ddec, outname_suffix=outname_suffix
        )
