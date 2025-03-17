from __future__ import annotations

from argparse import ArgumentParser

from cross_bones.images import shift_image_collections


def get_parser() -> ArgumentParser:
    parser = ArgumentParser()

    parser.add_argument(
        "offset_file",
        type=str,
        help="File of offsets - must be in the same order as `images`.",
    )

    parser.add_argument(
        "images",
        type=str,
        nargs="*",
        help="List of images to shift. Must match order of offsets in `offset_file`.",
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
        "--suffix",
        type=str,
        default=".shift",
        help="Suffix to append to output shifted image name. Default '.shift'.",
    )

    return parser


def cli() -> None:
    parser = get_parser()
    args = parser.parse_args()

    shift_image_collections(
        images=args.images,
        offset_file=args.offset_file,
        outname_suffix=args.outname_suffix,
        dxs_key=args.dra_key,
        dys_key=args.ddec_key,
    )


if __name__ == "__main__":
    cli()
