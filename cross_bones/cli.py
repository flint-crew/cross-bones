"""Cross-Bones CLI"""

from __future__ import annotations

import sys
from argparse import ArgumentParser

from cross_bones import (
    __version__,
    align_catalogues,
    align_catalogues_external,
    apply_shifts,
    shift_stats,
    unwise_align_catalogues,
)


def get_parser() -> ArgumentParser:
    """Create the CLI argument parser

    Returns:
        ArgumentParser: Constructed argument parser
    """
    parser = ArgumentParser(description="☠ cross-BONES ☠")
    subparsers = parser.add_subparsers(title="Available commands", dest="command")
    subparsers.required = True

    subparsers.add_parser(
        "internal",
        help="Align catalogues using internal cross-matching",
        parents=[align_catalogues.get_parser(parent_parser=True)],
    )
    subparsers.add_parser(
        "external",
        help="Align catalogues using external cross-matching",
        parents=[align_catalogues_external.get_parser(parent_parser=True)],
    )
    subparsers.add_parser(
        "apply",
        help="Apply offset shifts to images",
        parents=[apply_shifts.get_parser(parent_parser=True)],
    )
    subparsers.add_parser(
        "stats",
        help="Calculate statistics of offset shifts",
        parents=[shift_stats.get_parser(parent_parser=True)],
    )
    subparsers.add_parser(
        "unwise",
        help="Align catalogues using unWISE catalogues",
        parents=[unwise_align_catalogues.get_parser(parent_parser=True)],
    )
    parser.add_argument(
        "--version",
        "-V",
        action="version",
        version=f"%(prog)s {__version__}",
        help="Show version information",
    )
    return parser


def cli() -> None:
    """Main CLI function

    Parses command line arguments and executes the appropriate command
    based on the provided subcommand.
    """
    parser = get_parser()
    args = parser.parse_args()

    command_map = {
        "internal": align_catalogues.cli,
        "external": align_catalogues_external.cli,
        "apply": apply_shifts.cli,
        "stats": shift_stats.cli,
        "unwise": unwise_align_catalogues.cli,
    }

    # Execute the selected command
    if args.command in command_map:
        # Modify sys.argv to make it look like the command was called directly
        # This is because your cli() functions likely parse args themselves
        sys.argv = [sys.argv[0]] + sys.argv[2:]
        command_map[args.command]()
    else:
        parser.print_help()


if __name__ == "__main__":
    cli()
