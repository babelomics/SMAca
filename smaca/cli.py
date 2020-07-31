# coding: utf-8
"""
author: Daniel López
email: daniel.lopez.lopez@juntadeandalucia.es

author: Carlos Loucera
email: carlos.loucera@juntadeandalucia.es

SMA carrier Test program logic.
"""

import atexit
import cProfile
import io
import pstats

import click
import numpy as np

import smaca.constants as C
from smaca.sma import SmaCalculator


@click.command()
@click.option("--profile",
              is_flag=True,
              help="execution statistics (only for debug purposes)")
@click.option('--output',
              default="output.csv",
              type=click.Path(writable=True),
              help='output file')
@click.option('--ncpus', default=1, type=int, help='number of cores to use')
@click.option('--reference_version',
              default=C.REF_HG19,
              type=click.Choice([C.REF_HG19, C.REF_HG38]),
              help='reference genome version that was used for alignment')
@click.option('--reference_file',
              default=None,
              type=click.Path(exists=True),
              help='reference FASTA file. Highly recommended for CRAM files')
@click.option("--force_no_ref",
              default=False,
              is_flag=True,
              help="don't use reference genome file for CRAM alignment. (Not recommended)")
@click.argument("bam_list",
                type=click.Path(exists=True),
                nargs=-1,
                required=True)
def main(profile, output, bam_list, ncpus, reference_version, reference_file, force_no_ref):
    """
    Spinal Muscular Atrophy Carrier Analysis tool. Detect putative SMA carriers
    and estimate the absolute SMN1 copy-number in a population.

    """

    if not bam_list:
        ctx = click.get_current_context()
        ctx.get_help()
        ctx.exit()

    if np.any([str(b).lower().endswith("cram") for b in bam_list]):
        if not reference_file and not force_no_ref:
            print("ERROR: The reference FASTA file must be used for the analysis of CRAM files. "
                  "Please, use --reference_file to specify the full FASTA reference path (RECOMMENDED) or "
                  "use --force_no_ref to try it out online (not recommended)")
            ctx.exit()

    if profile:
        print("Profiling...")
        prf = cProfile.Profile()
        prf.enable()

        def exit():
            prf.disable()
            print("Profiling completed")
            ios = io.StringIO()
            pstats.Stats(prf,
                         stream=ios).sort_stats("cumulative").print_stats()
            print(ios.getvalue())

        atexit.register(exit)

    res = SmaCalculator(bam_list=bam_list, ref_version=reference_version, ref_file=reference_file, n_jobs=ncpus)
    res.write_stats(output)


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    main()
