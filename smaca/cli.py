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
from smaca.sma import SmaCalculator


@click.command()
@click.option("--profile", is_flag=True)
@click.option('--output',
              default="output.csv",
              type=click.Path(writable=True),
              help='output file')
@click.option('--ncpus', default=1, type=int, help='number of cores to use')
@click.argument("bam_list",
                type=click.Path(exists=True),
                nargs=-1,
                required=True)
def main(profile, output, bam_list, ncpus):
    """
    Predict proportion of SMN1:SMN2 for a set of BAM files.

    bam_list: input BAM files list
    """
    if not bam_list:
        ctx = click.get_current_context()
        ctx.get_help()
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

    res = SmaCalculator(bam_list, n_jobs=ncpus)
    res.write_stats(output)


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    main()
