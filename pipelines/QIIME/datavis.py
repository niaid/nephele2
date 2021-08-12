#!/usr/bin/env python3

import argparse
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import traceback
from collections import defaultdict
import os
from rpy2.robjects.vectors import IntVector
from rpy2.robjects.packages import importr
import rpy2.rinterface
from biom.cli.table_summarizer import summarize_table
from nephele2 import config
from nephele2.pipelines.pipebase import PipeBase
from nephele2.pipelines import pipeline_error

class DadaWrapper(PipeBase):
    """
    **Written by**:
    Poorani Subramanian & Philip MacMenamin

    From within R this works:

    trycomputewrapper('/mnt/EFS/user_uploads/88bce3b496d9',
    '/mnt/EFS/user_uploads/88bce3b496d9',
    '/mnt/EFS/user_uploads/88bce3b496d9/mapfile.txt',
    logfilename='DADA_logfile.txt',
    refdb='/mnt/EFS/user_uploads/88bce3b496d9',
    refdb_species='/mnt/EFS/user_uploads/88bce3b496d9/silva_species_assignment_v128.fa')

    See :any:`dada2nephele.readme`

    This script is just a wrapper for the above R package.
    The above R package attempts to run function ``trycomputewrapper``, returns 0 on success and
    rpy2.rinterface.RRuntimeError on error.

    """
    def __init__(self, args):
        super().__init__(args)

    @staticmethod
    def run_datavis(kwargs):
        """
        run_datavis runs R function ``allgraphs`` from datavis16s R package.

        Parameters
        ----------
        kwargs : list
            See package documentation for information on ``allgraphs`` arguments:
            :any:`Reference_Manual_datavis16s`

        Returns
        -------
        exit_code : rpy2.robjects.vectors.IntVector
            `allgraphs` returns a rpy2.robjects.vectors.IntVector, which is 1 long.
            we test type (just to be safe) and return the first elt of this vector.
            else assume fail.

        Raises
        ------
        :any:`nephele2.pipelines.pipeline_error.RPipelineError`
        """
        try:
            datavis16s = importr('datavis16s')
            exit_code = datavis16s.trygraphwrapper(
                kwargs['datafile'],
                kwargs['out_dname'],
                kwargs['map_fname'],
                logfilename=kwargs['log_fname'],
                FUN='allgraphs',
                sampdepth=kwargs['sampdepth'])

            if isinstance(exit_code, IntVector):
                return exit_code[0]

            return 1

        except rpy2.rinterface.RRuntimeError as rpy2_err:
            raise pipeline_error.RPipelineError(rpy2_err, kwargs['job_id']) from None


def main(args):
    pipe = DadaWrapper(args)

    ##pipe_args are slightly different to args (coming in off the GUI)
    pipe_args = defaultdict()

    r_mod_name = 'DADA2/dada2nephele'
    pipe.log.info('Loading R module: {r_mod_name}.'.format(r_mod_name=r_mod_name))
    pipe.load_R_mod(config.PIPELINES_LOC_ON_WRKR + r_mod_name)

    pipe_args['job_id'] = args.job_id
    pipe_args['input_dname'] = "/home/admin/inputs"
    pipe_args['out_dname'] = "/home/admin/outputs"
    pipe_args['map_fname'] = args.map_file
    pipe_args['log_fname'] = "logfile.txt"
    pipe_args['datafile'] = args.biom
    pipe_args['sampdepth'] = 1000

    r_mod_name = 'datavis16s'
    pipe.load_R_mod(config.PIPELINES_LOC_ON_WRKR + r_mod_name)
    pipe.log.info('Running data visualization pipeline.')
    exit_status = pipe.run_datavis(pipe_args)

    pipe.log.info('DADA2 pipeline complete.')

if __name__ == '__main__':
    PARSER = argparse.ArgumentParser()
    PARSER.add_argument('--map_file', type=str, required=True, help="required")  # type=argparse.FileType('r')
    PARSER.add_argument('--job_id', type=str, help="if running inside Nephele")
    PARSER.add_argument('--sampling_depth', type=int, help="user option. sampling depth for downstream analysis. if none will calculate automatically.")
    PARSER.add_argument('--biom', type=str, help="biom file")
    ARGS = PARSER.parse_args()
    main(ARGS)
