#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" Not a real pipe - for demo purposes only. """

import argparse
import traceback
import sys
import sh
from nephele2.pipelines import pipeline_error
from nephele2.pipelines.pipebase import PipeBase


class DemoPipe(PipeBase):
    """Start pipe class definition."""

    def __init__(self, args):
        """Call into super on init, and set base args to args"""
        super().__init__(args, inputs_dir=args.inputs_dir)

    class Conf:
        """ Putting this inside the pipe class defn.  the names in these pipes
        can be a bit grim. Try to keep them in one place in here.  Caps denote
        fnames are constants.  """
        PLOT_DIR = 'loads_of_plots/'
        FILE_A = 'a_plot.png'
        FILE_B = 'x_plot.png'
        FILE_C = 'y_plot.png'
        ERROR_MSGS = {
            'UNIQUE_SEQS_OUT_FA':
            "There is an insufficient number of contigs/reads for analysis.",
            'SILVA_PCR_OUT_FNAME':
            "Contigs are too short; most likely due to low quality sequences.",
            'FILTER_UNIQUE_SEQS_OUT_FA':
            "There is an insufficient number of unique contigs.",
            'VSEARCH_COUNT_TABLE':
            "de novo clustering failed because of low sequence identity (< 97%).",
            'FILTER_PRECLUSTER_FA':
            "No chimeric sequences can be found.",
            'FILTER_GOOD_UNIQUE':
            "Sequences were completely filtered out because of low quality alignment."}
        # END OF Conf #

    @staticmethod
    def bad_thing(do_bad_thing):
        """simulates an error"""
        if do_bad_thing:
            raise pipeline_error.UnknownPipeError(
                msg='A bad thing happened in do_bad_thing() '
                'resulting from arg:' + str(do_bad_thing))

    @staticmethod
    def create_plots(dname):
        """mocking some outputs to move later..."""
        try:
            sh.mkdir('-p', dname)
            for fname in [dname+DemoPipe.Conf.FILE_A,
                          dname+DemoPipe.Conf.FILE_B,
                          dname+DemoPipe.Conf.FILE_C]:
                sh.touch(fname)
        except:
            raise

    @staticmethod
    def gen_demo_cmd(var_name):
        """Notice:
            the command is just a string and
            any relevant vars are passed to method (ie var_name)
            """
        return 'touch ' + var_name

    @staticmethod
    def final_step(is_a_success):
        """simulating successful operation"""
        if is_a_success:
            return
        raise pipeline_error.UnknownPipeError(
            msg='A bad thing happened in final_step() '
            'resulting from arg:' + str(is_a_success))


def main(args):
    """Called below by __main__
    """
    exit_status = 0

    try:
        pipe = DemoPipe(args)
        pipe.log.info('Starting pipe...')

        pipe.log.info('Attempting to make plots')
        pipe.create_plots(pipe.outputs_dir+DemoPipe.Conf.PLOT_DIR)
        for plot in [pipe.outputs_dir + DemoPipe.Conf.PLOT_DIR
                     + DemoPipe.Conf.FILE_A,
                     pipe.outputs_dir + DemoPipe.Conf.PLOT_DIR
                     + DemoPipe.Conf.FILE_B]:
            pipe.ensure_file_exists(plot)

        pipe.log.info('A bad thing might happen...')
        pipe.bad_thing(pipe.args.exit_1)

        try:
            cmnd = pipe.gen_demo_cmd('file.txt')
            pipe.log.info(cmnd)
            pipe.exec_cmnd(cmnd)
        except Exception:
            pipe.log.error('Core diversity failed.')
            exit(0)

        pipe.log.info('Running final step...')
        pipe.final_step(pipe.args.run_successfully)

    except pipeline_error.NecessaryFileNotFound as fnf:
        if fnf.fname is not None and fnf.fname in pipe.Conf.ERROR_MSGS.keys():
            fnf.msg = fnf.msg + pipe.Conf.ERROR_MSGS[fnf.fname]
            pipe.log.error(fnf.msg)
        else:
            pipe.log.error('Something strange happened??')
            pipe.log.error(fnf)
            sys.stderr.write(traceback.format_exc())
            pipe.log_to_db(job_id=pipe.job_id,
                           stack=traceback.format_exc(), msg=fnf.msg)
            exit_status = 1

    except Exception:
        sys.stderr.write(traceback.format_exc())
        pipe.log.error('Unknown Error arose:')
        pipe.log.error(traceback.format_exc())
        pipe.log_to_db(job_id=pipe.job_id,
                       stack=traceback.format_exc(), msg=traceback.format_exc())

        exit_status = 1
        exit(exit_status)


if __name__ == '__main__':
    PARSER = argparse.ArgumentParser()
    # Call from this set of args:
    MUTEX = PARSER.add_mutually_exclusive_group(required=True)
    MUTEX.add_argument("--job_id", type=str)
    MUTEX.add_argument("--inputs_dir", type=str)

    PARSER.add_argument(
        '--map_file', type=argparse.FileType('r'), required=True)
    VALID_DBS = ['SILVA', 'Greengenes']
    PARSER.add_argument('--database', choices=VALID_DBS, default='SILVA')
    PARSER.add_argument('--exit_1', action="store_true")
    PARSER.add_argument('--run_successfully', action="store_true")

    ARGS = PARSER.parse_args()
    main(ARGS)
