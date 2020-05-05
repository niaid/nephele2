#!/usr/bin/env python3
import traceback
import argparse
import os
import gc
import re
import warnings
from collections import defaultdict
import rpy2.rinterface
from rpy2.robjects.vectors import IntVector
from rpy2.robjects.packages import importr
import rpy2.robjects as robjects
from biom.cli.table_summarizer import summarize_table
from nephele2 import config
from nephele2.pipelines.pipebase import PipeBase
from nephele2.pipelines import pipeline_error
from nephele2.pipelines import r_pipeline_error
from nephele2.pipelines.DADA2 import dada2_config
warnings.simplefilter(action='ignore', category=FutureWarning)


class DadaWrapper(PipeBase):
    """

    From within R this works:

    .. code-block:: R

      trycomputewrapper('/mnt/EFS/user_uploads/88bce3b496d9',
      '/mnt/EFS/user_uploads/88bce3b496d9',
      '/mnt/EFS/user_uploads/88bce3b496d9/mapfile.txt',
      logfilename='DADA_logfile.txt',
      refdb='/mnt/EFS/user_uploads/88bce3b496d9',
      refdb_species='/mnt/EFS/user_uploads/88bce3b496d9/silva_species_assignment_v128.fa')

    See :ref:`trycomputewrapper <Reference_Manual_dada2nephele:dada2compute>`.

    This script is just a wrapper for the `dada2nephele R package <pipelines.DADA2.dada2nephele.readme>`.
    The above R package attempts to run function ``trycomputewrapper``, returns 0 on success and
    ``rpy2.rinterface.RRuntimeError`` on error.

    """

    def __init__(self, args):
        super().__init__(args)

    @staticmethod
    def convert_none(pyobj):
        """Convert python ``None`` to R ``NULL``

        Parameters
        ----------
        pyobj : any python object

        Returns
        -------
        pyobj : the same python object if not ``None`` or an ``rpy2.robjects.rinterface.NULL`` object.
        """
        if pyobj is None:
            return robjects.rinterface.NULL
        return pyobj

    @staticmethod
    def run_dada2(kwargs):
        """
        run_dada2 runs R function ``trycomputewrapper`` from dada2nephele R package.

        Parameters
        ----------
        kwargs : list
            See package documentation for information on ``trycomputewrapper`` arguments:
            `Reference_Manual_dada2nephele:dada2compute`

        Returns
        -------
        exit_code : rpy2.robjects.vectors.IntVector
            ``trycomputewrapper`` returns a ``rpy2.robjects.vectors.IntVector``, which is 1 long.
            we test type (just to be safe) and return the first elt of this vector.
            else assume fail.

        Raises
        ------
        `nephele2.pipelines.r_pipeline_error.RPipelineError`
        """
        try:
            kwargs = {k: DadaWrapper.convert_none(v) for k, v in kwargs.items()}
            dada2nephele = importr('dada2nephele')
            exit_code = dada2nephele.trycomputewrapper(
                kwargs['input_dname'], kwargs['out_dname'], kwargs['map_fname'],
                logfilename=kwargs['log_fname'], refdb=kwargs['ref_db'],
                refdb_species=kwargs['ref_db_species'], chimera=kwargs['chimera'],
                trimLeft=kwargs['trimLeft'], trimOverhang=kwargs['trimOverhang'], nthread=os.cpu_count() - 1,
                data_type=kwargs['data_type'], maxEE=kwargs['maxEE'], truncQ=kwargs['truncQ'],
                truncLen=kwargs['truncLen'], maxMismatch=kwargs['maxMismatch'],
                justConcatenate=kwargs['justConcatenate'], taxmethod=kwargs['taxmethod'],
                band_size=kwargs['band_size'], homopolymer_gap_penalty=kwargs['homopolymer_gap_penalty'])

            if isinstance(exit_code, IntVector):
                return exit_code[0]

            return 1

        except rpy2.rinterface.RRuntimeError as rpy2_err:
            raise r_pipeline_error.RPipelineError(rpy2_err, kwargs['job_id']) from None

    @staticmethod
    def run_datavis(kwargs):
        """
        run_datavis runs R function ``allgraphs`` from datavis16s R package.

        Parameters
        ----------
        kwargs : list
            See `package documentation <pipelines.datavis16s.readme:R package manual>` for information on ``allgraphs`` arguments.

        Returns
        -------
        exit_code : rpy2.robjects.vectors.IntVector
            ``allgraphs`` returns a rpy2.robjects.vectors.IntVector, which is 1 long.
            we test type (just to be safe) and return the first elt of this vector.
            else assume fail.

        Raises
        ------
        `nephele2.pipelines.r_pipeline_error.RPipelineError`
        """
        try:
            datavis16s = importr('datavis16s')
            exit_code = datavis16s.trygraphwrapper(kwargs['datafile'], kwargs['out_dname'],
                                                   kwargs['map_fname'], logfilename=kwargs['log_fname'],
                                                   FUN='allgraphs', tsvfile=True,
                                                   sampdepth=kwargs['sampdepth'])

            if isinstance(exit_code, IntVector):
                return exit_code[0]

            return 1

        except rpy2.rinterface.RRuntimeError as rpy2_err:
            raise r_pipeline_error.RPipelineError(rpy2_err, kwargs['job_id']) from None

    class Conf:
        """ Contains constants and error messages.
        When `RPipelineError` is thrown,
        the error message is checked for a known error in ERROR_REGEX, and if it exists
        the user message in ERROR_MSGS is set as msg in logging to the database.
        """
        #: maps from database option to database filename
        REFDB = {
            'sv99': 'dada2_silva_v132/silva_nr_v132_train_set.fa',
            'homd': 'dada2_homd/HOMD_16S_rRNA_RefSeq_V15.1.train_set.fa',
            'idtaxa': 'dada2_silva_v132/SILVA_SSU_r132_March2018.RData',
            'gg97': 'dada2_greengenes/gg_13_8_train_set_97.fa.gz'
        }
        #: maps from database option to species database filename
        REFDB_SPECIES = {
            'sv99': 'dada2_silva_v132/silva_species_assignment_v132.fa',
            'homd': 'dada2_homd/HOMD_16S_rRNA_RefSeq_V15.1.species_assignment.fa',
            'idtaxa': None,
            'gg97': None
        }
        #: maps from regex to message for user
        ERROR_MSGS = {
            '(Mismatched forward and reverse sequence files)':
            'you may want to check that the forward and reverse files have the same number of reads, '
            'and the reads are properly paired in matching order.  Sometimes, this indicates files '
            'were improperly demultiplexed prior to job submission.',

            '(Error in names\(seqs\) <- paste0\(\"seq\", 1:length\(seqs\)\) : attempt to set an '
            'attribute on NULL)':
            'indicates no sequence variants were produced after denoising and merging reads (for PE).  '
            'You may want to examine the dataset quality and modify your filterAndTrim or mergePairs '
            '(for PE) parameters.',

            '(Input must be a valid sequence table).+removeBimeraDenovo':
            'indicates sequence table is empty because no sequence variants were produced after denoising '
            'and merging reads (for PE).  You may want to examine the dataset quality and modify your '
            'filterAndTrim or mergePairs (for PE) parameters.'
        }
        #: contains compiled regex statements for efficiency
        ERROR_REGEX = []

        for raw_errmsg in ERROR_MSGS:
            ERROR_REGEX.append(re.compile(raw_errmsg))


def main(args):
    """main pipeline
    """
    exit_status = 0

    try:

        ## initialize pipeline
        pipe = DadaWrapper(args)

        ##pipe_args are slightly different to args (coming in off the GUI)
        pipe_args = defaultdict()

        pipe_args['job_id'] = args.job_id
        pipe_args['input_dname'] = pipe.inputs_dir
        pipe_args['out_dname'] = pipe.outputs_dir
        pipe_args['map_fname'] = args.map_file.name
        pipe_args['log_fname'] = str(pipe.log.name)

        ## DADA2 computational parameters
        pipe_args['maxEE'] = args.maxee
        pipe_args['truncQ'] = args.truncq
        pipe_args['maxMismatch'] = args.maxmismatch
        pipe_args['justConcatenate'] = args.just_concatenate
        pipe_args['chimera'] = args.chimera
        pipe_args['trimOverhang'] = args.trim_overhang
        pipe_args['data_type'] = args.data_type

        if args.data_type == 'PE':
            pipe_args['trimLeft'] = [args.trimleft_fwd, args.trimleft_rev]
            pipe_args['truncLen'] = [args.trunclen_fwd, args.trunclen_rev]
        else:
            pipe_args['trimLeft'] = args.trimleft_fwd
            pipe_args['truncLen'] = args.trunclen_fwd

        #Ion Torrent https://benjjneb.github.io/dada2/faq.html#can-i-use-dada2-with-my-454-or-ion-torrent-data
        if args.ion_torrent:
            pipe_args['band_size'] = 32
            pipe_args['homopolymer_gap_penalty'] = -1
        else:
            pipe_args['band_size'] = None
            pipe_args['homopolymer_gap_penalty'] = None

        # Databases
        pipe_args['taxmethod'] = args.taxmethod
        if args.taxmethod == 'idtaxa':
            args.ref_db = 'idtaxa'

        pipe_args['ref_db'] = pipe.db_dir + pipe.Conf.REFDB[args.ref_db]
        if pipe.Conf.REFDB_SPECIES[args.ref_db] is not None:
            pipe_args['ref_db_species'] = pipe.db_dir + pipe.Conf.REFDB_SPECIES[args.ref_db]
        else:
            pipe_args['ref_db_species'] = None

        ## graphs parameter
        pipe_args['datafile'] = pipe.outputs_dir + dada2_config.OTUTABLE

        ## Load DADA2 R package
        if pipe.args.job_id is not None:
            # using this to see if we are running inside Nephele env
            r_mod_name = 'DADA2/dada2nephele'
            pipe.log.info('Loading R module: %s.', r_mod_name)
            pipe.load_R_mod(config.PIPELINES_LOC_ON_WRKR + r_mod_name)

        pipe.log.info('Running DADA2.')
        exit_status = pipe.run_dada2(pipe_args)

        pipe.log.info('Summarizing biom file to %s%s.', pipe.outputs_dir, dada2_config.BIOMSUMMARY)
        pipe.ensure_file_exists(pipe.outputs_dir + dada2_config.BIOMFILE)
        summarize_table.callback(pipe.outputs_dir + dada2_config.BIOMFILE,
                                 pipe.outputs_dir + dada2_config.BIOMSUMMARY, False, False)

        ## do garbage collection in python
        ## https://stackoverflow.com/questions/5199334/clearing-memory-used-by-rpy2
        gc.collect()

        pipe.log.info('Checking output file from dada2 pipeline required by data visualization pipeline.')
        pipe.ensure_file_exists(pipe.outputs_dir + dada2_config.OTUTABLE)

        ## Run visualization pipeline
        ## use logic in pipeline to provide partial output based on sampling depth.
        g_depth = pipe.get_depth(pipe.outputs_dir + dada2_config.BIOMFILE, args.sampling_depth)

        if args.sampling_depth is not None:
            pipe_args['sampdepth'] = args.sampling_depth
        elif g_depth > 0:
            pipe_args['sampdepth'] = g_depth
        else:
            pipe_args['sampdepth'] = 10000

        if pipe.args.job_id is not None:
            r_mod_name = 'datavis16s'
            pipe.log.info('Loading R module: {r_mod_name}.'.format(r_mod_name=r_mod_name))
            pipe.load_R_mod(config.PIPELINES_LOC_ON_WRKR + r_mod_name)
        pipe.log.info('Running data visualization pipeline.')
        exit_status = pipe.run_datavis(pipe_args)

        pipe.log.info('DADA2 pipeline complete.')

    except r_pipeline_error.RPipelineError as rerr:
        pipe.log.error('R Pipeline Error:')
        pipe.log.error(rerr.msg)

        ## search for error_regex patterns in rerr.msg
        for idx, error_regex in enumerate(pipe.Conf.ERROR_REGEX):
            match = error_regex.search(rerr.stack)
            if match:
                rerr.msg = match.group(1) + ' - ' + list(pipe.Conf.ERROR_MSGS.values())[idx]
                # non_smart_err = False
                break

        rerr.msg = rerr.msg + ' Please refer to logfile.txt for more information.'
        pipe.log_to_db(job_id=rerr.job_id, stack=rerr.stack, msg=rerr.msg)
        exit_status = 1

    except pipeline_error.NecessaryFileNotFound as file_err:
        pipe.log.error('File Not Found Error:')
        pipe.log.error(file_err.msg)
        stack = traceback.format_exc()
        pipe.log_to_db(job_id=args.job_id, stack=stack, msg=file_err.msg)
        exit_status = 1

    except pipeline_error.RefDBError as db_err:
        pipe.log.error('Reference DB Error:')
        stack = traceback.format_exc()
        pipe.log.error(db_err.msg)
        pipe.log_to_db(job_id=args.job_id, stack=stack, msg=db_err.msg)
        exit_status = 1

    except pipeline_error.PipelineError as pp_err:
        pipe.log.error('Pipeline Error:')
        stack = traceback.format_exc()
        pipe.log.error(stack)
        pipe.log.error(pp_err.msg)
        pipe.log_to_db(job_id=args.job_id, stack=stack, msg=pp_err.msg)
        exit_status = 1

    except Exception:
        pipe.log.error(traceback.format_exc())
        exit_status = 1
    finally:
        pipe.log.info(exit_status)
        exit(exit_status)


if __name__ == '__main__':

    PARSER = argparse.ArgumentParser()
    PARSER.add_argument('--job_id', type=str, help="if running inside Nephele")
    PARSER.add_argument('-i', '--inputs_dir', type=str, help="for running outside of Nephele")
    PARSER.add_argument('-o', '--outputs_dir', type=str, help="for running outside of Nephele")
    PARSER.add_argument('--map_file', type=argparse.FileType('r'), required=True, help="required")  # type=str
    PARSER.add_argument('--data_type', type=str, choices=['SE', 'PE'], default='PE',
                        help="single end or paired end.")
    USER_OPT = PARSER.add_argument_group('user options')
    USER_OPT.add_argument('--ion_torrent', action='store_true', help="dada2::dada use recommended "
                        "parameters for ion torrent data")
    USER_OPT.add_argument('--trimleft_fwd', type=int,
                        help="dada2::filterAndTrim trimLeft bp of fwd read.",
                        default=int(dada2_config.TRIMLEFT))
    USER_OPT.add_argument('--trimleft_rev', type=int,
                        help="dada2::filterAndTrim trimLeft bp of rev read.",
                        default=int(dada2_config.TRIMLEFT))
    USER_OPT.add_argument('--maxee', type=int, help="dada2::filterAndTrim discard reads "
                        "with EE higher than this value.", default=int(dada2_config.MAXEE))
    USER_OPT.add_argument('--trunclen_fwd', type=int,
                        help="dada2::filterAndTrim truncate fwd reads at this length.",
                        default=int(dada2_config.TRUNCLEN))
    USER_OPT.add_argument('--trunclen_rev', type=int,
                        help="dada2::filterAndTrim truncate rev reads at this length.",
                        default=int(dada2_config.TRUNCLEN))
    USER_OPT.add_argument('--truncq', type=int, help="dada2::filterAndTrim truncate "
                        "reads at first bp with qual <= this value.", default=int(dada2_config.TRUNCQ))
    USER_OPT.add_argument('--just_concatenate', action="store_true",
                        help="dada2::mergePairs concatenates instead of merging reads")
    USER_OPT.add_argument('--maxmismatch', type=int,
                        help="dada2::mergePairs max mismatches allowed.",
                        default=int(dada2_config.MAXMISMATCH))
    USER_OPT.add_argument('--trim_overhang', action="store_true",
                        help="dada2::mergePairs trims overhanging seq.")
    USER_OPT.add_argument('--chimera', action="store_true", help="run dada2::removeBimeraDenovo.")
    USER_OPT.add_argument('--ref_db', type=str, default="sv99", choices=["sv99", "homd", "gg97"],
                        help="reference database")
    USER_OPT.add_argument('--taxmethod', type=str, default=dada2_config.TAXMETHOD, choices=["rdp", "idtaxa"],
                        help='taxonomic assignment method')
    USER_OPT.add_argument('--sampling_depth', type=int,
                        help="sampling depth for downstream analysis. optional")

    ARGS = PARSER.parse_args()
    main(ARGS)
