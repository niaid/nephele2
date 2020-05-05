#!/usr/bin/env python

"""
Pipeline for performing quality control on fastq reads.

- See spec and usage here: :doc:`nephele2.pipelines.QC_reads_readme`
- You can't call this script directly due to qiime2 dependency; instead use qc_reads.sh_ , but the arguments and usage are otherwise the same.

.. _qc_reads.sh: https://github.com/niaid/nephele2/blob/master/pipelines/QC_reads/qc_reads.sh
"""

import argparse
import traceback
import os
import sys
import re
import subprocess
import sh

from qiime2.plugins import cutadapt
import qiime2.plugin.model.base

from nephele2.pipelines.pipebase import PipeBase
from nephele2.pipelines import pipeline_error
from nephele2.pipelines.QC_reads.trim import Q2Trim
from nephele2.pipelines.QC_reads.sample import Sample

# NEEDS TO RUN:
# source activate qiime2-2018.6
# PYTHONPATH=$PYTHONPATH:/usr/local/lib/python3.5/dist-packages/
# https://coderwall.com/p/ssuaxa/how-to-make-a-jar-file-linux-executable


class QCReads(PipeBase):
    """Start pipe class definition."""
    _limits_conf = os.path.dirname(os.path.realpath(__file__)) + '/limits.txt'
    _fastqc = sh.fastqc.bake('-t', str(os.cpu_count()), '-q', '-l', _limits_conf)
    _pe_trim = sh.trimmomatic.bake('PE', '-threads', str(os.cpu_count()), '-phred33')
    _se_trim = sh.trimmomatic.bake('SE', '-threads', str(os.cpu_count()), '-phred33')
    _flash2 = sh.flash2.bake('-D', '-z', '-t', str(os.cpu_count()))

    _yaml_conf = os.path.dirname(os.path.realpath(__file__)) + '/multiqc_config.yaml'
    _multiqc = sh.multiqc.bake('-c', _yaml_conf)

    def __init__(self, args):
        """Call into super on init, and set base args to args"""
        super().__init__(args, inputs_dir=args.inputs_dir)
        self.cutadapt_out_dir = self.outputs_dir + 'cutadapt_trimmed/'
        self.trim_out_dir = self.outputs_dir + 'qtrimmed_seqs/'
        self.flash2_d = self.outputs_dir + 'merged/'
        self.multiqc_in_dir = self.outputs_dir + 'multiqc_input/'

    @classmethod
    def fastqc(cls, fq_fp, out_fp):
        """executes the fastqc binary on fq_fp.

        Args:
          fq_fp: path to fastq file
          out_fp: output directory outputs/multiqc_input

        Raises:
          `nephele2.pipelines.pipeline_error.shError`
        """
        try:
            result = cls._fastqc('-o', out_fp, fq_fp)
        except sh.ErrorReturnCode as err:
            raise pipeline_error.shError(err.stdout, err.stderr)
        if result.exit_code != 0:
            raise pipeline_error.shError(result.stdout, result.stderr)
        return result

    @classmethod
    def pe_trim(cls, args):
        """executes the trimmomatic binary."""
        try:
            result = cls._pe_trim(args, _err_to_out=True)
        except sh.ErrorReturnCode as err:
            raise pipeline_error.shError(err.stdout, err.stderr)
        if result.exit_code != 0:
            raise pipeline_error.shError(result.stdout, result.stderr)
        return result

    @classmethod
    def se_trim(cls, args):
        """executes the trimmomatic binary."""
        try:
            result = cls._se_trim(args, _err_to_out=True)
        except sh.ErrorReturnCode as err:
            raise pipeline_error.shError(err.stdout, err.stderr)
        if result.exit_code != 0:
            raise pipeline_error.shError(result.stdout, result.stderr)
        return result

    @classmethod
    def flash2(cls, args):
        """executes the fastqc binary on fq_fp."""
        try:
            result = cls._flash2(args, _err_to_out=True)
        except sh.ErrorReturnCode as err:
            raise pipeline_error.shError(err.stdout, err.stderr)
        if result.exit_code != 0:
            raise pipeline_error.shError(result.stdout, result.stderr)
        return result

    @classmethod
    def multiqc(cls, inputs_d, out_dir, job_desc):
        """executes multiqc on its inputs directory."""
        try:
            if (job_desc and job_desc.strip()):
                result = cls._multiqc('--title', job_desc.strip(), '-o', out_dir, inputs_d)
            else:
                result = cls._multiqc('-o', out_dir, inputs_d)
        except sh.ErrorReturnCode as err:
            raise pipeline_error.shError(err.stdout, err.stderr)
        if result.exit_code != 0:
            raise pipeline_error.shError(result.stdout, result.stderr)
        return result

    @staticmethod
    def equal_pe_reads(fastqc_dir, sample):
        """
        Checks number of reads in fastqc html for paired end files.

        Args:
          fastqc_dir (str): directory contain fastqc output
          sample (sample.Sample): paired end sample to check

        Returns:
          (bool): ``True`` if R1 & R2 have equal numbers of reads; ``False`` otherwise.
        """
        pat = re.compile('<td>Total Sequences</td><td>(\d+)</td>')
        def find_total_seqs(fastqc_file, pat):
            with open (fastqc_file, 'rt') as myfile:
                contents = myfile.read()
            return(int(pat.search(contents).group(1)))
        fwd_fastqc=os.path.join(fastqc_dir, os.path.splitext(os.path.basename(sample.fwd_fp))[0]) + '_fastqc.html'
        rev_fastqc=os.path.join(fastqc_dir, os.path.splitext(os.path.basename(sample.rev_fp))[0]) + '_fastqc.html'
        return(find_total_seqs(fwd_fastqc, pat) == find_total_seqs(rev_fastqc, pat))

    @staticmethod
    def gen_cutadapt_args(cli_args):
        """
        Takes ALL cli args, because those args listed in str_arg : ['overlap', 'p-adapter-f'...
        can be strings or not exist. I look up the dict inside the argparse namespace,
        ask if the dict elt exists, and if it does I load it.
        I'm not sure if this is nasty or not (it's not great), ultimately I'd like
        a dict that looks like this:

        | cutadapt_args = {'indels':args.indels,
        | 'match_read_wildcards':args.match_read_wildcards,
        | 'match_adapter_wildcards':args.match_adapter_wildcards,
        | 'overlap':args.overlap,
        | 'error_rate':args.error_rate,
        | 'adapter_f':args.adapter_f,
        | 'adapter_r':args.adapter_r,
        | 'front_f':args.front_f,
        | 'front_r':args.front_r,
        | 'anywhere_f':args.anywhere_f,
        | 'anywhere_r':args.anywhere_r}, with the above logic
        | (eg if args.front_r is None, don't add it to the dict)

        """
        cutadapt_args = dict()
        cutadapt_args['cores'] = os.cpu_count()
        for bool_arg in ['indels', 'match_read_wildcards', 'match_adapter_wildcards']:
            cutadapt_args[bool_arg] = cli_args[bool_arg]
        for num_arg in ['overlap', 'error_rate']:
            if cli_args.get(num_arg) is not None:
                cutadapt_args[num_arg] = cli_args[num_arg]
        for str_arg in ['adapter_f', 'adapter_r', 'front_f', 'front_r', 'anywhere_f', 'anywhere_r']:
            if cli_args.get(str_arg):
                if cli_args.get('data_type') == 'PE':
                    cutadapt_args[str_arg] = [cli_args[str_arg]]  ## must be list
                else:
                    ## for single end we remove _f at the end of the arg for trim_single function
                    cutadapt_args[str_arg[:-2]] = [cli_args[str_arg]]

        return cutadapt_args

    @staticmethod
    def gen_flash_args(sample, out_dir, f2_min_overlap, f2_max_overlap, f2_min_overlap_outie,
                       f2_max_mismatch_density):
        """Generate arguments for flash2

        Args:
          sample (sample.Sample): for which to generate the flash2 command
          out_dir (str): pipe.flash2_d - outputs/merged
          f2_min_overlap (int): args.f2_min_overlap
          f2_max_overlap (int): args.f2_max_overlap
          f2_min_overlap_outie (int): args.f2_min_overlap_outie
          f2_max_mismatch_density (float): args.f2_max_mismatch_density

        Returns:
          flash_args (list): list of arguments to feed to flash2 - last 2 elements are files to be merged
        """
        if sample.trimmed_fwd_fp:
            fwd_fp = sample.trimmed_fwd_fp
            rev_fp = sample.trimmed_rev_fp
        elif sample.cutadpt_fwd_fp:
            fwd_fp = sample.cutadpt_fwd_fp
            rev_fp = sample.cutadpt_rev_fp
        else:
            fwd_fp = sample.fwd_fp
            rev_fp = sample.rev_fp

        flash_args = [
            '--min-overlap=' + str(f2_min_overlap), '--max-overlap=' + str(f2_max_overlap),
            '--min-overlap-outie=' + str(f2_min_overlap_outie),
            '--max-mismatch-density=' + str(f2_max_mismatch_density), '-o', sample.id + '_merged', '-d',
            out_dir
        ]
        flash_args.extend([fwd_fp, rev_fp])
        return flash_args

    @staticmethod
    def gen_trimmo_args(sample, out_d, lead_qual, trail_qual, window_size, req_qual, minlen, avg_qual):
        """generate arguments for trimmomatic

        Args:
          sample (sample.Sample): for which to generate the trimmomatic command
          out_d (str): path to trimmomatic output directory - outputs/qtrimmed_seqs
          input args (int): args.lead_qual, args.trail_qual, args.window_size, args.req_qual, args.minlen, args.avg_qual

        Returns:
          trimmo_args (list): list of arguments for trimmomatic.  first 1-2 (SE/PE) elements are the input fastq files.
        """
        lead = 'LEADING:' + str(lead_qual)
        trail = 'TRAILING:' + str(trail_qual)
        window = 'SLIDINGWINDOW:' + str(window_size) + ':' + str(req_qual)
        minlen = 'MINLEN:' + str(minlen)
        avg_qual = 'AVGQUAL:' + str(avg_qual)
        b_out = out_d + sample.id + '.fastq.gz'
        if sample.format == 'PE':
            if sample.cutadpt_fwd_fp:
                fwd_fp = sample.cutadpt_fwd_fp
            else:
                fwd_fp = sample.fwd_fp
            if sample.cutadpt_rev_fp:
                rev_fp = sample.cutadpt_rev_fp
            else:
                rev_fp = sample.rev_fp
            return [fwd_fp, rev_fp, '-baseout', b_out, lead, trail, window, minlen, avg_qual]
        else:
            if sample.cutadpt_fwd_fp:
                fwd_fp = sample.cutadpt_fwd_fp
            else:
                fwd_fp = sample.fwd_fp
            return [fwd_fp, b_out, lead, trail, window, minlen, avg_qual]


def main(args):
    try:
        exit_status = 0
        pipe = QCReads(args)
        pipe.log.info('Starting pipe...')
        sh.mkdir('-p', pipe.multiqc_in_dir)
        pipe.log.info('Trying to generate samples')
        samples = Sample.load_samples(pipe.map_fp, pipe.inputs_dir, args.data_type)

        link_name = pipe.multiqc_in_dir + os.path.basename(pipe.log.name)
        pipe.create_link(pipe.log.name, link_name)

        ##############
        ### FASTQC ###
        pipe.log.info('Starting Fastqc, outputs in %s.', pipe.multiqc_in_dir)
        for sample in samples:
            try:
                if sample.fwd_fp:
                    pipe.log.info('trying to run:\nfastqc -o %s %s',
                                  pipe.multiqc_in_dir, sample.fwd_fp)
                    pipe.fastqc(sample.fwd_fp, pipe.multiqc_in_dir)
                if sample.rev_fp:
                    pipe.log.info('trying to run:\nfastqc -o %s %s',
                                  pipe.multiqc_in_dir, sample.rev_fp)
                    pipe.fastqc(sample.rev_fp, pipe.multiqc_in_dir)
            except pipeline_error.shError as err:
                pipe.log.error(err.stdout.decode('utf-8'))
                pipe.log.error(err.stderr.decode('utf-8'))
                pipe.log.warning('%s will not be processed.\n', sample.id)
                sample.failed_step = 'FastQC'
                continue
            except:
                raise Exception

        ## check that some samples passed through FastQC ok ###
        if not Sample.filter_failed(samples):
            raise pipeline_error.PipelineError(job_id = args.job_id, msg='All samples failed FastQC.')

        pipe.log.info('Finished Fastqc, see %s for outputs.', pipe.multiqc_in_dir)
        ### END FASTQC ###
        ##################

        ##########################
        ### QIIME TRIM PRIMERS ###
        # see: https://docs.qiime2.org/2018.6/plugins/available/cutadapt/
        if args.run_cutadapt:
            sh.mkdir('-p', pipe.cutadapt_out_dir)
            cutadapt_log = pipe.multiqc_in_dir + 'cutadapt.log'
            stdoutlog = open(pipe.log.name, 'a')
            tee = subprocess.Popen(["tee", '-a', stdoutlog.name], stdin=subprocess.PIPE)
            os.dup2(tee.stdin.fileno(), sys.stdout.fileno())
            os.dup2(tee.stdin.fileno(), sys.stderr.fileno())

            pipe.log.info('Run cutadapt trim selected by user, '\
                          'see %s for outputs, see %s for cutadapt log.',
                          pipe.cutadapt_out_dir, cutadapt_log)

            ## generate qiime2 manifest for samples that passed FastQC
            pipe.log.info('Trying to generate QIIME2 manifest.')
            q2_mfest_fp = Q2Trim.gen_manifest(Sample.filter_failed(samples), pipe.cutadapt_out_dir)
            pipe.log.info('Done. See %s', q2_mfest_fp)

            pipe.log.info('Trying to generate input data.')
            q2_input_data = Q2Trim.gen_input_data(q2_mfest_fp, args.data_type)
            # Uncomment to HERE if you'd like to have the qza intermediate file.
            # q2_input_data_fp = pipe.cutadapt_out_dir + 'q2_input_demux.qza'
            # q2_input_data.save(q2_input_data_fp)
            # pipe.log.info('Done. See {fp}'.format(fp=q2_input_data_fp))
            # HERE

            cutadapt_args = pipe.gen_cutadapt_args(vars(args))
            if args.data_type == 'PE':
                # Note the IN MEMORY q2_input_data is used here, not the file (q2_input_data_fp)
                pipe.log.info('Inputs are Paired End, trying to run: qiime cutadapt trim-paired')
                pipe.log.info(cutadapt_args)
                trimmed_seqs = cutadapt.methods.trim_paired(demultiplexed_sequences=q2_input_data,
                                                            **cutadapt_args)
            else:
                pipe.log.info('Inputs are Single End, trying to run: qiime cutadapt trim-single')
                pipe.log.info(cutadapt_args)
                trimmed_seqs = cutadapt.methods.trim_single(demultiplexed_sequences=q2_input_data,
                                                            **cutadapt_args)

            tee.terminate()
            stdoutlog.close()

            pipe.log.info('Trying to write output to %s', pipe.cutadapt_out_dir)
            trimmed_seqs.trimmed_sequences.export_data(pipe.cutadapt_out_dir)

            for sample in Sample.filter_failed(samples):
                sample.load_sample_mfest_data_cutadpt(pipe.cutadapt_out_dir)
            pipe.log.info('Finished cutadapt.')
        ### END QIIME TRIM PRIMERS ###
        ##############################

        # ### Check that paired end files have equal numbers of reads ###
        # if args.data_type == 'PE' and (args.run_qual_trimming or args.run_flash2_merge) and not args.run_cutadapt:
        #     pipe.log.info("Checking that paired end files have equal numbers of reads.")
        #     for sample in Sample.filter_failed(samples):
        #         if not (pipe.equal_pe_reads(pipe.multiqc_in_dir, sample)):
        #             pipe.log.warning('%s will not be processed, as R1 and R2 have unequal numbers of reads.\n', sample.id)
        #             sample.failed_step = 'R1 and R2 have unequal numbers of reads.'

        #     if not Sample.filter_failed(samples):
        #         raise pipeline_error.PipelineError(job_id = args.job_id, msg='No samples have equal numbers of reads in R1 and R2. ' +
        #                                            'You can check the FastQC/multiQC reports to see the read counts.')



        ########################
        ### QUALITY TRIMMING ###
        if args.run_qual_trimming:
            pipe.log.info('Run Trimmomatic quality trimming selected.')
            pipe.log.info('Trying to create %s.', pipe.trim_out_dir)
            sh.mkdir('-p', pipe.trim_out_dir)
            trimo_log_fp = pipe.multiqc_in_dir + 'trimmo.log'

            for sample in Sample.filter_failed(samples):
                try:
                    trimo_args = pipe.gen_trimmo_args(sample, pipe.trim_out_dir, args.lead_qual,
                                                      args.trail_qual, args.window_size, args.req_qual,
                                                      args.minlen, args.avg_qual)
                    if args.data_type == 'PE':
                        pipe.log.info('Trying to run\ntrimmomatic PE -phred33 %s',
                                      ' '.join(trimo_args))
                        result = pipe.pe_trim(trimo_args)
                        sample.set_trimmed_pair(pipe.trim_out_dir + sample.id, '.fastq.gz')
                    else:
                        pipe.log.info('Trying to run\ntrimmomatic SE -phred33 %s',
                                      ' '.join(trimo_args))
                        result = pipe.se_trim(trimo_args)

                    pipe.log.info(result.stdout.decode('utf-8'))
                except pipeline_error.shError as err:
                    pipe.log.error(err.stdout.decode('utf-8'))
                    pipe.log.error(err.stderr.decode('utf-8'))
                    pipe.log.warning('%s will not be trimmed for quality.\n', sample.id)
                    sample.failed_step = 'Quality trimming'
                    continue
                except:
                    raise Exception

            pipe.log.info('Quality trimming complete, see %s.', trimo_log_fp)
        ### END QUALITY TRIMMING   ###
        ##############################

        ######################
        ### FLASH PE MERGE ###
        if args.run_flash2_merge:

            pipe.log.info('Run paired end merging selected.')
            if args.data_type != 'PE':
                pipe.log.error('Flash Merging is only Runnable on Paired End data!')
                raise ValueError('Flash Merging is only Runnable on Paired End data!')

            sh.mkdir('-p', pipe.flash2_d)
            for sample in Sample.filter_failed(samples):
                try:
                    flash_args = pipe.gen_flash_args(sample, pipe.flash2_d, args.f2_min_overlap,
                                                     args.f2_max_overlap, args.f2_min_overlap_outie,
                                                     args.f2_max_mismatch_density)
                    ## check if files exist and are nonempty
                    pipe.ensure_file_exists(flash_args[-1])
                    pipe.ensure_file_exists(flash_args[-2])

                    pipe.log.info('Trying to run:\n flash2 %s', " ".join(flash_args))
                    fl2 = pipe.flash2(flash_args)
                    pipe.log.info(fl2.stdout.decode('UTF-8'))
                except pipeline_error.NecessaryFileNotFound as fnf_err:
                    pipe.log.warning(fnf_err.msg)
                    pipe.log.warning('%s will not be merged.\n', sample.id)
                    sample.failed_step = 'Quality trimming'
                    continue
                except pipeline_error.shError as err:
                    pipe.log.error(err.stdout.decode('utf-8'))
                    pipe.log.error(err.stderr.decode('utf-8'))
                    pipe.log.warning('%s will not be merged.\n', sample.id)
                    sample.failed_step = 'Read merging'
                    continue
                except:
                    raise Exception

            pipe.log.info('Finished paired end merging.')
            pipe.log.info('Trying to link hist files for multiqc analysis.')
            for file in os.listdir(pipe.flash2_d):
                if file.endswith('.hist'):
                    mqc_link_name = pipe.multiqc_in_dir + os.path.basename(file)
                    pipe.create_link(pipe.flash2_d + file, mqc_link_name)
            pipe.log.info('Finished linking.')
        ### END FLASH PE MERGE ###
        ##########################

    ## Catch qiime2 exceptions ##
    except (qiime2.plugin.model.base.ValidationError, subprocess.CalledProcessError) as q2e:
        try:
            tee.terminate()
            stdoutlog.close()
        except:
            pass
        stack = traceback.format_exc()
        pipe.log.error(stack)
        pipe.log.error("QIIME2 Cutadapt Error:")
        pipe.log.error(q2e)
        msg = 'q2-cutadapt error: May occur because of incorrectly specified adapter/primer sequences or paired-end files with differing numbers of reads.\n'
        pipe.log.error(msg)
        msg += 'See logfile.txt for more information.'
        pipe.log_to_db(job_id=args.job_id, stack=stack, msg=msg)
        exit_status = 1

    ## known pipeline errors ##
    except pipeline_error.PipelineError as err:
        msg = 'PIPELINE SUMMARY\n' + err.msg + '\n'
        pipe.log.error("Pipeline Error:")
        stack = traceback.format_exc()
        pipe.log.error(stack)
        pipe.log.info(msg)
        exit_status = 1
        pipe.log.error('QC pipeline did not complete.')
        pipe.log.info(exit_status)
        pipe.log_to_db(job_id=args.job_id, stack=stack, msg=err.msg)

    ## Catch other super bad exceptions ##
    except Exception as exc:
        pipe.log.error("Pipeline Error:")
        stack = traceback.format_exc()
        pipe.log.error(stack)
        pipe.log.error(exc)
        pipe.log.error('QC pipeline did not complete.')
        exit_status = 1
        pipe.log.info(exit_status)
        msg = str(exc) + '\nUnknown pipeline error.  Please see logfile.txt for more information.'
        pipe.log.error(msg)
        pipe.log_to_db(job_id=args.job_id, stack=stack, msg=msg)

    try:
        ###############
        ### multiQC ###
        pipe.log.info('Trying to run multiqc on %s.', pipe.multiqc_in_dir)
        pipe.log.info('Linking logfile to multiqc inputs directory, %s.', pipe.multiqc_in_dir)
        link_name = pipe.multiqc_in_dir + os.path.basename(pipe.log.name)
        pipe.create_link(pipe.log.name, link_name)
        pipe.log.info('Finished linking.\n')
        mqc = pipe.multiqc(pipe.multiqc_in_dir, pipe.outputs_dir, pipe.job_description)
        pipe.log.info(mqc.stderr.decode('UTF-8'))
        pipe.log.info('Finished multiqc.\n')
        ### END multiQC ###
        ###################

       ### Print table of failed samples ###
        if Sample.filter_failed(samples, 'any'):
            bsamples = Sample.filter_failed(samples, 'any')
            # :-/
            # this ain't perl. At least let's not have magic numbers.

            scol = max(
                list(map(len, [s.id for s in bsamples])) + list(map(len, [s.failed_step
                                                                          for s in bsamples])) + [13]) + 2
            bstable = 'PIPELINE SUMMARY\n'
            bstable += 'The following samples failed a step in the pipeline:\n'
            bstable += 'Sample'.ljust(scol) + '\t' + 'Pipeline Step'.ljust(scol) + '\n'
            for bs in bsamples:
                bstable += bs.id.ljust(scol) + '\t' + bs.failed_step.ljust(scol) + '\n'
            pipe.log.info(bstable)



        ## Print additional error messages
        if exit_status == 1:
            pipe.log.error('q2-cutadapt failed.')
        else:
            ### Check proportion of successful samples ###
            if len(Sample.filter_failed(samples)) / len(samples) < 0.5:
                pipe.log.error('Greater than half the samples submitted failed.')
                pipe.log_to_db(job_id=args.job_id, msg='Greater than half the samples submitted failed a step '\
                           'in the pipeline.  Please see logfile.txt for more information.')
                exit_status = 1


    except:
        print(traceback.format_exc())
        exit_status = 1
        try:
            if hasattr(pipe, 'log'):
                pipe.log.error(traceback.format_exc())
            else:
                print(traceback.format_exc())

        except NameError as ne:
            pipe.log.error(ne)
            print(ne)
            exit_status = 1
        except:
            print(traceback.format_exc())
            pipe.log.error(traceback.format_exc())
            exit_status = 1

    finally:
        if exit_status == 1:
            pipe.log.error('QC pipeline completed with errors.')
        else:
            pipe.log.info('QC pipeline complete.')
        pipe.log.info(exit_status)
        exit(exit_status)

class FixDataType(argparse.Action):
    """If the arg --data_type is ^QC_ the QC_ is rm'd """

    def __call__(self, parser, namespace, values, option_string=None):
        if values == 'QC_PE':
            values = 'PE'
        elif values == 'QC_SE':
            values = 'SE'
        setattr(namespace, self.dest, values)




if __name__ == '__main__':

    PARSER = argparse.ArgumentParser()
    # Call from this set of args:
    REQ_PARAMS = PARSER.add_argument_group('required arguments')
    REQ_PARAMS.add_argument(
        "--job_id", type=str,
        help="job_id when running in Nephele.  either job_id or inputs_dir should be specified.")
    REQ_PARAMS.add_argument(
        "--inputs_dir", type=str,
        help="input directory for running outside Nephele.  either job_id or inputs_dir should be specified.")
    REQ_PARAMS.add_argument("--outputs_dir", type=str,
                            help="output directory for running outside Nephele. optional (I know ~PS)")
    REQ_PARAMS.add_argument('--map_file', type=argparse.FileType('r'), required=True, help='required')
    REQ_PARAMS.add_argument('--data_type', required=True, choices=['PE', 'SE', 'QC_PE',
                                                                   'QC_SE'], action=FixDataType,
                            help='required. QC_PE = PE and QC_SE = SE - either works.', type=str)

    # run_cutadapt params
    # Note, these args are grouped. If run_cutadapt is not set cutadapt args ARE IGNORED.
    CUTADAPT_PARAMS = PARSER.add_argument_group('cutadapt')
    CUTADAPT_PARAMS.add_argument('--run_cutadapt', action="store_true",
                                 help='Master Switch for running QIIME2 cudadapt. If this flag is set, '
                                 'also need to spec one or more of below params')
    CUTADAPT_PARAMS.add_argument('--error_rate', type=float, help='max allowed error rate')
    CUTADAPT_PARAMS.add_argument('--indels', action="store_true", help='allow indels')
    CUTADAPT_PARAMS.add_argument('--overlap', type=int, help='min overlap betw adapter and read')
    CUTADAPT_PARAMS.add_argument('--match_read_wildcards', action="store_true",
                                 help='interpret IUPAC wildcards in reads')
    CUTADAPT_PARAMS.add_argument('--match_adapter_wildcards', action="store_true",
                                 help='interpret IUPAC wildcards in adapters')

    CUTADAPT_PARAMS.add_argument('--adapter_f', type=str, help='trim 3\' end')
    CUTADAPT_PARAMS.add_argument('--adapter_r', type=str, help='trim 3\' end')
    CUTADAPT_PARAMS.add_argument('--front_f', type=str, help='trim 5\' end')
    CUTADAPT_PARAMS.add_argument('--front_r', type=str, help='trim 5\' end')
    CUTADAPT_PARAMS.add_argument('--anywhere_f', type=str, help='trim anywhere')
    CUTADAPT_PARAMS.add_argument('--anywhere_r', type=str, help='trim anywhere')
    # end run_cutadapt params

    # qual trim params
    QUAL_TRIM_PARAMS = PARSER.add_argument_group('trimmomatic')
    QUAL_TRIM_PARAMS.add_argument('--run_qual_trimming', action="store_true",
                                  help='Master Switch for running Trimmomatic. If this flag is used, '
                                  'also need to spec one or more of below params')
    QUAL_TRIM_PARAMS.add_argument('--window_size', type=int, default=4,
                                  help='specifies the number of bases to average across')
    QUAL_TRIM_PARAMS.add_argument('--req_qual', type=int, default=12,
                                  help='specifies the average quality required')
    QUAL_TRIM_PARAMS.add_argument('--lead_qual', type=int, default=3,
                                  help='Cut bases off the start of a read, if below a threshold quality')
    QUAL_TRIM_PARAMS.add_argument('--trail_qual', type=int, default=3,
                                  help='Cut bases off the end of a read, if below a threshold quality')
    QUAL_TRIM_PARAMS.add_argument('--minlen', type=int, default=30,
                                  help='Drop the read if it is below a specified length')
    QUAL_TRIM_PARAMS.add_argument('--avg_qual', type=int, default=0, metavar='AVGQUAL',
                                  help='Drop the read if the average quality is below the specified level')
    # end qual trim params

    # FLASH PE MERGE
    FLASH_PARAMS = PARSER.add_argument_group('flash2')
    FLASH_PARAMS.add_argument('--run_flash2_merge', action="store_true",
                              help='Can be run if data type is Paired End (PE)')
    FLASH_PARAMS.add_argument(
        '--f2_min_overlap', type=int, default=10,
        help='The minimum required overlap length between two reads to '
        'provide a confident overlap.')
    FLASH_PARAMS.add_argument(
        '--f2_max_overlap', type=int, default=315, help='Maximum overlap length expected in '
        'approximately 90%% of read pairs. Overlaps longer than the maximum'
        'overlap parameter are still considered as good overlaps, but '
        'the mismatch density (explained below) is calculated over the '
        'first max_overlap bases in the overlapped region rather than '
        'the entire overlap.')
    # FLASH_PARAMS.add_argument('--f2_allow_outies', action="store_true",
    #                     help='Also try combining read pairs in the "outie" orientation')
    FLASH_PARAMS.add_argument(
        '--f2_min_overlap_outie', type=int, default=35,
        help='The minimum required overlap length between two reads to '
        'provide a confident overlap in an outie scenario.')
    FLASH_PARAMS.add_argument(
        '--f2_max_mismatch_density', type=float, default=0.25,
        help='Maximum allowed ratio between the number of mismatched '
        'base pairs and the overlap length. ')
    # END FLASH PE MERGE

    ARGS = PARSER.parse_args()
    main(ARGS)
