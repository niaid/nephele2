#!/usr/bin/env python3

"""
Pipeline for running wmgx and wmgx_vis bioBakery workflows.

"""

import argparse
import traceback
import os
import re
import csv
from collections import namedtuple
import sh
from nephele2.pipelines.pipebase import PipeBase
from nephele2.pipelines import pipeline_error


class BiobakeryPipe(PipeBase):
    def __init__(self, args):
        super().__init__(args)
        self._visoutputs_dir = os.path.join(self._outputs_dir, 'wmgx_vis')
#        os.makedirs(self._visoutputs_dir, exist_ok=True)

    @property
    def visoutputs_dir(self):
        return self._visoutputs_dir

    Sample = namedtuple('Sample', ['sample_id', 'fwd_file', 'rev_file'])

    @staticmethod
    def gen_samples(map_file, data_type):
        """Generate list of Sample namedtuples from mapping file.
        ``Sample = namedtuple('Sample', ['sample_id', 'fwd_file', 'rev_file'])``

        Parameters
        ----------
        map_file : file handle
            full path of mapping file from argparse
        data_type : string
            either WGS_SE or WGS_PE

        Returns
        -------
        samples : list
            list of Sample objects
        """
        samples = list()
        with open(map_file) as f_in:
            reader = csv.DictReader(f_in, delimiter='\t')
            for row in reader:
                if row['#SampleID'].strip():
                    if data_type == 'WGS_SE':
                        rev_file = ''
                    else:
                        rev_file = row['ReverseFastqFile']
                s = BiobakeryPipe.Sample(sample_id=row['#SampleID'],
                                         fwd_file=row['ForwardFastqFile'],
                                         rev_file=rev_file)
                samples.append(s)
        return samples

    @classmethod
    def rename_paired_end_files(cls, inputs_dir, outputs_dir, samples, file_ext):
        """Creates renamed_inputs subdirectory in outputs_dir.  Then, creates
        symbolic links in the subdirectory to the original sequence files.  The
        link names follow the format needed by biobakery_workflows - *SampleID.R1.fastq* and *SampleID.R2.fastq*

        Parameters
        ----------
        inputs_dir : path object
            full path of input directory which contains sequence files
        outputs_dir : path object
            full path of outputs directory where renamed links subdir will go
        samples : list
            list of namedtuple Sample objects
        file_ext : string
            file extension for sequence files

        Returns
        -------
        renamed_dir : path object
            path to renamed_inputs directory which contains renamed symbolic
            links to original sequence files
        """
        renamed_dir = os.path.join(outputs_dir, 'renamed_inputs') + '/'
        os.makedirs(renamed_dir, exist_ok=True)

        for s in samples:
            new_fwd = s.sample_id + '.R1.' + file_ext
            new_rev = s.sample_id + '.R2.' + file_ext
            cls.create_link(os.path.join(inputs_dir, s.fwd_file),
                            os.path.join(renamed_dir, new_fwd))
            cls.create_link(os.path.join(inputs_dir, s.rev_file),
                                os.path.join(renamed_dir, new_rev))

        return renamed_dir

    @staticmethod
    def gen_wmgx_cmd(inputs_dir, outputs_dir, strainphlan, threads, file_ext, local_jobs, keep):
        """Generate biobakery_workflows command for the wmgx pipeline.  The pipeline will run all the steps
        from here:
        https://bitbucket.org/biobakery/biobakery/wiki/biobakery_workflows#rst-header-metagenome-profiling

        Parameters
        ----------
        inputs_dir : path object
            full path of input directory which contains sequence files for processing
        outputs_dir : path object
            full path of output directory for the pipeline
        strainphlan : logical
            should strainphlan be run
        threads : string
            number of threads to run wmgx pipeline.  must be string to pass to biobakery_workflows
        file_ext : string
            file extension of input sequence files
        local_jobs : string
            number of wmgx tasks to run at once.  must be string to pass to biobakery_workflows
        keep : logical
            should we keep intermediate output

        Returns
        -------
        cmnd : string
            biobakery_workflows wmgx command
        """
        cmnd = 'biobakery_workflows wmgx'\
            + ' --input-extension ' + file_ext\
            + ' --threads ' + threads\
            + ' --input ' + inputs_dir\
            + ' --output ' + outputs_dir\
            + ' --skip-nothing'\
            + ' --local-jobs ' + local_jobs
        if not keep:
            cmnd += ' --remove-intermediate-output'
        if not strainphlan:
            cmnd += ' --bypass-strain-profiling'
        return cmnd

    def check_wmgx_outputs(self, outputs_dir):
        """Check if outputs of wmgx which are required by wmgx_vis exist.  If not, raise an exception.

        Parameters
        ----------
        outputs_dir : path object
            full path of output directory containing the output of the wmgx pipeline

        Raises
        ------
        `nephele2.pipelines.pipeline_error.NecessaryFileNotFound`

        """
        wmgx_files = [
            'kneaddata/merged/kneaddata_read_count_table.tsv',
            'metaphlan2/merged/metaphlan2_taxonomic_profiles.tsv', 'humann2/merged/pathabundance_relab.tsv',
            'humann2/counts/humann2_read_and_species_count_table.tsv',
            'humann2/counts/humann2_feature_counts.tsv'
        ]
        for f in wmgx_files:
            self.ensure_file_exists(os.path.join(outputs_dir, f))

    @staticmethod
    def gen_wmgx_vis_cmd(outputs_dir, visoutputs_dir, project_name, strainphlan, threads):
        """Generate biobakery_workflows command for the wmgx_vis pipeline.  The pipeline
        will run all of the steps from here:
        https://bitbucket.org/biobakery/biobakery/wiki/biobakery_workflows#rst-header-visualization
        and produce visualizations, including an html report.

        Parameters
        ----------
        outputs_dir : path object
            full path of output directory containing the output of the wmgx pipeline
        visoutputs_dir : path object
            full path of visoutputs subdirectory where the visualizations will go
        project_name : string
            project name to go at the top of the output html report
        strainphlan : logical
            was strainphlan run
        threads : string
            number of threads to run wmgx_vis pipeline.  must be string to pass to biobakery_workflows

        Returns
        -------
        cmnd : string
            biobakery_workflows wmgx_vis command
        """
        # sanitize project_name
        project_name = repr(project_name)
        cmnd = 'biobakery_workflows wmgx_vis'\
               + ' --input ' + outputs_dir\
               + ' --project-name ' + project_name\
               + ' --format html'\
               + ' --output ' + visoutputs_dir\
               + ' --local-jobs ' + threads\
               + ' --introduction-text \"The data was run through the standard workflow for whole metagenome shotgun sequencing'

        if not strainphlan:
            cmnd += '  with the exception of strain profiling (StrainPhlAn)'

        cmnd += '.  Details of the pipelines can be found in the <a href=https://bitbucket.org/biobakery/biobakery/wiki/biobakery_workflows#rst-header-metagenome-profiling>bioBakery Workflows Tutorial</a>.\"'
        return cmnd

    def check_wmgx_vis_outputs(self, visoutputs_dir):
        """Check if wmgx_vis report file exists.  If not, raise an exception.

        Parameters
        ----------
        visoutputs_dir : path object
            full path of directory containing the output of the wmgx_vis pipeline

        Raises
        ------
        `nephele2.pipelines.pipeline_error.NecessaryFileNotFound`
        """
        self.ensure_file_exists(os.path.join(
            visoutputs_dir, 'wmgx_report.html'))

    @staticmethod
    def cleanup_files(outputs_dir):
        """
        Clean up FASTQ files output by kneaddata and sam files output by
        metaphlan2 using os.remove.

        Parameters
        ----------
        outputs_dir : path object
            full path of output directory

        Returns
        -------
        cleanup_files_log : string
            log message with list of files which are removed
        """
        cleanup_files_log = ''
        # kneaddata files
        kneadmain = os.path.join(outputs_dir, 'kneaddata', 'main')
        if os.path.exists(kneadmain):
            for f in os.listdir(kneadmain):
                if re.search("\.fastq$", f):
                    rfile = os.path.join(kneadmain, f)
                    cleanup_files_log += 'Removing ' + str(rfile) + '\n'
                    os.remove(rfile)
        # metaphlan files
        metaphlanmain = os.path.join(outputs_dir, 'metaphlan2', 'main')
        if os.path.exists(metaphlanmain):
            for f in os.listdir(metaphlanmain):
                if re.search("\.sam$", f):
                    rfile = os.path.join(metaphlanmain, f)
                    cleanup_files_log += 'Removing ' + str(rfile) + '\n'
                    os.remove(rfile)
        return cleanup_files_log


def main(args):
    """Called below by __main__ """
    container_name = 'biobakery/nephele2'
    exit_status = 0

    try:
        pipe = BiobakeryPipe(args)
        # read in mapping file
        samples = pipe.gen_samples(args.map_file.name, args.data_type)

        # check for paired end data and merge files, if needed
        if pipe.args.data_type == 'WGS_PE':
            pipe.log.info('Renaming paired end files.')
            inputs_dir = pipe.rename_paired_end_files(pipe.inputs_dir, pipe.outputs_dir, samples,
                                                      pipe.args.file_ext)
        else:
            inputs_dir = pipe.inputs_dir

        pipe.log.info('Inputs directory: ' + inputs_dir)

        pipe.log.info('Running Whole Metagenome Shotgun Workflow (wmgx).')
        cmnd = pipe.gen_wmgx_cmd(inputs_dir, pipe.outputs_dir, pipe.args.strainphlan, pipe.args.threads,
                                 pipe.args.file_ext, pipe.args.local_jobs, pipe.args.keep)
        docker_cmnd = pipe.gen_docker_cmnd(cmnd, pipe.base_dir, container_name)
        pipe.log.info(docker_cmnd)
        pipe.exec_docker_cmnd(docker_cmnd)

        # check we have at least 3 samples to run wmgx_vis
        if len(samples) > 2:
            # make the visualization outputs directory
            pipe.log.info('Create wmgx_vis output directory: ' +
                          pipe.visoutputs_dir)
            sh.mkdir('-p', pipe.visoutputs_dir)
            sh.chmod('777', pipe.visoutputs_dir)
            sh.chown('www-data:www-data', pipe.visoutputs_dir)
            # os.makedirs(pipe.visoutputs_dir, mode=0o777, exist_ok=True)

            pipe.log.info(
                'Checking output files from wmgx workflow that are required by wmgx_vis workflow.')
            pipe.check_wmgx_outputs(pipe.outputs_dir)

            pipe.log.info(
                'Running Visualization for Whole Metagenome Shotgun Workflow (wmgx_vis).')
            cmnd = pipe.gen_wmgx_vis_cmd(pipe.outputs_dir, pipe.visoutputs_dir,
                                         pipe.args.project_name, pipe.args.strainphlan,
                                         pipe.args.threads)

            docker_cmnd = pipe.gen_docker_cmnd(
                cmnd, pipe.base_dir, container_name)
            pipe.log.info(docker_cmnd)
            pipe.exec_docker_cmnd(docker_cmnd)

            pipe.log.info('Checking output files from wmgx_vis pipeline.')
            pipe.check_wmgx_vis_outputs(pipe.visoutputs_dir)
        else:
            pipe.log.info('The Visualization for Whole Metagenome Shotgun Workflow '
                          '(wmgx_vis) will not be run, as at least 3 samples are needed.')

        pipe.log.info('bioBakery WGS pipeline done.')

    except pipeline_error.NecessaryFileNotFound as pp_err:
        pipe.log.error('Pipeline Error:')
        pipe.log.error(pp_err)
        pipe.log.error(
            'A step in the biobakery workflows may have failed. Check anadama.log files.')
        pipe.log.error('')
        exit_status = 1

    except Exception:
        pipe.log.error('Error:')
        pipe.log.error(traceback.format_exc())
        exit_status = 1

    finally:
        if not pipe.args.keep:
            try:
                pipe.log.info('Cleaning up intermediate files.')
                cleanup_files_log = pipe.cleanup_files(pipe.outputs_dir)
                pipe.log.info(cleanup_files_log)
            except Exception:
                pipe.log.error('Error:')
                pipe.log.error(traceback.format_exc())
                exit_status = 1
        exit(exit_status)


if __name__ == '__main__':
    PARSER = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # Call from this set of args:
    PARSER.add_argument('--job_id', type=str,
                        required=True, help="nephele job id")
    PARSER.add_argument('--map_file', action='store', type=argparse.FileType('r', encoding='latin-1'),
                        required=True, help="full path to mapping file")
    PARSER.add_argument('--data_type', type=str, choices=['WGS_SE', 'WGS_PE'], default='WGS_SE',
                        help="input sequence data type")
    PARSER.add_argument('--threads', type=str, help="number of processors per task (local job).  threads*local_jobs <= "
                        + str(os.cpu_count()) + ".  If none, will use total available processors/local_jobs.")
    PARSER.add_argument('--file_ext', choices=['fastq.gz', 'fastq', 'fq.gz', 'fq', 'fasta', 'fasta.gz'],
                        default='fastq', help="the input sequence file extension")
    PARSER.add_argument('--local_jobs', type=str, default=str(4),
                        help="number of bb tasks to run at a time.")
    PARSER.add_argument('--strainphlan', action='store_true',
                        help="user option. run strainphlan.")
    PARSER.add_argument('--keep', action='store_true',
                        help="debug option.  keep intermediate output.")
    PARSER.add_argument(
        '--project_name',
        help="user option. project name for visualization pipeline. if none, will use job_id.")
    PARSER.add_argument('-i', '--inputs_dir', type=str,
                        help="for running outside of Nephele")

    ARGS = PARSER.parse_args()
    ARGS.project_name = ARGS.project_name if ARGS.project_name else ARGS.job_id
    ARGS.threads = ARGS.threads if ARGS.threads else str(
        int(os.cpu_count()/int(ARGS.local_jobs)))
    main(ARGS)
