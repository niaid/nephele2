#!/usr/bin/env python3

import shlex
import subprocess
import traceback
import argparse

from nephele2.pipelines.pipebase import PipeBase
from nephele2.pipelines import file_utils

DIRS_TO_KEEP = ["consensus", "filtered_snps", "filtered_indels", "reports"]
NEXTFLOW_PROJECT_NAME = "nephele2_covid19"
NEXTFLOW_WORK_DIR = "/tmp/nextflow_work_dir/"
SARS2_PIPELINE_PATH = "/usr/local/src/nephele2/pipelines/SARS-CoV2/"

class Covid19Pipeline(PipeBase):
    def __init__(self, args):
        super().__init__(args)

class Covid19PipelineError(Exception):
    pass

def main(args):
    try:
        exit_status = 0
        pipe = Covid19Pipeline(args)
        """
        Add the pipeline's dependencies to the logfile
        Hard code this here, we need to comeback and figure out a way to make Nextflow print out the version the way we want
        """
        depedencies_log = [
            "Starting SARS-CoV-2 pipeline...",
            "Dependencies:",
            "--------------------------------",
            "TRIMMOMATIC 0.39",
            "BWA 0.7.17",
            "PICARD 2.23.8",
            "GATK 4.1.9.0",
            "SAMTOOLS 1.9",
            "HTSLIB 1.11",
            "BCFTOOLS 1.11",
            "DEEPTOOLS 3.5.0",
            "PILON 1.23",
            "BEDTOOLS 2.27.1",
            "PYSAM 0.16.0.1",
            "PYPAIRIX 0.3.7",
            "SNPEFF 5.0",
            "--------------------------------",
        ]
        pipe.log.info("\n".join(depedencies_log))

        # Pick nextflow script based on data type
        if args.data_type == "COVID19_PE":
            script_name = "main_paired_end.nf"
        elif args.data_type == "COVID19_SE":
            script_name = "main_single_end.nf"
        else:
            raise Covid19PipelineError("Data type {} is not supported.".format(args.data_type))

        # Command to execute the Nextflow pipeline
        command = "nextflow run -name {} -work-dir {} {} --inputs_dir {} --out {} --map_file {}".format(
            NEXTFLOW_PROJECT_NAME,
            pipe.base_dir + NEXTFLOW_WORK_DIR,
            SARS2_PIPELINE_PATH + script_name,
            pipe.inputs_dir,
            pipe.outputs_dir,
            pipe.map_fp
        )

        pipe.log.info(command)

        command_args = shlex.split(command)
        result = subprocess.run(command_args, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        pipe.log.info(result.stdout.decode())
        if result.returncode == 1:
            # Nextflow somehow does not produce stderr, need the following to capture stderr
            stderr_capture_command = "nextflow log {} -f stderr".format(NEXTFLOW_PROJECT_NAME)
            stderr_capture_result = subprocess.run(shlex.split(stderr_capture_command), stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            stderr_msg = stderr_capture_result.stdout.decode()
            raise Covid19PipelineError(stderr_msg)
    except Covid19PipelineError as exc:
        exit_status = 1
        pipe.log.error("There is an error while executing the pipeline:")
        pipe.log.error(str(exc))
        pipe.log_to_db(job_id=args.job_id, stack=str(exc), msg=str(exc))
    except Exception as exc:
        exit_status = 1
        pipe.log.error("Pipeline Error:")
        pipe.log.error(traceback.format_exc())
        stack = traceback.format_exc()
        msg = str(exc) + '\nUnknown pipeline error.  Please see logfile.txt for more information.'
        pipe.log_to_db(job_id=args.job_id, stack=stack, msg=msg)
    finally:
        # Removing intermediate folders
        pipe.log.info('Removing intermediate dirs...')
        try:
            if args.get_merged_bam_files:
                DIRS_TO_KEEP.append("merged")
            file_utils.remove_intermediate_dirs(pipe.outputs_dir, DIRS_TO_KEEP)
        except Exception:
            pass

        if exit_status == 1:
            pipe.log.error('SARS-CoV-2 pipeline completed with errors.')
        else:
            pipe.log.info('SARS-CoV-2 pipeline complete.')
        pipe.log.info(exit_status)
        exit(exit_status)
    


if __name__ == '__main__':
    PARSER = argparse.ArgumentParser()
    PARSER.add_argument("--job_id", type=str,
                        help="job_id when running in Nephele. either "
                        "job_id or inputs_dir should be specified.")
    PARSER.add_argument("-i", "--inputs_dir", type=str,
                        help="input directory for running outside Nephele. "
                        "either job_id or inputs_dir should be specified.")
    PARSER.add_argument("-o", "--outputs_dir", type=str,
                        help="output directory for running outside Nephele. optional")
    PARSER.add_argument('--data_type', help='Ignore - placeholder.')
    PARSER.add_argument('--get_merged_bam_files', action='store_true')

    REQ_PARAMS = PARSER.add_argument_group('required arguments')

    REQ_PARAMS.add_argument('--map_file', type=argparse.FileType('r'),
                            required=True, help='required')

    ARGS = PARSER.parse_args()
    main(ARGS)
