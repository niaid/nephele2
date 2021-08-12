#!/usr/bin/env python3

import traceback
import argparse
import os
import shlex
import subprocess

from nephele2.pipelines.pipebase import PipeBase

SNAKEFILE_PATH = "/usr/local/src/nephele2/pipelines/picrust2_nephele/Snakefile"

class Picrust2Pipeline(PipeBase):
    def __init__(self, args):
        super().__init__(args)

class Picrust2PipelineError(Exception):
    pass

def main(args):
    try:
        exit_status = 0
        pipe = Picrust2Pipeline(args)
        logfile = pipe.outputs_dir + "logfile.txt"
        output_dir = pipe.outputs_dir + "picrust2_result/"
        tmp_dir = os.environ.get("TMPDIR", pipe.base_dir + "tmp/")
        num_cpus = os.cpu_count()

        pipe.log.info("Starting Picrust2 pipeline...")
        pipe.log.info(args)
        base_command = [
            f"snakemake -s {SNAKEFILE_PATH} -p --cores --keep-going",
            f"--directory {tmp_dir}"
        ]
        pipeline_config = [
            f"--config num_cpus={num_cpus} fasta_fp={args.fasta_fp.name} biom_fp={args.biom_fp.name} map_file={args.map_file.name}",
            f"output_dir={output_dir} tmp_dir={tmp_dir}",
            f"max_nsti={args.max_nsti} min_reads={args.min_reads} min_samples={args.min_samples} stratified={args.stratified}"
        ]
        base_command.extend(pipeline_config)

        main_pipeline = " ".join(base_command)
        pipe.log.info(main_pipeline)
        main_pipeline_args = shlex.split(main_pipeline)
        with open(logfile, "a") as f:
            result = subprocess.run(
                main_pipeline_args,
                stderr=f if args.stdout_redirect == "file" else None,
                stdout=f if args.stderr_redirect == "file" else None
            )
        if result.returncode == 1:
            error_message = (
                "There is an error while executing the pipeline. Please see logfile.txt for more information. "
                "Tip: search for keyword, i.e., Error in rule, to figure out which step(s) cause the failure."
            )
            raise Picrust2PipelineError(error_message)
    except Picrust2PipelineError as exc:
        exit_status = 1
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
        if exit_status == 1:
            pipe.log.error('Picrust2 pipeline completed with errors.')
        else:
            pipe.log.info('Picrust2 pipeline completed.')
        pipe.log.info(exit_status)
        exit(exit_status)
    

if __name__ == '__main__':
    PARSER = argparse.ArgumentParser()
    PARSER.add_argument("--job_id", type=str,
                        help="job_id when running in Nephele. Either job_id or inputs_dir should be specified.")
    PARSER.add_argument("-i", "--inputs_dir", type=str,
                        help="input directory for running outside Nephele. Either job_id or inputs_dir should be specified.")
    PARSER.add_argument("-o", "--outputs_dir", type=str,
                        help="output directory for running outside Nephele. optional")
    PARSER.add_argument("--fasta_fp", type=argparse.FileType('r'), required=True, help='Fasta file')
    PARSER.add_argument("--biom_fp", type=argparse.FileType('r'), required=True, help='Biom file')
    PARSER.add_argument("--max_nsti", type=int, default=2, help="Sequences with NSTI values above this value will be excluded")
    PARSER.add_argument("--min_reads", type=int, default=1, help="Minimum number of reads across all samples for each input ASV")
    PARSER.add_argument("--min_samples", type=int, default=1, help="Minimum number of samples that an ASV needs to be identfied within")
    PARSER.add_argument("--stratified", type=str, default="none", choices=["none", "stratified_only", "stratified_per_sequence_contrib"])
    PARSER.add_argument("--data_type", help="Ignore - placeholder.")
    PARSER.add_argument("--map_file", type=argparse.FileType('r'), required=True, help="Mapping file")
    PARSER.add_argument("--stdout_redirect", type = str, default = "file", choices = ["file", "std"], help = "choose where to redirect stdout. Using std to see it from terminal")
    PARSER.add_argument("--stderr_redirect", type = str, default = "file", choices = ["file", "std"], help = "choose where to redirect stderr. Using std to see it from terminal")
    ARGS = PARSER.parse_args()
    main(ARGS)
