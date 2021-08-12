#!/usr/bin/env python3

import traceback
import argparse
import os
import shlex
import subprocess
import shutil

from psutil import virtual_memory

from nephele2.pipelines.pipebase import PipeBase
from nephele2.infra.utils import s3utils
from nephele2.pipelines.wgsa import wgsa_config


class WgsaPipeline(PipeBase):
    def __init__(self, args):
        super().__init__(args)

class WgsaPipelineError(Exception):
    pass

def main(args):
    try:
        exit_status = 0
        pipe = WgsaPipeline(args)
        logfile = pipe.outputs_dir + "logfile.txt"
        tmp_dir = os.environ.get("TMPDIR", pipe.base_dir + "tmp/")
        num_cpus = os.cpu_count()
        # for bbmap (avoid heap memory issue)
        mem = virtual_memory()
        mem_gb = int(mem.total/1024/1024/1024)

        pipe.log.info("Starting WGSA pipeline...")
        pipe.log.info(args)
        base_command = [
            f"snakemake -s {wgsa_config.SNAKEFILE_PATH} -p --cores --keep-going --nocolor",
            f"--resources mem_gb={mem_gb}",
            f"--directory {tmp_dir}"
        ]
        pipeline_config = [
            f"--config num_cpus={num_cpus} decontaminate={args.decontaminate} map_file={pipe.map_fp}",
            f"input_dir={pipe.inputs_dir} output_dir={pipe.outputs_dir} bbmap_mem_gb={mem_gb} tmp_dir={tmp_dir}",
            f"trim_filter={args.trim_filter} include_ted_files={args.include_ted_files} ted_file_name={wgsa_config.TED_FILE_NAME}",
            f"biom_file_name={wgsa_config.BIOMFILE}"
        ]
        if args.trim_filter:
            pipeline_config.extend(
                [
                    f"average_read_quality={args.average_read_quality} min_read_length={args.min_read_length}",
                    f"trim_of_5={args.trim_of_5} trim_of_3={args.trim_of_3}"
                ]
            )
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
            raise WgsaPipelineError(error_message)
    except WgsaPipelineError as exc:
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
        # Copy output diagram if pipeline completed
        shutil.copy2(wgsa_config.OUTPUT_DIAGRAM_PATH, pipe.outputs_dir)
        # Transfer TED file to S3
        TED_FILE_PATH = pipe.outputs_dir + wgsa_config.TED_FILE_NAME
        if os.path.exists(TED_FILE_PATH):
            s3utils.transfer_file_to_s3(TED_FILE_PATH, pipe.job_id + "/" + wgsa_config.TED_FILE_NAME)
            os.remove(TED_FILE_PATH)

        if exit_status == 1:
            pipe.log.error('WGSA pipeline completed with errors.')
        else:
            pipe.log.info('WGSA pipeline completed.')
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
    PARSER.add_argument("--map_file", type=argparse.FileType('r'), required=True, help='Mapping file')
    PARSER.add_argument("--decontaminate", type=str, default="human", choices=["human", "house_mouse"])
    PARSER.add_argument("--trim_filter", action="store_true")
    PARSER.add_argument("--average_read_quality", type=int)
    PARSER.add_argument("--min_read_length", type=int)
    PARSER.add_argument("--trim_of_5", type=int)
    PARSER.add_argument("--trim_of_3", type=int)
    PARSER.add_argument("--include_ted_files", action="store_true")
    PARSER.add_argument("--data_type", help="Ignore - placeholder.")
    PARSER.add_argument("--stdout_redirect", type = str, default = "file", choices = ["file", "std"], help = "choose where to redirect stdout. Using std to see it from terminal")
    PARSER.add_argument("--stderr_redirect", type = str, default = "file", choices = ["file", "std"], help = "choose where to redirect stderr. Using std to see it from terminal")
    ARGS = PARSER.parse_args()
    main(ARGS)
