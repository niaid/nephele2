#!/usr/bin/env python3

import os
import sys
import time
import argparse
import subprocess
import csv
import sh
from nephele2 import config
from nephele2.infra.utils.email import email
from nephele2.infra.utils import fs_utils
from nephele2.infra.job import Job
from nephele2.rds.db_utils import DBUtils


def get_file_list(map_fp):
    """creates a list of Sample objects"""
    samples = list()
    reader = csv.DictReader(open(map_fp), delimiter='\t')
    for row in reader:
        if row['#SampleID'] == '':
            continue
        if row.get('ForwardFastqFile'):
            samples.append(row['ForwardFastqFile'])
        if row.get('ReverseFastqFile'):
            samples.append(row['ReverseFastqFile'])
    return samples


def process_job_on_worker(job_id):
    """
    DAEMONIZED ON WORKER
    sets job to Running in DB
    gets job arguments from DB (& converts them to --arg_name format)
    starts job and waits for result
    marks that job result in DB (Succeeded || Failed)
    runs compress_results
    transfer data to S3
    rm results
    Note, if an exception is thrown by local_exec, it will be logged in the
    DB
    """
    job = Job(job_id)
    try:
        # learn script name, args, mark job as Running
        DBUtils.set_job_status(job_id, 'Running')
        args = DBUtils.get_job_arguments(job_id)
        script_name = DBUtils.get_script_name(job_id)

        # TODO this shouldn't be here - move to job?
        if 'ftp' in args.keys():
            files = get_file_list(args['map_file'])
            fs_utils.ftp_get_files(args['ftp'], files, job.inputs)
            del args['ftp']
        job.transfer_inputs()
        cli_args = _argify_dict(args)

        script_path = config.PIPELINES_LOC_ON_WRKR + script_name
        _local_exec(job_id, script_path, cli_args)
        # mark job as Succeeded. If fails to transfer, overwrite status.
        DBUtils.set_job_status(job_id, "Succeeded")
    except Exception as excpn:
        DBUtils.set_job_status(job_id, "Failed", stack_trace=str(excpn))
        email.send_infra_failure_email(str(excpn), job_id=job_id)

    finally:
        try:
            job.compress_results()
            job.transfer_results()
            DBUtils.set_data_transferred(job_id, True)
        except Exception as excpn:
            DBUtils.set_job_status(job_id, "Failed", stack_trace=str(excpn))
            email.send_infra_failure_email(str(excpn), job_id=job_id)
        try:
            job.remove_efs_dir()
        except Exception as excpn:
            email.send_infra_failure_email(str(excpn), job_id=job_id)
        email.notify_user_job_completion(job_id)


def _local_exec(job_id, script_name, args):
    """_local_exec
    :param job_id:
    :param script_name:
    :param args:
    """
    args.insert(0, script_name)
    # subprocess.run returns a CompletedProcess obj, but need to wire
    # it to grab outs I think we only need to capture the stderr, but
    # if we want stdout too add: stdout=subprocess.PIPE,
    result = subprocess.run(args, stderr=subprocess.PIPE)
    if result.returncode != 0:
        raise RuntimeError(result.stderr)


def _argify_dict(dictnry):
    """
    Used in calling scripts with --arg_name arg type args
    converts a dict, eg:
    { "criteria": "90", "database": "SILVA", "body_site": "Saliva",
    "nearest_n_samples": "7", "hmp_database": "SILVA_99", "optimize":
    "start-end", "region_dacc": "v1v3", "maxlength": "300",
    "comp_with_dacc": "False", "difference_rank": "2"} to a list, eg:
    ['--region_dacc v1v3', '--hmp_database SILVA_99', '--body_site Saliva',
    '--database SILVA', '--criteria 90', '--difference_rank 2',
    '--maxlength 300', '--optimize start-end', '--nearest_n_samples 7',
    '--comp_with_dacc False'] """
    arg_list = list()
    for key, val in dictnry.items():
        if val is True:
            arg_list.append('--' + key.lower())
        elif val is not False:
            arg_list.append('--' + key.lower())
            arg_list.append(str(val))
    return arg_list


def main(args):
    """main
    entry point for worker.

    :param args:
    """
    JOB_ID = os.environ.get('JOB_ID')
    if JOB_ID is None:
        raise Exception('EC2 started without Job ID - exiting now.')

    if args.process_job:
        try:
            print('Received job ID, starting job.')
            sys.stdout.flush()
            process_job_on_worker(JOB_ID)

        except Exception as excpn:
            email.send_infra_failure_email(str(excpn), job_id=JOB_ID)

        finally:
            print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
            print('SHUTDOWN SOON!')
            print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
            sys.stdout.flush()
            time.sleep(120)

    elif args.chk_efs:
        try:
            print('Checking EFS...')
            fs_utils.chk_efs()
        except IOError as ioerror:
            email.send_infra_failure_email(str(ioerror), job_id=JOB_ID)
            print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
            print('SHUTDOWN SOON!')
            print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
            sys.stdout.flush()
            time.sleep(120)
            sh.shutdown('-h', 'now')


if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(description='CLI Interface to ')
    PARSER.add_argument("--process_job", help="Starts job.",
                        action="store_true")
    PARSER.add_argument("--chk_efs",
                        help='Asks if the EFS is reachable. Times out if not.',
                        action="store_true")
    ARGS = PARSER.parse_args()
    main(ARGS)
