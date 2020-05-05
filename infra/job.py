#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Stores info about a nephele job."""

import os
import botocore
from nephele2 import tfvars, config
from nephele2.infra.utils import s3utils, fs_utils


class Job():
    """
    The output of a neph pipe ends up in S3.
    Record info about this output here!
    """

    def __init__(self, job_id):
        self._job_id = job_id
        self._bucket = tfvars.S3_JOBS_BUCKET
        self._prefix = config.UPLOAD_PATH + job_id
        self._results_file = self._prefix+"/PipelineResults."+job_id+".tar.gz"
        self._inputs = self._prefix+"/inputs/"
        self._outputs = self._prefix+"/outputs/"
        self._log_file = self._outputs+"logfile.txt"
        self._s3_results = self.id + '/' + os.path.basename(self._results_file)
        self._s3_log = self.id + '/' + os.path.basename(self.log_file)
        self._s3_inputs = self.id + '/inputs/'

    @property
    def id(self):
        """id
        returns the job id"""
        return self._job_id

    @property
    def outputs(self):
        """outputs
        the job output dir name"""
        return self._results_file

    @property
    def results_file(self):
        """results_file
        the compressed outputs"""
        return self._results_file

    @property
    def inputs(self):
        """inputs
        inputs dir"""
        return self._inputs

    @property
    def log_file(self):
        """log_file
        path to log file"""
        return self._log_file

    @property
    def bucket(self):
        """bucket
        the S3 bucket the results will be transferred to"""
        return self._bucket

    def exists(self):
        """
        Checks to ensure that results files exists.
        """
        return s3utils.prefix_exists(self.bucket, self.id)

    def retrieve_prev_job(self, prev_id):
        """retrieve_prev_job

        :param prev_id:
        """
        try:
            s3utils.retrieve_prev_job_inputs_gz(self.id, prev_id)
        except botocore.exceptions.ClientError:
            pass
        s3utils.retrieve_prev_job_inputs(self.id, prev_id)

    def get_job_expiration(self):
        """get_job_expiration
        returns an exp date for job."""
        date = s3utils.get_obj_expiration(self.bucket, self._s3_results)
        if not date:
            date = '?'
        return date

    def get_logfile_url(self):
        """get_logfile_url
        returns a signed URL for logfile"""
        return s3utils.get_presigned_url(self.bucket, self._s3_log)

    def get_results_url(self):
        """get_results_url
        returns a signed URL for results"""
        return s3utils.get_presigned_url(self.bucket, self._s3_results)

    def get_results_size(self):
        """get_results_size
        returns size of result set."""
        size = s3utils.get_file_size(self.bucket, self._s3_results)
        if not size:
            size = '?'
        return size

    def transfer_inputs(self):
        """
        Does an S3 Transfer, job data from EFS to S3.
        uses key_base to create the path in the S3 bucket.
        """
        s3utils.transfer_dir_to_s3(self.inputs, self._s3_inputs)

    def compress_results(self):
        """generates standard Nephele Outputs. Called at end of pipeline."""
        fs_utils.compress_dir(self.results_file, self._prefix, 'outputs')

    def transfer_results(self):
        """
        Does an S3 Transfer, job data from EFS to S3.
        If the transfer fails, sets the job status to failure in the DB.
        """
        s3utils.transfer_file_to_s3(self.log_file, self._s3_log)
        s3utils.transfer_file_to_s3(self.results_file, self._s3_results)

    def remove_efs_dir(self):
        """
        For some reason NFS can bork upon big rm -rf's, and you will see an
        error similar to:
        cannot remove: Directory not empty
        or:
        /bin/rm: cannot remove':
        Device or resource busy
        This attempts to just catch those errors, delay x*2 seconds (?), retry
        the delete.
        If it KEEPS failing THEN it will raise expn
        """
        fs_utils.remove_efs_dir(self._prefix)
