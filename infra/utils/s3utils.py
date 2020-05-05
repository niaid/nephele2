#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import tarfile
import os
import boto3
from botocore.exceptions import ClientError
import humanize
from nephele2 import tfvars, config


S3_CLIENT = boto3.client('s3', region_name='us-east-1')

def get_obj_expiration(bucket, key):
    try:
        header = S3_CLIENT.head_object(Bucket=bucket, Key=key)
        expiry_info = header.get('Expiration')
        # NOTE: if AWS syntax ever changes, this regex may break
        match = re.search('expiry-date=\".*GMT\"', expiry_info)
        date_match = re.search('\"(.*GMT)\"', match.group(0))
        expiry_date = date_match.group(0)
        return re.sub(r'"', '', expiry_date)
    except ClientError:
        return None


def prefix_exists(bucket, key):
    """ensure_object_exists

    :param bucket:
    :param key:
    """
    objects = S3_CLIENT.list_objects(Bucket=bucket, Prefix=key)
    contents = objects.get('Contents')
    if contents:
        return len(contents) > 0
    return False


def ensure_object_exists(bucket, key):
    """ensure_object_exists

    :param bucket:
    :param key:
    """
    try:
        return S3_CLIENT.head_object(Bucket=bucket, Key=key)
    except ClientError as e:
        # If a client error is thrown, then check that it was a 404 error.
        # If it was a 404 error, then the bucket does not exist.
        if e.response['ResponseMetadata']['HTTPStatusCode'] == '404':
            return False
        else:
            # TODO notify and then return false
            # if not 404 then something is weird
            #  N2Manager.notify_admin_neph_exception(e)
            return False


def get_presigned_url(bucket, key):
    # check for existence because AWS will return a URL for non-existent
    # objects this is not a bug, it's a feature that allows you to use
    # pre-signed URLs for uploads
    return S3_CLIENT.generate_presigned_url(
        ClientMethod='get_object',
        Params={'Bucket': bucket, 'Key': key})


def get_file_size(bucket, key):
    """get_file_size

    :param bucket:
    :param key:
    """
    header = ensure_object_exists(bucket, key)
    if 'ContentLength' not in header:
        return None
    byte_size = header.get('ContentLength')
    if not byte_size:
        return None
    return humanize.naturalsize(byte_size)


def retrieve_prev_job_inputs(new_job_id, orig_job_id):
    """retrieve_prev_job_inputs
    gets the inputs of a previous run local to the EC2
    note, previous inputs may no long be available.

    :param new_job_id:
    :param orig_job_id:
    """
    transfer = boto3.s3.transfer.S3Transfer(S3_CLIENT)
    paginator = S3_CLIENT.get_paginator('list_objects')
    page_iterator = paginator.paginate(Bucket=tfvars.S3_JOBS_BUCKET,
                                       Prefix=orig_job_id+'/inputs')
    new_inputs_loc = config.UPLOAD_PATH + new_job_id + '/inputs/'
    for page in page_iterator:
        contents = page.get('Contents')
        if not contents:
            return
        for content in contents:
            key = content.get('Key')
            fname = os.path.basename(key)
            transfer.download_file(tfvars.S3_JOBS_BUCKET,
                                   key,
                                   new_inputs_loc+fname)


def _extract_archive(archive_name, path):
    """_extract_archive

    :param archive_name:
    :param path:
    """
    archive = tarfile.open(archive_name, 'r:gz')
    archive.extractall(path)


def retrieve_prev_job_inputs_gz(new_job_id, orig_job_id):
    """retrieve_prev_job_inputs
    delete this function 90 days after whole input dir ship
    gets the inputs of a previous run local to the EC2
    need to check if previous run inputs were compressed
    also, previous inputs may no long be available.
    If file not found, ignore and carry on
    (will look for dir)

    :param new_job_id:
    :param orig_job_id:
    """
    transfer = boto3.s3.transfer.S3Transfer(S3_CLIENT)
    zipped_inputs = "inputs.{id}.tar.gz".format(id=orig_job_id)
    new_inputs_loc = config.UPLOAD_PATH + new_job_id + '/'\
        + zipped_inputs
    transfer.download_file(tfvars.S3_JOBS_BUCKET,
                           orig_job_id+'/'+zipped_inputs,
                           new_inputs_loc)
    _extract_archive(new_inputs_loc,
                     config.UPLOAD_PATH + new_job_id)


def ensure_job_exists_in_s3(job_id):
    """_ensure_job_exists_in_s3

    :param job_id:
    """
    try:
        inputs = S3_CLIENT.list_objects(Bucket=tfvars.S3_JOBS_BUCKET,
                                        Prefix=job_id+'/inputs')
        if inputs.get('Contents'):
            return True
        return False
    except ClientError:
        # assume that if an error is thrown it's a 404.
        # in any event, data cannot be found.
        # TODO send infra fail if ClientError is NOT a 404
        #  if e.response['ResponseMetadata']['HTTPStatusCode']) == '404':
        return False

def transfer_dir_to_s3(dname, key_base):
    transfer = boto3.s3.transfer.S3Transfer(S3_CLIENT)
    for fname in os.listdir(dname):
        key = key_base + fname
        transfer.upload_file(dname+fname, tfvars.S3_JOBS_BUCKET, key)

def transfer_file_to_s3(fname, key):
    transfer = boto3.s3.transfer.S3Transfer(S3_CLIENT)
    transfer.upload_file(fname, tfvars.S3_JOBS_BUCKET, key)
