#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Utilities for interacting with the local file system"""
import os
import uuid
import time
import signal
from ftplib import FTP
import sh
from nephele2 import config


def curl_check_url(url, target_dir):
    """
    tries to get .listing file at url and write to target_dir.

    Args:
        url (str): the URL location of the remote file directory
        target_dir (str): the filepath to write the listing file to

    Returns:
        str: the .listing filepath or None if one wasn't generated

    Raises:
        Exception: whatever exception was raised when things went wrong
    """
    listing_fname = target_dir+"/.listing"
    sh.curl('-l', url, _out=listing_fname)
    if os.path.isfile(listing_fname):
        return listing_fname
    return None


def gen_uniq_dname(base):
    """
    Generates a uuid to use for a job ID and creates a directory with
    the uuid as it's name.

    Args:
        base (str): the base filepath where we want to create the job dir

    Returns:
        str: a new job ID
    """
    while True:
        job_id = str(uuid.uuid4()).split('-')[-1]
        base_d = base + job_id
        in_d = base_d + '/inputs'
        out_d = base_d + '/outputs'
        # TODO: this only makes sure we don't have an active job with this
        # id, check the db instead?
        if not os.path.lexists(base_d):
            sh.mkdir('--mode=u+rwx,g+rs,g+rwx,o+rx', base_d)
            #  sh.chown('www-data:www-data', base_d)
            sh.mkdir('-p', in_d)
            sh.mkdir('-p', out_d)
            #  sh.chown('-R', 'www-data:www-data', base_d)
            return job_id


def compress_dir(compressed_fname, prefix, dname):
    """compress_d is something like inputs
    fname_base is something like inputs.{job_id}
    base_d is something like /mnt/EFS/user_uploads/{job_id}/
    The command uses this trick:
    sh.tar('-czf', '/mnt/EFS/user_uploads/tmp/outs.tar.gz', '-C', \
    '/mnt/EFS/user_uploads/tmp/', 'outputs/')
    where you -C change into the dir, spec the dir you want to tar,
    then create the fname outs.tar.gz"""

    # rm .gz from compressed_fname for tar name
    tar_fname = compressed_fname[:-3]
    sh.tar('-cf', tar_fname, '-C', prefix, dname)
    if os.path.exists(compressed_fname):
        sh.rm(compressed_fname)
    sh.pigz(tar_fname)


def remove_efs_dir(dname):
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
    max_num_rm_retrys = 3
    i = 0
    while True:
        try:
            sh.rm('-rf', dname)
            break
        except sh.ErrorReturnCode_1:
            i += 1
            if i > max_num_rm_retrys:
                raise
            time.sleep(2)


def chk_efs():
    """chk_efs
    The EFS can hang up at times. This function tries to sense if this
    has occured, and raise an IO error if it occurs"""
    def efs_excpn_handler(signum, frame):
        raise IOError('EFS Unreachable')
    # job_id
    num_seconds_wait = 1
    try:
        signal.signal(signal.SIGALRM, efs_excpn_handler)
        signal.alarm(num_seconds_wait)
        #  + job_id
        if not os.path.isdir(config.UPLOAD_PATH):
            raise IOError('EFS Unreachable')
    finally:
        signal.alarm(0)


def ftp_get_files(srvr_and_dname, files, target_dir):
    """
    gets a list of files from FTP URI srvr_and_dname, places into
    target_dir

    Args:
        srvr_and_dname (str): eg ftp://helix.nih.gov/pub/bcbb/sjogren
        files (list): list of string filenames to look for
        target_dir (str): where to copy the files above to

    Raises:
        Exception: FTP exception, filesystem write exception
    """
    if srvr_and_dname.lower().startswith('ftp://'):
        srvr_and_dname = srvr_and_dname[6:]
    if not target_dir.endswith('/'):
        target_dir += '/'
    ftp_loc = srvr_and_dname.split('/')
    ftp_dir = ''
    if ftp_loc:
        ftp_dir = '/'.join(ftp_loc[1:])
    with FTP(ftp_loc[0]) as ftp:
        ftp.login()
        ftp.cwd(ftp_dir)
        for file in files:
            # sleep some random amount of time to keep CIT happy.
            time.sleep(10)
            local_fname = target_dir + file
            retr_str = 'RETR ' + file
            ftp.retrbinary(retr_str, open(local_fname, 'wb').write)
            os.chmod(local_fname, 0o666)
