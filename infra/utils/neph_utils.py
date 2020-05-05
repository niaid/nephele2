import os
import json
import traceback
import argparse
from datetime import datetime, timedelta
from botocore.exceptions import ClientError
from nephele2 import config, tfvars
from nephele2.infra.utils import ec2, s3utils, fs_utils, sqs
from nephele2.infra.utils.email import email
from nephele2.rds.db_utils import DBUtils
from nephele2.infra.job import Job


class Serverd:
    """
    Code that runs on the daemonized on Server.
    """

    @staticmethod
    def manage_pending_q():
        """
        gets jobs that are due to be run,
        tries to start an EC2,
        if unable to start EC2 adds 30 mins to the backoff time
        (the time before which the job will be attempted again)
        Notifies the user,
        and then adds the job back to the pending queue.
        """
        try:
            jobs = sqs.receive_pending()
            for job in jobs:
                try:
                    instance = Serverd._launch_job(job.job_id)
                except ClientError as err:
                    backoff = datetime.now() + timedelta(minutes=30)
                    if not job.user_notified:
                        email.notify_user_job_queued(job.job_id)
                    sqs.add_to_pending(job.job_id,
                                       user_notified='true',
                                       backoff=backoff)
                    N2Manager.notify_admin_neph_exception(err)
                if instance:
                    DBUtils.set_machine_info(job.job_id,
                                             instance.id,
                                             instance.image_id,
                                             instance.instance_type)
                    email.send_started_email(job.job_id)
        except ClientError as err:
            N2Manager.notify_admin_neph_exception(err)

    @staticmethod
    def _launch_job(job_id):
        """
        start a new worker instance & record machine details in DB
        mark job as XXXX in db (PROCESSING??)
        """
        job = DBUtils.get_job(job_id)
        if not job:
            err = ("Failure at _start_ec2: no such job {}".format(job_id))
            N2Manager.notify_admin_neph_exception(err)
            return None
        if job.status == 'Pre-Processing':
            # implies that the job has already started being processed
            # below.
            return None
        job_details = DBUtils.get_job_type_details(job_id)
        instance = ec2.start_worker_EC2(
            job_id,
            job_details['ami_id'],
            job_details['default_instance_type'])
        # ONLY set to pre Processing once we have an instance.
        # otherwise jobs put back into pending are ignored.
        DBUtils.set_job_status(job_id, 'Pre-Processing')
        return instance

    @staticmethod
    def rm_get_old_unsubmitted_jobs():
        """
        DAEMONIZED ON SERVER MACHINE.
        This is currently not doing anything.
        RM's jobs from DB
        """
        print('Running rm_get_old_unsubmitted_jobs()')
        job_ids = DBUtils.get_old_unsubmitted_jobs()
        for job_id in job_ids:
            DBUtils.delete_job(job_id)

    @staticmethod
    def rm_expired_subsribers():
        pass
        """
        DAEMONIZED ON SERVER MACHINE.
        A user initially adds their name to the DB, and once they
        reply to an email their address is considered valid. If they fail
        to validate within 1 day, the address is removed."""
        #  print('Running rm_expired_subsribers()')
        #  try:
        #      users = DBUtils.get_old_unconfirmed_users()
        #      for addr in users:
        #          DBUtils.delete_user(addr)
        #  except Exception as unknown_expn:
        #      print(unknown_expn)
        #      N2Manager.notify_admin_neph_exception(unknown_expn)

    @staticmethod
    def lookup_efs_mnt(path):
        """lookup_efs_mnt

        sends warning if vital file path doesn't exist

        :param fp:
        """
        try:
            return os.listdir(path)
        except FileNotFoundError as fnf:
            N2Manager.notify_admin_neph_exception(fnf)

    @staticmethod
    def _job_expired_in_db(job_id):
        two_hours_ago = datetime.now() - timedelta(hours=2)
        day_ago = datetime.now() - timedelta(hours=24)
        job = DBUtils.get_job(job_id)
        if not job:
            print('Unable to find DB Entry for dir: ' + job_id)
            return False
        too_long_init = (job.status == "Initializing" and
                         job.created < two_hours_ago)
        day_old_not_rmd = (job.status in ['Failed', 'Succeeded'] and
                           job.completed < day_ago)

        if too_long_init or day_old_not_rmd:
            return True
        return False

    @staticmethod
    def rm_old_stale_efs_dirs():
        """
        DAEMONIZED ON SERVER MACHINE.
        Removes dead job folders from EFS.
        rm anything that's been INIT-ing for > 2 hours.
        """
        try:
            job_dirs = Serverd.lookup_efs_mnt(config.UPLOAD_PATH)
            if job_dirs:
                for job_dir in job_dirs:
                    if Serverd._job_expired_in_db(job_dir):
                        fs_utils.remove_efs_dir(job_dir)
        except Exception as ex:
            N2Manager.notify_admin_neph_exception(ex)


class N2Manager:
    """
    interface for nephele2

    runs by view:
    -------------
    init_job(user_id)  # used in views, returns job_id, uniq dname in efs &
    adds job to db

    get_job_type(job_id)  # eg 'mothur miseq' etc. return job_type (see
    job_type table)

    get_job_type_details(job_id)  # looks up details about the job type
    associated with the job

    get_default_args(job_name)  # looks up defulat args for some *job type*

    get_default_job_arguments(input_type)    # look up details of pipeline by
    *input data type*

    submit_job(job_id)          # check job ok, adds to pending. used in views,


    """

    @staticmethod
    def get_job_base_dir(job_id):
        """used in views indirectly."""
        return config.UPLOAD_PATH + job_id + '/'

    @staticmethod
    def get_job_inputs_d(job_id):
        """Used in views."""
        return N2Manager.get_job_base_dir(job_id) + 'inputs/'

    @staticmethod
    def get_job_outputs_d(job_id):
        """Used in views."""
        return N2Manager.get_job_base_dir(job_id) + 'outputs/'

    @staticmethod
    def notify_admin_neph_exception(exception):
        """
        A failure that's greater than a job failure. Indicates something
        systemic is wrong.
        Send email to admins.
        """
        from nephele2 import nephele
        if nephele.app.email_timer.ready_to_send():
            try:
                email.send_infra_failure_email(str(exception))
                print('aws_infrastructure_failure')
                print(str(exception))
            except Exception:
                # at this point everything has failed, log and carry on...
                print(traceback.format_exc())

    @staticmethod
    def get_authenticated_user(email):
        """
        Returns the user associated with the given email address.
        """
        return DBUtils.get_user_by_email(email)

    @staticmethod
    def get_user(email):
        """
        Returns the user associated with the given email address.
        """
        return DBUtils.get_user_by_email(email)

    @staticmethod
    def create_user(email, fname, lname, affiliation, affiliation_category,
                    ref, subscribe, analysis=None):
        """create_user

        :param email:
        :param fname:
        :param lname:
        :param affiliation:
        :param affiliation_category:
        :param ref:
        :param subscribe:
        :param analysis:
        """
        DBUtils.create_user(email,
                            fname,
                            lname,
                            affiliation,
                            affiliation_category,
                            ref,
                            subscribe,
                            analysis=analysis)

    @staticmethod
    def confirm_user(user_email, is_confirmed):
        """confirm_user

        :param user_email:
        :param is_confirmed:
        """
        DBUtils.set_confirm_email(user_email, is_confirmed)

    @staticmethod
    def init_job(user_id):
        """
        TODO CHECKS user_id is valid
        Creates a job id
        Creates the EFS location to store data (using job id)"""
        job_id = fs_utils.gen_uniq_dname(config.UPLOAD_PATH)
        DBUtils.create_job(user_id, job_id)
        return job_id

    @staticmethod
    def submit_job(job_id, job_name, job_desc, job_args):
        # form.instance_type <- FIXME: we need to pass this in and do something
        # with it
        """
        Sets the job details in the database and puts the job ID in the
        pending queue to trigger job start.  Raises an exception on error.

        Args:
            job_id (str): a valid job ID
            job_name (str): a valid job type
            job_args (JSON): a dict or JSON str of the job arguments

        Raises:
        """
        # check job looks ok and try to add it to pending Q
        if isinstance(job_args, dict):
            job_args = json.dumps(job_args)
        DBUtils.set_job_type(job_id, job_name)
        DBUtils.set_job_description(job_id, job_desc)
        DBUtils.set_job_arguments(job_id, job_args)
        DBUtils.set_job_status(job_id, "Pending")
        sqs.add_to_pending(job_id)

    @staticmethod
    def get_map_file_template(data_type):
        """get_map_file_template

        :param data_type:
        """
        if data_type == "SE":
            key = "fastq_single_end_mapping_template_qiime.xlsx"
        elif data_type == "WGS_PE":
            key = "WGS_bioBakery_template_N2.xlsx"
        elif data_type == "WGS_SE":
            key = "WGS_bioBakery_template_single_end_N2.xlsx"
        elif data_type == "QC_PE":
            key = "QC_pipeline_template_paired_end_N2.xlsx"
        elif data_type == "QC_SE":
            key = "QC_pipeline_template_single_end_N2.xlsx"
        else:
            # default, in case data_type is None
            # will also be used for PE
            key = "Demultiplexed_QIIME_mothur_DADA2_template.xlsx"
        return s3utils.get_presigned_url(
            tfvars.RESOURCES_BUCKET_ID, key)

    @staticmethod
    def get_resource_files(filename):
        """get_resource_files

        :param filename:
        """
        return s3utils.get_presigned_url(
            tfvars.RESOURCES_BUCKET_ID, filename)

    @staticmethod
    def get_test_files(filepath):
        """get_test_files

        :param filepath:
        """
        return s3utils.get_presigned_url(
            tfvars.TESTDATA_BUCKET_ID, filepath)

    @staticmethod
    def send_registration_email(user_email, confirm_url):
        """
        Sends an email confirmation email to the user.

        Args:
            user_email (str): the email address to send the message to
        """
        try:
            #  salted_email = N2Manager.salt_string(user_email)
            email.send_registration_email(user_email, confirm_url)
        except ClientError as aws_err:
            N2Manager.notify_admin_neph_exception(aws_err)
        except Exception as unknown_expn:
            N2Manager.notify_admin_neph_exception(unknown_expn)

    @staticmethod
    def save_user_file(uniq_dname, file):
        """
        Downloads the contents of the file from user's local to the
        specified location in the local system.

        Args:
            uniq_dname (str): the directory to place the file in.
            file (str): the file handle to the file in the user's local system.

        Returns:
            str: the path to the directory where the content was downloaded.
        """
        from werkzeug.utils import secure_filename
        filename = secure_filename(file.filename)
        filepath = os.path.join(uniq_dname, filename)
        file.save(filepath)
        # set files to read only
        os.chmod(filepath, 0o666)
        return filepath

    @staticmethod
    def load_job_type_table():
        """load_job_type_table"""
        DBUtils.load_job_types()

    @staticmethod
    def retrieve_prev_job(user_id, prev_id):
        """retrieve_prev_job

        called by transfer_data, via JS.
        At this point we've already ensured the job is available on S3. It
        should not fail, but if it does for whatever reason we should notify
        admin. This is the reason for the broad Exception catch.

        :param user_id:
        :param prev_id:
        """
        try:
            job_id = N2Manager.init_job(user_id)
            job = Job(job_id)
            job.retrieve_prev_job(prev_id)
            return job_id
        except Exception as e:
            N2Manager.notify_admin_neph_exception(e)
            raise e

    @staticmethod
    def start_polling_pending_q():
        """start_polling_pending_q

        Called by poll_pending_q_d.service.
        Catching all exceptions here since this is a daemon running on remote
        server. Should never happen, but need to notitfy admin if does."""
        while True:
            try:
                Serverd.manage_pending_q()
            except Exception as unknown_expn:
                N2Manager.notify_admin_neph_exception(unknown_expn)


def main(args):

    if args.poll_pending_q:
        N2Manager.start_polling_pending_q()

    if args.rm_non_subscribed_users is True:
        Serverd.rm_expired_subsribers()

    if args.rm_old_unsubmitted_jobs is True:
        Serverd.rm_get_old_unsubmitted_jobs()

    if args.clean_efs is True:
        Serverd.rm_old_stale_efs_dirs()

    if args.load_job_type:
        N2Manager.load_job_type_table()


if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(description='CLI Interface to N2.')
    PARSER.add_argument(
        "--load_job_type",
        help="Loads contents of rds/conf.yaml to DB.",
        action="store_true")
    PARSER.add_argument(
        "--poll_pending_q",
        help="Starts polling the pending queue.",
        action="store_true")
    PARSER.add_argument(
        "--rm_non_subscribed_users",
        help="Removes users who never finish subscribing",
        action="store_true")
    PARSER.add_argument(
        "--rm_old_unsubmitted_jobs",
        help="Removes jobs which are never submitted.",
        action="store_true")
    PARSER.add_argument(
        "--clean_efs", help="Flag to clean EFS.", action="store_true")
    ARGS = PARSER.parse_args()
    main(ARGS)
