"""
Class for updating and retrieving data from the DB.
"""

import json
from datetime import datetime, timedelta
import yaml
import sqlalchemy
from nephele2.rds.db_models import UserEmail, User, MachineInfo, Job,\
    JobType
from nephele2 import tfvars
import nephele2


class DBUtils(object):
    """
    Contains:
    load_job_types()
    create_user()
    get_user_by_job_id()
    get_user_by_email()
    get_old_unconfirmed_users()
    set_confirm_email()
    set_bad_email()
    get_subscribed_user_email()
    delete_user()
    create_job()
    get_job()
    get_machine_info_by_job_id()
    get_machine_info()
    set_machine_info()
    set_job_type()
    set_job_arguments()
    get_data_location()
    get_job_id()

    """
    @staticmethod
    def delete_job(job_id):
        job = DBUtils.get_job(job_id)
        if job:
            with nephele2.db_write() as db:
                db.delete(job)

    @staticmethod
    def load_job_types():
        """ Loads the JobType table, which will be considered read-only by the
        system. This data must be present in order for the application to
        function.

        Args:
            db (SQLClient): a SQLAlchemy database connection with write
            permissions.
        """
        with open('rds/conf.yaml') as yaml_file:
            conf = yaml.safe_load(yaml_file)
        with nephele2.db_write() as db:
            for entry in conf['job_types']:
                try:
                    db.query(JobType).filter_by(name=entry['name']).one()
                    db.query(JobType).filter_by(
                        name=entry['name']).update(entry)
                except sqlalchemy.orm.exc.NoResultFound:
                    db.add(JobType(entry))

    @staticmethod
    def create_user(email_addr,
                    fname,
                    lname,
                    affiliation,
                    affiliation_category,
                    ref,
                    subscribe,
                    analysis=None):
        """
        Creates a new user, adding rows to the users, user_email, and user_info
        tables.

        Args:
            db (SQLClient): a SQLAlchemy database write connection
            email_addr (str): a validly formed email address
            fname (str): the user's first name
            lname (str): the user's last name
            affiliation (str): name of the user's employer
            affiliation_category (str): one of: NIH, Government (non-NIH),
            University, Research Institute (Private), Citizen Scientist
            (no-affiliation), Other ref (str): who referred the user. One of:
            NA, NIH Bioinformatics listserve, NIH Microbiome listserve,
            Internet, Colleague, Other subscribe (bool): indicates whether the
            user wants to subscribe to the newsletter analysis (str): the type
            of analysis the user intends to use the system for

        """
        # Apparently the database considers an empty string to not be NULL,
        # but we don't want to allow blank strings in the address field in
        # the DB, so make sure "empty" strings are considered NULL.

        # subscribe will come in as a boolean, but needs to be stored as an int
        if subscribe:
            subscribe = 1
        else:
            subscribe = 0

        user_data = {
            'first_name': fname,
            'last_name': lname,
            'user_address': {
                'address': email_addr,
                'subscribed': subscribe
            },
            'user_info': {
                'affiliation': affiliation,
                'affiliation_category': affiliation_category,
                'referrer': ref,
                'analysis_type': analysis
            }
        }
        with nephele2.db_write() as db:
            db.add(User(user_data))

    @staticmethod
    def get_user_by_job_id(job_id):
        """
        Returns a User object associated with the job_id.

        Args:
            db (SQLClient): a SQLAlchemy database read connection
            job_id (str): a valid job ID

        Returns:
            :obj:`User`: a User object

        """
        with nephele2.db_read() as db:
            try:
                job = db.query(Job).filter_by(job_id=job_id).one()
                return job.user
            except sqlalchemy.orm.exc.NoResultFound:
                return None

    @staticmethod
    def get_user_and_email_by_email(email_addr):
        with nephele2.db_read() as db:
            try:
                user = db.query(User).filter_by(email_address=email_addr).one()
                user_email = user.user_address
                return (user, user_email)
            except sqlalchemy.orm.exc.NoResultFound:
                return (None, None)

    @staticmethod
    def get_user_by_email(email_addr):
        """
        Retrieves the user by email address.

        Args:
            db (SQLClient): a SQLAlchemy database read connection
            email_addr (str): a validly formatted email address

        Returns:
            :obj:`User`: a User object

        """
        with nephele2.db_read() as db:
            try:
                return db.query(User).filter_by(email_address=email_addr).one()
            except sqlalchemy.orm.exc.NoResultFound:
                return None

    @staticmethod
    def get_old_unconfirmed_users():
        """Looks up users who registered > a day ago and have not confirmed.

        Args:
            db (SQLClient): a SQLAlchemy database read connection

        Returns:
            list(str): list of bad emails.

        Raises:
            Exception: any exception thrown
        """
        yesterday = datetime.now() - timedelta(hours=24)
        with nephele2.db_read() as db:
            users = db.query(UserEmail).filter(
                UserEmail.is_confirmed.is_(False),
                UserEmail.subscr_date < yesterday)
            return [u.address for u in users.all()]

    @staticmethod
    def set_confirm_email(email_addr, is_confirmed):
        """
        Sets the confirmed flag on a UserEmail.

        Args:
            db (SQLClient): a SQLAlchemy database write connection
            email_addr (str): an email address, must exist in the database or
            the method will throw a NepheleError.NepheleRowNotFound exception
            is_confirmed (bool): True if email_addr is confirmed, False
            otherwise

        """
        # Booleans are represented as ints in the DB, so we need to convert
        if is_confirmed:
            is_confirmed = 1
        else:
            is_confirmed = 0

        # TODO : possibly want to raise error if no such email_addr
        with nephele2.db_write() as db:
            try:
                user_email = db.query(UserEmail).filter_by(
                    address=email_addr).one()
                user_email.is_confirmed = is_confirmed
                db.add(user_email)
            except sqlalchemy.orm.exc.NoResultFound:
                return

    @staticmethod
    def set_bad_email(email_addr):
        """
        Sets the is_bad flag on a UserEmail.
        called by lambda
        Args:
            db (SQLClient): a SQLAlchemy database write connection
            email_addr (str): an email address, must exist in the database or
            the method will throw a NepheleError.NepheleRowNotFound exception
            is_bad (bool): True if the email address is bad, False otherwise

        """
        with nephele2.db_write() as db:
            try:
                user_email = db.query(UserEmail).filter_by(
                    address=email_addr).one()
                user_email.is_bad = 1
                db.add(user_email)
            except sqlalchemy.orm.exc.NoResultFound:
                return

    @staticmethod
    def get_subscribed_user_email():
        """
        Gets a list of all confirmed email address, that aren't flagged as bad,
        that are subscribed to the newsletter.

        Args:
            db (SQLClient): a SQLAlchemy database read connection

        Returns:
            list(str): a list of email addresses

        Raises:
        """
        with nephele2.db_read() as db:
            users = db.query(UserEmail).filter(
                UserEmail.subscribed.is_(True),
                UserEmail.is_confirmed.is_(True),
                UserEmail.is_bad.isnot(True)).all()
            addresses = []
            for user in users:
                addresses.append(user.address)
            return addresses

    @staticmethod
    def delete_user(email_addr):
        pass
        """
        Deletes the user identified by email_addr. Will also delete all user
        info and email info for the email_addr. The delete will fail due to
        foreign key constraints if there is a job associated with the user.

        Args:
            db (SQLClient): a SQLAlchemy database write connection
            email_addr(str): an email address

        """
        #  with nephele2.db_write() as db:
        #      try:
        #          user = db.query(User).filter_by(email_address=email_addr).one()
        #          db.delete(user.user_info)
        #          db.delete(user.user_address)
        #          db.delete(user)
        #      except sqlalchemy.exc.IntegrityError:
        #          db.rollback()
        #          return
        #      except sqlalchemy.orm.exc.NoResultFound:
        #          return
        #
    @staticmethod
    def create_job(user_email, job_id):
        """
        Verifies that the user is valid then creates new row in jobs table.

        Args:
            db (SQLClient): a SQLAlchemy database write connection
            user_email (str): a validly formatted email address. Must exist in
            the database and be marked as confirmed.
            job_id (str): a valid job id

        Raises:
        """
        with nephele2.db_write() as db:
            try:
                user = db.query(User).filter_by(email_address=user_email).one()
                job = Job(job_id=job_id, user=user, status="Initializing")
                db.add(job)
            except sqlalchemy.orm.exc.NoResultFound:
                return

    @staticmethod
    def get_job(job_id):
        """
        Gets a job by it's primary key.
        """
        with nephele2.db_read() as db:
            try:
                return db.query(Job).filter_by(job_id=job_id).one()
            except sqlalchemy.orm.exc.NoResultFound:
                return None

    @staticmethod
    def get_old_unsubmitted_jobs():
        """Looks up jobs that were created > a day ago and have not
        been submitted.

        Args:
            db (SQLClient): a SQLAlchemy database read connection

        Returns:
            list(str): a list of old job IDs.
        """
        yesterday = datetime.now() - timedelta(hours=24)
        with nephele2.db_read() as db:
            jobs = db.query(Job).filter(Job.status == "Initializing",
                                        Job.created < yesterday)
            return [j.job_id for j in jobs.all()]

    @staticmethod
    def get_machine_info_by_job_id(job_id):
        """
        Returns the MachineInfo object associated with the specified job, if
        one exists.

        Args:
            db (SQLClient): a SQLAlchemy database read connection
            job_id (str): a valid job ID (required)

        Returns:
            :obj:`MachinInfo`: a MachineInfo object

        Raises:
        """
        with nephele2.db_read() as db:
            try:
                return db.query(MachineInfo).filter_by(job_id=job_id).one()
            except sqlalchemy.orm.exc.NoResultFound:
                return None

    @staticmethod
    def get_machine_info(machine_id):
        """
        Returns a MachineInfo object.

        Args:
            db (SQLClient): a SQLAlchemy database read connection
            machine_id (str): an AWS machine ID

        Returns:
            :obj:`MachineInfo`: a MachineInfo object

        Raises:
        """
        with nephele2.db_read() as db:
            try:
                return db.query(MachineInfo).filter_by(
                    instance_id=machine_id).one()
            except sqlalchemy.orm.exc.NoResultFound:
                return None

    @staticmethod
    def set_machine_info(job_id, machine_id, ami_id, instance_type):
        """
        Creates a machine_info entry and adds it to the job identified by
        job_id.  If, for some reason, AWS returns the same machine ID twice
        (which is theoretically probable, but highly unlikely) this method will
        fail to create the second entry, since the machine ID is considered a
        primary key and must be unique.

        Args:
            db (SQLClient): a SQLAlchemy database write connection
            job_id (str): a valid job ID
            machine_id (str): an AWS instance identifier
            ami_id (str): an AWS AMI ID
            instance_type (str): the type of EC2 spun up

        """
        machine = MachineInfo(job_id=job_id,
                              instance_id=machine_id,
                              instance_type=instance_type,
                              ami=ami_id)
        with nephele2.db_write() as db:
            db.add(machine)

    @staticmethod
    def get_stack_trace(job_id):
        """
        Returns the stack trace for the job. Value may be NULL.

        Args:
            db (SQLClient): a SQLAlchemy database read connection
            job_id (str): a unique identifier for a job

        Returns:
            str: the message string for the job

        Raises:
        """
        with nephele2.db_read() as db:
            try:
                job = db.query(Job).filter_by(job_id=job_id).one()
                return job.error_msg
            except sqlalchemy.orm.exc.NoResultFound:
                return None

    @staticmethod
    def get_error_msg(job_id):
        """
        Returns the error message for the job. Value may be NULL.

        Args:
            db (SQLClient): a SQLAlchemy database read connection
            job_id (str): a unique identifier for a job

        Returns:
            str: the error message string for the job

        Raises:
        """
        with nephele2.db_read() as db:
            try:
                job = db.query(Job).filter_by(job_id=job_id).one()
                return job.error_msg
            except sqlalchemy.orm.exc.NoResultFound:
                return None

    @staticmethod
    def set_job_status(job_id, status, stack_trace=None, error=None):
        """
        Sets the job status to one of the enum values allowed by the database.

        Args:
            db (SQLClient): a SQLAlchemy database write connection
            job_id (str): a valid job ID
            status (str): must be one of
                ['Pending', 'Pre-Processing', 'Running', 'Failed', 'Succeeded']

        Raises:
            ValueError, RuntimeError
        """
        # TODO this is could be cleaned up.
        valid_statuses = ['Pending', 'Pre-Processing', 'Running', 'Failed',
                          'Succeeded']
        if status not in valid_statuses:
            raise ValueError('Bad Status on job:{}'.format(job_id),
                             status)
        with nephele2.db_write() as db:
            try:
                job = db.query(Job).filter_by(job_id=job_id).one()
                job.status = status
                if status == 'Pending':
                    job.submitted = datetime.now()
                elif status == 'Running':
                    job.started = datetime.now()
                elif status in ['Failed', 'Succeeded']:
                    job.completed = datetime.now()
                if error:
                    if job.error_msg:
                        # if there's already an error don't loose it
                        job.error_msg = job.error_msg + ", " + error
                    else:
                        job.error_msg = error
                if stack_trace:
                    if job.stack_trace:
                        job.stack_trace = job.stack_trace + ', ' + stack_trace
                    else:
                        job.stack_trace = stack_trace
                db.add(job)
            except sqlalchemy.orm.exc.NoResultFound:
                raise RuntimeError('Unable to set job status', job_id)

    @staticmethod
    def get_job_status(job_id):
        """
        Gets the status for the job with the requested ID.

        Args:
            db (SQLClient): a SQLAlchemy database read connection
            job_id (str): a valid job ID

        Returns:
            str: the status of the requested job, one of
            ['Pending','Pre-Processing', 'Running','Failed','Succeeded']

        """
        with nephele2.db_read() as db:
            try:
                job = db.query(Job).filter_by(job_id=job_id).one()
                return job.status
            except sqlalchemy.orm.exc.NoResultFound:
                return None

    @staticmethod
    def set_data_transferred(job_id, is_transferred):
        """
        Sets the transferred flag on a Job.

        Args:
            db (SQLClient): a SQLAlchemy database write connection
            job_id (str): a job id, must exist in the database or the method
            will throw a NepheleError.NepheleRowNotFound exception
            is_transferred (bool): True if data successfully transferred to S3,
            False otherwise

        Raises:
        """
        if is_transferred:
            is_transferred = 1
        else:
            is_transferred = 0

        with nephele2.db_write() as db:
            try:
                job = db.query(Job).filter_by(job_id=job_id).one()
                job.transferred = is_transferred
                db.add(job)
            except sqlalchemy.orm.exc.NoResultFound:
                raise RuntimeError(
                    'Cannot find job, unable to set data transferred', job_id)

    @staticmethod
    def get_data_transferred(job_id):
        """
        Gets the data transfer status for the job with the requested ID.

        Args:
            db (SQLClient): a SQLAlchemy database read connection
            job_id (str): a valid job ID

        Returns:
            bool: transfer status of the requested job

        Raises:
            RuntimeError: if job can't be found
        """
        with nephele2.db_read() as db:
            try:
                job = db.query(Job).filter_by(job_id=job_id).one()
                return job.transferred
            except sqlalchemy.orm.exc.NoResultFound:
                raise RuntimeError(
                    'Cannot find job, unable to get transferred status',
                    job_id)

    @staticmethod
    def set_job_type(job_id, job_name):
        """"
        Sets the job type for the requested job.

        Args:
            db (SQLClient): a SQLAlchemy database write connection
            job_id (str): a valid job ID
            job_name (str): one of the job types in job_type.name

        Raises:
            RuntimeError
        """
        with nephele2.db_write() as db:
            try:
                job = db.query(Job).filter_by(job_id=job_id).one()
                job.job_type = DBUtils.get_job_type_by_name(job_name)
                db.add(job)
            except sqlalchemy.orm.exc.NoResultFound:
                raise RuntimeError(
                    'Unable to find job, failed to set job type', job_id)

    # There are no setters for the JobType table because we shouldn't be doing
    # writes to it, it should be considered a READ ONLY table.
    @staticmethod
    def get_job_type_by_name(name):
        """
        Gets a JobType object by name.

        Args:
            db (SQLClient): a SQLAlchemy database read connection
            name (str): name of the job to run, from job_type.name

        Returns:
            :obj:`JobType`: a JobType object

        Raises:
        """
        with nephele2.db_read() as db:
            try:
                return db.query(JobType).filter(JobType.name == name).one()
            except sqlalchemy.orm.exc.NoResultFound:
                return None

    @staticmethod
    def make_email_subject(job_id):
        """
        Creates a subject line specific to a job.

        Args:
            job_id (str): a unique identifier for a job

        Returns:
            str: the subject line containing the job ID, description, and
            job type

        Raises:
            Exception: raises any exception resulting from the database
            fetch and formatting
        """
        with nephele2.db_read() as db:
            try:
                job = db.query(Job).filter_by(job_id=job_id).one()
                if tfvars.environment in ['dev', 'qa']:
                    return "Nephele Job: {id} ({type} - {desc}) -- {env}"\
                        .format(env=tfvars.environment,
                                id=job_id,
                                type=job.job_type.name,
                                desc=job.description)
                return "Nephele Job: {id} ({type} - {desc})"\
                    .format(id=job_id,
                            type=job.job_type.name,
                            desc=job.description)
            except sqlalchemy.orm.exc.NoResultFound:
                raise RuntimeError(
                    'Unable to find job, failed to make_email_subject', job_id)

    @staticmethod
    def get_job_type_details(job_id):
        """
        Gets details about the job type from the job ID.

        Args:
            db (SQLClient): a SQLAlchemy database read connection
            job_id (str): a valid job ID

        Returns:
            dict: a dict containing the job name and data type::

                {
                  'name': pipe_name,
                  'data_type': data_type,
                  'ami_id': ami_id,
                  'default_instance_type': default_instance_type
                }

        Raises:
        """
        with nephele2.db_read() as db:
            try:
                job = db.query(Job).filter_by(job_id=job_id).one()
                job_details = {}
                job_details['name'] = job.job_type.name
                job_details['data_type'] = job.job_type.data_type
                job_details['ami_id'] = job.job_type.ami_id
                job_details['default_instance_type'] = \
                    job.job_type.default_instance_type
                return job_details
            except sqlalchemy.orm.exc.NoResultFound:
                return None

    @staticmethod
    def get_script_name(job_id):
        """
        Gets the script location for the type of job associated with the
        job_id.

        Args:
            db (SQLClient): a SQLAlchemy database read connection
            job_id (str): a valid job ID

        Returns:
            str: script name

        Examples:
            >>> with nephele2.db_write() as db:
                    a = get_script_name(db,"1")
                    print(a)
            "QIIME/qiime.py"
        """
        with nephele2.db_read() as db:
            try:
                job = db.query(Job).filter_by(job_id=job_id).one()
                return job.job_type.package + '/' + job.job_type.script_name
            except sqlalchemy.orm.exc.NoResultFound:
                return None

    @staticmethod
    def set_job_arguments(job_id, args):
        """
        Adds the arguments to the Job as a JSON string.

        Args:
            job_id (str): a valid job ID
            args (str): a JSON string of job parameters

        Raises:
        """
        args = json.loads(args)
        with nephele2.db_write() as db:
            try:
                job = db.query(Job).filter_by(job_id=job_id).one()
                job.args = args
                db.add(job)
            except sqlalchemy.orm.exc.NoResultFound:
                raise RuntimeError(
                    'Unable to find job, failed to set job arguments',
                    job_id)

    @staticmethod
    def get_job_arguments(job_id):
        """
        Gets the arguments for the requested job.
        These are the arguments that will actually be used for running the job.

        Args:
            job_id (str): a unique identifier for a job

        Returns:
            str: the arguments for the job as a JSON string

        Raises:
        """
        with nephele2.db_read() as db:
            try:
                job = db.query(Job).filter_by(job_id=job_id).one()
                return job.args
            except sqlalchemy.orm.exc.NoResultFound:
                return None

    @staticmethod
    def get_job_id(instance_id):
        """
        Retrieves the job ID mapped to the given machine_id.

        Args:
            instance_id (str): an AWS machine identifier

        Returns:
            str: a job id

        Raises:
        """
        with nephele2.db_read() as db:
            try:
                machine = db.query(MachineInfo).filter_by(
                    instance_id=instance_id).one()
                return machine.job_id
            except sqlalchemy.orm.exc.NoResultFound:
                return None

    @staticmethod
    def set_job_description(job_id, desc):
        """
        Sets the job description field.

        Args:
            job_id (str): a unique job id
            desc (str): string describing the job

        Raises:
        """
        with nephele2.db_write() as db:
            try:
                job = db.query(Job).filter_by(job_id=job_id).one()
                job.description = desc
                db.add(job)
            except sqlalchemy.orm.exc.NoResultFound:
                raise RuntimeError('Unable to set job desc', job_id)

    @staticmethod
    def get_job_description(job_id):
        """
        Returns the job description, may be NULL.

        Args:
            job_id (str): a job id

        Returns:
            str: a string describing the job or NULL

        Raises:
        """

        with nephele2.db_read() as db:
            try:
                job = db.query(Job).filter_by(job_id=job_id).one()
                return job.description
            except sqlalchemy.orm.exc.NoResultFound:
                return None

    #  @staticmethod
    #  def set_checkpoint(job_id, msg):
    #      """
    #      Adds a checkpoint for a job with a timestamp of now.
    #
    #      Args:
    #          job_id (str): a valid job ID
    #          msg (str): the message detailing the checkpoint
    #
    #      Raises:
    #      """
    #      chck_pt = Checkpoints(
    #          job_id=job_id, checkpoint=msg, transition_in=datetime.now())
    #      with db.transaction():
    #          db.add(chck_pt)
