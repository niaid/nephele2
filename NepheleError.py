"""
Nephele specific exception classes.  NepheleError is the base Nephele exception
and it is subclassed by the other exceptions.
"""


class NepheleError(Exception):
    """
    Basic exception for errors raised by nephele2.
    All neph exceptions should inherit this exception!
    """

    VALID_JOB_STATES = ['Pending',
                        'Pre-Processing',
                        'Running',
                        'Failed',
                        'Succeeded']

    def __init__(self, job_id=None, stack_trace=None, msg=None):
        super().__init__(msg, job_id)
        # if job_id is None:
        #     job_id = 'Job ID unknown.'
        # if msg is None:
        #     msg = 'No Details.'
        self.job_id = job_id
        self.msg = msg
        self.stack_trace = stack_trace


class S3_404(NepheleError):
    """
    getting a 404 on a bucket lookup
    System wide.
    """
    pass


class S3Error(NepheleError):
    """ Unknown General S3 problem, eg bucket not exist."""
    pass


class S3JobError(NepheleError):
    """ Unknown Job specific S3 problem, eg prefix not exist."""
    pass


class S3TransferError(NepheleError):
    """sonething went wrong during file transfer"""
    pass


class PresignedUrlError(NepheleError):
    """ isuse with creating presigned URL"""
    pass


class NepheleEmailError(NepheleError):
    """There was a problem getting an email sent.
    """
    pass


class EC2TermError(NepheleError):
    """There was a problem with killing an EC2 machine"""
    pass


class NepheleRowNotFound(NepheleError):
    """There was no row with the requested query value in the database.
    """
    pass
    # def __init__(self, job_id=None, msg=None):
    #     super().__init__(msg, job_id)


class NepheleBadUserAddress(NepheleError):
    """The requested address is marked as bad.
    """
    pass


class NepheleUnregisteredUser(NepheleError):
    """The user was found but has not confirmed the address.
    """
    pass


class NepheleInvalidUser(NepheleError):
    """The requested user is either not allowed to perform
    actions, or does not exist.
    """
    pass


class NepheleInvalidJob(NepheleError):
    """The requested job was not found.
    """
    pass


class NepheleDatabaseWriteError(NepheleError):
    """Something went wrong when attempting to write to the database.
    """
    pass

# TODO: we may just want to use the BadArgumentException instead of this one


class NepheleMissingArgument(NepheleError):
    """A required argument was passes with a None value.
    """
    pass


class NepheleBadArgumentException(NepheleError):
    """A bad argument combination was passed"""
    pass


class NepheleNoDBConnectionException(NepheleError):
    """Unable to reach DB."""
    pass


class NepheleQueueLookupException(NepheleError):
    """something happened when we were trying to get a queue"""
    pass


class NepheleQueueCreationException(NepheleError):
    """Failed to create a queue."""
    pass


class PendingNotAvailableException(NepheleError):
    """ The Pending Queue is a special queue. If it doesn't exist it's bad."""
    pass


class WorkerQueueNotAvailableException(NepheleError):
    """ The worker queue is un reachable for some reason."""
    pass


class NepheleUnknownTypeException(NepheleError):
    """
    Catch all unknown type exception. All I'm doing
    is repackaging any unknown exception as an NepheleError
    """
    pass


class UnableToStartEC2Exception(NepheleError):
    """ something is up with AWS, unable to start EC2"""
    pass


class EC2LimitExceededException(NepheleError):
    """We have reached our limit for running EC2 instances"""
    pass


class BadJobException(NepheleError):
    """Checked job in DB prior to submit, and is malformed."""
    pass


class ScriptException(NepheleError):
    """A script called by Worker._local_exec through an exception"""
    pass


class InvalidJobStatus(NepheleError):
    """
    Job states must be one of config.VALID_JOB_STATES
    """
    pass


class InvalidDataLocation(NepheleError):
    """
    Data location must be either temp or retain.
    """
    pass
# CLASS NepheleWorkerQueueException(NepheleError):
#     """something happened when we were trying to get the worker queue"""
#     pass


class NepheleOSError(NepheleError):
    """
    Something happened while we were looking for files or moving files around.
    """
    pass


class FailedToFindExpectedFileError(NepheleError):
    """If we fail to find an expected file."""
    pass


class NoMoreComputeError(NepheleError):
    """< 0 compute time available."""
    pass


class NoSuchDirError(NepheleError):
    pass


class UnableToRmNFSException(NepheleError):
    pass
