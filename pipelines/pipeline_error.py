from nephele2.NepheleError import NepheleError

class PipelineError(NepheleError):
    """
    Basic exception for errors raised in pipelines.
    All pipelines exceptions should inherit this exception!
    """
    def __init__(self, msg, job_id=None):
        if job_id is None:
            job_id = 'Job ID Unknown'
        super().__init__(job_id=job_id, msg=msg)


class RefDBError(PipelineError):
    pass

class UnknownPipeError(PipelineError):
    pass

class NoSamplesInMappingFile(PipelineError):
    pass


class NecessaryFileNotFound(PipelineError):
    def __init__(self, msg=None, fname=None, job_id=None):
        super().__init__(msg, job_id=job_id)
        self.fname = fname #: filename not found

class InsufficientSamples(PipelineError):
    def __init__(self, job_id=None):
        if job_id is None:
            job_id = 'Job ID Unknown'
        super().__init__('There is an insufficient number of samples. '\
                         'The analysis cannot proceed with the core diversity analysis.', job_id)

class SummarySeqsError(PipelineError):
    def __init__(self, job_id=None):
        if job_id is None:
            job_id = 'Job ID Unknown'
        super().__init__('NEPHELE ERROR: The Output of summary.seqs seems to be wrong', job_id)

class NoDepth(PipelineError):
    def __init__(self, job_id=None):
        if job_id is None:
            job_id = 'Job ID Unknown'
        super().__init__('NEPHELE ERROR: No depth value was found. Cannot compute depth for core diversity.', job_id)

class NoTreatmentGroup(PipelineError):
    def __init__(self, job_id=None):
        if job_id is None:
            job_id = 'Job ID Unknown'
        super().__init__('NEPHELE ERROR: No TreatmentGroup column found in:', job_id)

class AverageLengthTooLow(PipelineError):
    def __init__(self, job_id=None):
        if job_id is None:
            job_id = 'Job ID Unknown'
        super().__init__('NEPHELE ERROR: Average length is too low, unable to continue. Please adjust your parameters to allow for better overlap of PE reads or use the QIIME PE pipeline.')

class UnableToLoadRModuleError(PipelineError):
    pass

class MissingMapFile(PipelineError):
    pass

class NepheleBadArgumentException(PipelineError):
    pass

class BadPipeArgsError(PipelineError):
    pass

class shError(PipelineError):
    """Parse errors raised by `sh.ErrorReturnCode`_ and output of
    commands from sh library to allow addition of usual pipeline_error msg.

    .. _sh.ErrorReturnCode: https://amoffat.github.io/sh/sections/command_class.html#errorreturncode
    """
    def __init__(self, stdout=None, stderr=None, msg=None, job_id=None):
        super().__init__(msg=msg, job_id=job_id)
        self.stdout = stdout  #: sh.ErrorReturnCode stdout output
        self.stderr = stderr  #: sh.ErrorReturnCode stderr output

class FilesizeTooBig(PipelineError):
    """Exception for when files or directories are too big
    for the pipeline to handle.
    """
    def __init__(self, job_id=None, msg=None, fname=None, fsize=None):
        super().__init__(job_id=job_id, msg=msg)
        self.fname = fname #: file/directory that is too big
        self.fsize = fsize #: size of file/directory
