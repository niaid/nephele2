from nephele2.pipelines import pipeline_error
import re


class RPipelineError(pipeline_error.PipelineError):
    """| Handle errors from R pipelines.
    | **Written by:**
    | Poorani Subramanian

    Parameters:
        rerr (rpy2.rinterface.RRuntimeError) object
        job_id (str): job_id if run inside Nephele

    Attributes:
        msg (str): the error message.
        job_id (str): Unknown if not passed
        stack (str):
            - Error : the error message.
            - Call : function call which actually generated the error.
            - Pipeline Step : function/step within the pipeline where it failed/where the error originated.
            - Pipeline : the top level R function which runs the pipeline.

    Example
    --------
    ::

        ('Error: error error chicken dinner, Call: dada2compute(datadir, outdir, mapfile, ...),
        Pipeline Step: stop("error error chicken dinner"),
        Pipeline: dada2compute', 'Job ID Unknown')

    """
    def __init__(self, rerr=None, job_id=None):
        stack = rerr.args[0]
        msg = re.sub("^Error: \\n ", "", stack)
        msg = re.sub("\\nCall.+", " ", msg)
        stack = re.sub("\\s+", " ", stack)
        stack = re.sub("\\sCall", ", Call", stack)
        msg = re.sub("\\s+", " ", msg)
        super().__init__(msg=msg, job_id=job_id)
        self.stack = stack
