import os
import csv

from errno import ENOENT
from nephele2.pipelines import pipeline_error


class Sample:
    """A Sample must have: id, format (PE|SE), extant fwd fastq file.
    If the format is PE, then it must also have an extant rev fastsq file.
    raises ValueErrors if something is logically like the above
    raises IOError for file not found.

    init method creates Sample.
        * sets : _format, _id, _fwd_fp
        * maybe sets : rev_fp
        * initializes to ``None`` :
            * set if cutadapt is run : cutadpt_fwd_fp, cutadpt_rev_fp
            * set if trimmomatic is run: trimmed_fwd, trimmed_rev
            * set if sample fails pipeline step : failed_step"""
    VALID_SAMPLE_FORMATS = ['PE', 'SE']

    def __init__(self, **kwargs):
        if not (kwargs.get('format') or kwargs.get('id')):
            raise ValueError('Sample must contain an id and format.')

        if not (kwargs.get('fwd_fp') or kwargs.get('rev_fp')):
            raise ValueError('Sample must contain a fwd_fp or an rev_fp.')

        self._format = kwargs['format']
        if self._format not in Sample.VALID_SAMPLE_FORMATS:
            raise ValueError('Sample "format" must be either PE or SE; not "{val}"'\
                             .format(val=kwargs['format']))

        self._id = kwargs['id']
        if kwargs.get('fwd_fp') and kwargs['fwd_fp'] != '':
            if os.path.exists(kwargs['fwd_fp']):
                self._fwd_fp = kwargs['fwd_fp']
            else:
                raise IOError(ENOENT, 'Not a file', self._fwd_fp)

        # All the paired end stuff goes here:
        if self._format == 'PE':
            if kwargs.get('rev_fp') and kwargs['rev_fp'] != '':
                self._rev_fp = kwargs['rev_fp']
            else:
                raise ValueError('Paired End data must have both fwd and rev files.')

            if os.path.exists(self._rev_fp):
                self._rev_fp = kwargs['rev_fp']
            else:
                raise IOError(ENOENT, 'Not a file', self._rev_fp)

        # if cutadapt is run then these are set
        self._cutadpt_fwd_fp = None
        self._cutadpt_rev_fp = None
        self._trimmed_fwd = None
        self._trimmed_rev = None

        # if sample fails in one step in the pipeline, then this is set
        self._failed_step = None

    @property
    def id(self):
        return self._id

    @property
    def format(self):
        return self._format

    @property
    def fwd_fp(self):
        return self._fwd_fp

    @property
    def rev_fp(self):
        if self._format == 'PE':
            return self._rev_fp

    @property
    def cutadpt_fwd_fp(self):
        return self._cutadpt_fwd_fp

    @property
    def cutadpt_rev_fp(self):
        return self._cutadpt_rev_fp

    def _set_cutadpt_fwd_fp(self, fp):
        self._cutadpt_fwd_fp = fp

    def _set_cutadpt_rev_fp(self, fp):
        self._cutadpt_rev_fp = fp

    def load_sample_mfest_data_cutadpt(self, mfest_dname):
        try:
            with open(mfest_dname + 'MANIFEST', 'r') as mfest:
                reader = csv.DictReader(mfest)
                for row in reader:
                    if row['sample-id'] != self.id:
                        continue
                    if row['direction'] == 'forward':
                        self._cutadpt_fwd_fp = mfest_dname + row['filename']
                    if row['direction'] == 'reverse':
                        self._cutadpt_rev_fp = mfest_dname + row['filename']
        except:
            raise

    def set_trimmed_pair(self, fp, ext):
        self._trimmed_fwd = fp + '_1P' + ext
        self._trimmed_rev = fp + '_2P' + ext

    @property
    def trimmed_fwd_fp(self):
        return self._trimmed_fwd

    @property
    def trimmed_rev_fp(self):
        return self._trimmed_rev

    @property
    def failed_step(self):
        """str: step in the pipeline in which the sample failed, otherwise None."""
        return self._failed_step

    @failed_step.setter
    def failed_step(self, step):
        self._failed_step = step

    @staticmethod
    def load_samples(map_fp, in_dir, sample_type):
        """creates a list of Sample objects"""
        samples = list()
        reader = csv.DictReader(open(map_fp), delimiter='\t')
        for row in reader:
            if row['#SampleID'] == '':
                continue
            if sample_type == 'PE':
                smpl = Sample(format=sample_type, id=row['#SampleID'],
                              fwd_fp=os.path.join(in_dir, row['ForwardFastqFile']),
                              rev_fp=os.path.join(in_dir, row['ReverseFastqFile']))
            else:
                smpl = Sample(format=sample_type, id=row['#SampleID'],
                              fwd_fp=os.path.join(in_dir, row['ForwardFastqFile']))
            samples.append(smpl)
        if len(samples) == 0:
            raise pipeline_error.NoSamplesInMappingFile(msg='No Samples in Mapping file.')
        return samples

    @staticmethod
    def filter_failed(sample_list, failed_step=None):
        """filters list of Sample objects based on failure reason

        Args:
          sample_list (list): list of Sample objects
          failed_step (str): value of `failed_step` property to filter on

            * if set to None, will return list with only successful samples
            * if set to 'any', will return list of samples which failed in any step

        Returns:
          list: list of Sample objects with `failed_step` property matching argument `failed_step`.
        """
        if failed_step == 'any':
            return ([s for s in sample_list if s.failed_step is not None])
        return ([s for s in sample_list if s.failed_step == failed_step])
