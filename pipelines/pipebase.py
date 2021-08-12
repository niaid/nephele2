# -*- coding: utf-8 -*-

"""PipeBase Class
=================
Functionality useful for python pipelines.
Inherit this class to generate python pipelines runnable inside Nephele2.


Pipelines are typically long scripts, and are usually time consuming to test.
They also usually need to be changed, fixed and modified by others. If
idiosyncratic styles or tricks are used this generally increases the cognitive
load on the maintainer, which increases the difficulty of maintenance. Because
of this a prototype pipeline is being put forward, namely demo_pipe/demo.py.

Importantly:
   * The prototype is not final, if there is an aspect of it that you think is
   silly this should be talked over (and changes be made if agreed upon).
   * It is not mandatory. If you don't want to code this way, that's fine.
   Having said that, some work and experience has lead to this, and you should
   be polite and considerate of those that have to work with your code. They
   will more than likely not thank you for doing things in some special unique
   clever way. That is, generally "Strange" is synonymous with "Bad."


"""

from concurrent import futures
import csv
import traceback
import subprocess
import shlex
import logging
import io
import os
import warnings
import numpy
import sh
import biom
from nephele2.pipelines import pipeline_error, utils
from nephele2 import tfvars, config
warnings.simplefilter(action='ignore', category=FutureWarning)


class PipeBase(utils.PipeUtils):
    """Inherit this class to generate python pipelines runnable inside Nephele2.
    """

    DB_DIR = '/mnt/EFS/dbs/'

    def __init__(self, args, inputs_dir=None):
        """
        EITHER Pipeline IS A "nephele" job:
        (Is running from within Nephele => has a job_id)
        OR It's NOT:
        (Is running outside of Nephele infrastructure => has an inputs_dir)

        That is, if a job has a job_id, this implies it's a Nephele Job.
        Else it must have an inputs_dir.

        Makes the outputs directory:
        * if job_id, /path/to/job_id/outputs
        * if inputs_dir, ./outputs in current directory
        * if inputs_dir and outputs_dir, outputs_dir

        gz_col_names is defined in config, it contains a list column header
        names that may contain files which are gzipped.  """
        self._args = args
        if not args.job_id and not args.inputs_dir:
            raise pipeline_error.BadPipeArgsError(
                msg='Must give a namespace with job ID or an inputs dir!')
        if args.job_id:
            self._is_nephele_job = True
            self._job_id = args.job_id
            self._base_dir = config.UPLOAD_PATH + args.job_id + '/'
            self._inputs_dir = self._base_dir + 'inputs/'
            if hasattr(args, 'outputs_dir') and args.outputs_dir is not None:
                if not args.outputs_dir.endswith('/'):
                    args.outputs_dir += '/'
                self._outputs_dir = args.outputs_dir
            else:
                self._outputs_dir = self._base_dir + 'outputs/'
        if args.inputs_dir:
            print('\n--inputs_dir specified => Running outside of Nephele.\n')
            if not hasattr(config, 'SHA_SHORT'):
                try:
                    config.SHA_SHORT = sh.git('--git-dir', os.path.dirname(config.__file__) + '/.git', 'describe').strip()
                except:
                    config.SHA_SHORT = None
            self._is_nephele_job = False

            if not args.inputs_dir.endswith('/'):
                args.inputs_dir += '/'

            self._inputs_dir = args.inputs_dir
            (_, self._job_id) = os.path.split(args.inputs_dir)

            if hasattr(args, 'outputs_dir') and args.outputs_dir is not None:
                if not args.outputs_dir.endswith('/'):
                    args.outputs_dir += '/'
                self._outputs_dir = args.outputs_dir
            else:
                self._outputs_dir = './outputs/'

        if not os.path.isdir(self._inputs_dir):
            raise pipeline_error.NecessaryFileNotFound(
                job_id=args.job_id,
                msg='Unable to find input dir {}'.format(self._inputs_dir))

        if not os.path.isdir(self._outputs_dir):
            sh.mkdir(self._outputs_dir)

        #: logfile :py:mod:`logging` object - full logfile path accessed with ``pipe.log.name``
        self.log = self.setup_logger(self._outputs_dir + 'logfile.txt')

        # adds version number to log
        self.log.info('Nephele, developed by BCBB/OCICB/NIAID/NIH version: '
                      '{}, tag: {}, commit: {}'.format(
                          config.NEPHELE_VERSION,
                          config.NEPHELE_TAG,
                          config.SHA_SHORT))
        self.log.info('Job Description : {}'.format(self.job_description))
        # sorts params alphabetically
        params_log = 'Job parameters\n'
        for key in args.__dict__.keys():
            params_log += '{} : {}\n'.format(key, args.__dict__[key])
        self.log.info(params_log.strip())


        if args.map_file:
            if isinstance(args.map_file, io.TextIOWrapper):
                #: map filename
                self.map_fp = args.map_file.name
            else:
                self.map_fp = args.map_file
            self.log.info('Checking Mapfile for Gzipped inputs.')
            try:
                zipped_files = PipeBase.check_map_file_for_gz(self.map_fp)
                if zipped_files:
                    self.log.info('Gzipped files listed in map file, '
                                  'attempting to rm .gz extension.')
                    self.map_fp = PipeBase.change_zip_fnames_in_map_file(
                        self.map_fp, zipped_files, self.outputs_dir)
                    args.map_file = open(self.map_fp, "r", encoding="utf-8")
                    self.log.info('Done. Attempting file decompression.')
                    utils.PipeUtils.decompr_files(
                        self._inputs_dir, zipped_files)
                    self.log.info('Finished decompression.')
                else:
                    self.cp_if_exists(self.map_fp, self.outputs_dir)
            except FileNotFoundError as fnf:
                self.log.warn(fnf)
            except BaseException:
                self.log.info(
                    'Error while decompressing gzipped files in Mapfile!')
                self.log.info(traceback.format_exc())
                self.log.info('Exiting pipeline')
                exit(1)

    @property
    def job_description(self):
        """job_description
        returns the job_description string entered from the web form
        """
        if self._is_nephele_job:
            try:
                from nephele2.rds.db_utils import DBUtils
                return DBUtils.get_job_description(self.job_id)
            except:
                return 'Unable to find job: {}'.format(self.job_id)


    @property
    def job_id(self):
        """str: The Nephele job ID."""
        return self._job_id

    @property
    def db_dir(self):
        """str: The root of directory containing Nephele Databases. EG
        /mnt/EFS/dbs/"""
        return PipeBase.DB_DIR

    @property
    def base_dir(self):
        """str: The base dir of a job. EG /mnt/EFS/user_data/the_job_id/ """
        return self._base_dir

    @property
    def inputs_dir(self):
        """str: EG /mnt/EFS/user_data/the_job_id/inputs"""
        return self._inputs_dir

    @property
    def outputs_dir(self):
        """str: EG /mnt/EFS/user_data/the_job_id/outputs"""
        return self._outputs_dir

    @property
    def args(self):
        """Namespace: All arguments passed to this pipe via the CLI."""
        return self._args

    @staticmethod
    def get_dir_size(dname):
        """Get size of directory in kilobytes using unix ``du -s``"""
        size = sh.du('-s', dname).split()[0]
        return int(size)

    @staticmethod
    def change_zip_fnames_in_map_file(map_fp, fnames_to_change, out_dname):
        """Modifies map file, rm's gz exts listed in fnames_to_change
        mvs original map file to <fname>.orig
        Args:
            map_fp (str): fully pathed map file
            fnames_to_change (list): ['a.gz', 'b.gz']
        """
        new_fp = out_dname + os.path.basename(map_fp) + '.no_gz'
        with open(map_fp, 'r') as in_f, open(new_fp, 'w') as out_f:
            reader = csv.DictReader(in_f, delimiter='\t')
            writer = csv.DictWriter(
                out_f, fieldnames=reader.fieldnames, delimiter='\t')
            writer.writeheader()
            for row in reader:
                for key in row.keys():
                    if row[key] in fnames_to_change:
                        row[key] = row[key].replace('.gz', '')
                writer.writerow(row)
        return new_fp

    @staticmethod
    def add_to_params_file(fname, lines):
        """adds line to file fname.
        a+ Open for reading and appending (writing at end of file).
        The file is created if it does not exist.
        """
        with open(fname, 'a+') as f_out:
            for line in lines:
                print(line, file=f_out)

    @staticmethod
    def cp_if_exists(f_to_cp, dest):
        """Copies file/folder to dest if dest exists. Does
        nothing if input file doesn't exist.

        Args:
          f_to_cp (str): file or folder to copy
          dest (str): destination
        Raises:
          pipeline_error.NecessaryFileNotFound` """
        if not os.path.lexists(dest):
            msg = 'Destination {0} does not exist.'.format(dest)
            raise pipeline_error.NecessaryFileNotFound(msg=msg)
        if os.path.lexists(f_to_cp):
            sh.cp('-r', f_to_cp, dest)

    @staticmethod
    def ensure_file_exists(fname):
        """Raises a pipeline_error.NecessaryFileNotFound if file can't be
        found."""
        if isinstance(fname, io.TextIOWrapper):
            fname = fname.name
        if not os.path.isfile(fname):
            msg = '{fname} does not exist.\n'.format(fname=fname)
            raise pipeline_error.NecessaryFileNotFound(msg=msg, fname=fname)

    @staticmethod
    def group_files(files, output_dir):
        """attempts to move files into output_dir."""
        try:
            sh.mkdir('-p', output_dir)
            for output in files:
                if os.path.exists(output):
                    sh.mv('-f', output, output_dir)
        except BaseException:
            raise

    @staticmethod
    def setup_logger(log_name):
        """creates a logger with filename log_name and returns it"""
        formatter = logging.Formatter(
            fmt='[%(asctime)s - %(levelname)s] %(message)s')
        handler = logging.FileHandler(log_name)
        handler.setLevel(logging.INFO)
        handler.setFormatter(formatter)
        logger = logging.getLogger(log_name)
        logger.setLevel(logging.INFO)
        logger.addHandler(handler)
        return logger

    @staticmethod
    def exec_cmnd_parallel(cmnds):
        """give the func a list of commands, and it will farm these out to
        cores appropriately.  Will use the number of cpus on the machine.
        if you want to change this, you'd do soimething like:

        with futures.ProcessPoolExecutor(max_workers=7) as pool ...
        """
        with futures.ProcessPoolExecutor() as pool:
            pool.map(PipeBase.exec_cmnd, cmnds)

    @staticmethod
    def exec_cmnd(cmd):
        try:
            args = shlex.split(cmd)
            result = subprocess.run(args, stderr=subprocess.PIPE)
            if result.returncode is not 0:
                raise pipeline_error.UnknownPipeError(msg=result.stderr)
        except BaseException:
            raise

    @staticmethod
    def gen_docker_cmnd(cmd, mnt_pt, cntr_name, db_loc='/mnt/EFS/dbs', db_target='/opt'):
        """Accepts a cmd, ie a (str) command which would work if you were
        running inside an interactive docker session, and a fully pathed input
        dir and output dir (output dir is created if not exists).  This command
        gets prepended with the relevant docker directives, and run.
        For example:

        .. code-block:: python

           cmnd = 'biobakery_workflows  wmgx --input-extension fastq \
           -i /mnt/EFS/user_uploads/bbjob/paired_end/inputs \
           -o /mnt/EFS/user_uploads/bbjob/paired_end/outputs \
           --threads 16'

           PipeBase.docker_exec_cmnd(cmnd, '/mnt/EFS/user_uploads/bbjob', 'biobakery/nephele')

        will generate:

        .. code-block:: bash

           /usr/bin/docker run \
           --mount type=bind,source=/mnt/EFS/dbs,target=/opt \
           --mount type=bind,source=/mnt/EFS/user_uploads/bbjob,target=/mnt/EFS/user_uploads/bbjob \
           biobakery/nephele \
           biobakery_workflows  wmgx --input-extension fastq \
           -i /mnt/EFS/user_uploads/bbjob/paired_end/inputs \
           -o /mnt/EFS/user_uploads/bbjob/paired_end/outputs \
           --threads 16
        """
        if not isinstance(cmd, str):
            raise pipeline_error.NepheleBadArgumentException(
                msg='gen_docker_cmnd() only accepts strings as commands')
        pre = 'run --mount type=bind,source={db_loc},target={db_target} '\
              '--mount type=bind,source={mnt_pt},target={mnt_pt} '\
              '--user www-data '\
              '{cntr_name} '\
              .format(db_loc=db_loc, db_target=db_target, mnt_pt=mnt_pt, cntr_name=cntr_name)
        return pre + cmd

    @staticmethod
    def exec_docker_cmnd(docker_cmnd):
        if not isinstance(docker_cmnd, str):
            raise pipeline_error.NepheleBadArgumentException(
                msg='docker_exec_cmnd() only accepts strings.')

        args = shlex.split(docker_cmnd)
        try:
            docker = sh.docker(args)
        except BaseException:
            raise
        if docker.exit_code is not 0:
            raise pipeline_error.UnknownPipeError(
                msg=docker.stderr + "\n\n" + args)

    @staticmethod
    def load_R_mod(module_fp):
        # '/usr/local/src/nephele2/pipelines/DADA2/dada2nephele'
        try:
            r_cmnd = sh.R('CMD', 'INSTALL', '--no-help', '--use-vanilla', '--preclean', module_fp)
        except BaseException:
            raise
        if r_cmnd.exit_code is not 0:
            raise pipeline_error.UnableToLoadRModuleError(msg=r_cmnd.stderr)

    def scan_dir(self, x):
        """
        scan a directory and generate the list of files under the directory
        written by Conrad Shyu, 3/26/2018
        why is self getting passed here? self is never used.
        This should be a static method or a class method.
        """
        output = {}
        for (r, d, f) in os.walk(x):
            output.update({i: "%s/%s" % (r, i) for i in f})
        return(output)

    @staticmethod
    def count_samples(fname):
        """
        counts the number of samples in a single file (FA or FQ)
        Args:
            fname (file path)
        """
        counter = 0
        with open(fname, 'r') as f_in:
            for line in f_in:
                if line.startswith('>'):
                    counter += 1
        return counter

    @staticmethod
    def create_link(source, target):
        """
        Creates a relative softlink *overwriting any existing path at target*.

        Args:
            source (str): the path of the source file
            target (str): the path for the soft link

        Raises:
            `pipeline_error.NecessaryFileNotFound`: if source file not found or
            insufficient perms at target directory.
        """
        if not os.path.exists(source):
            msg = ("File {0} does not exist, so cannot to make link in "
                   "target {1} directory.").format(
                       os.path.basename(source),
                       os.path.dirname(target))
            raise pipeline_error.NecessaryFileNotFound(msg=msg, fname=source)
        try:
            if (os.path.lexists(target)):
                os.unlink(target)
            os.symlink(os.path.relpath(source, os.path.dirname(target)),
                       target)
        except Exception as pe:
            pe.strerror = "Error linking {0} at target {1}."\
                .format(os.path.basename(source), target)
            msg = str(pe)
            raise pipeline_error.NecessaryFileNotFound(msg=msg,
                                                       fname=target) from pe

    def get_depth(self, biomfile, sampling_depth, mincount=10000):
        """Calculates sampling depth for 16S pipelines
        Poorani Subramanian & Conrad Shyu

        Args:
            biomfile (str):  path of biom file
            sampling_depth (int): desired minimum sample count for downstream
            analysis. can be None
            mincount (int): default minimum sample count used in calculating if
            sampling_depth is None

        Returns:
            sampling_depth (int): If sampling_depth arg is None, set to
            mincount.  If less than 3 sample counts >= sampling_depth, logs
            warning and returns -1.
            """
        if sampling_depth is not None:
            mincount = sampling_depth
        table = biom.load_table(biomfile)
        sample_counts = table.sum(axis='sample')
        good_sc_index = (sample_counts >= mincount)
        good_sc = numpy.extract(good_sc_index, sample_counts)
        if len(good_sc) < 3:
            msg = '{} samples are above the sampling depth of {}, which is '\
                'insufficient.  At least 3 are needed.'.format(
                    str(len(good_sc)), str(mincount))
            self.log.warning(msg)
            return -1
        if sampling_depth is None:
            sampling_depth = numpy.min(good_sc)
        return sampling_depth

    def log_to_db(self, job_id, stack=None, msg=None):
        """ Logs error messages to the db

        Args:
            job_id (`PipeBase.job_id`): run in Nephele
            stack (str): stack/error trace like from :py:func:`traceback.format_exc`
            msg (str): message to send the user in email
        """
        if self._is_nephele_job:
            from nephele2.rds.db_utils import DBUtils
            DBUtils.set_job_status(job_id, 'Failed', stack_trace=stack,
                                   error=msg)
        else:
            print('Emailing:')
            print('Stack : ')
            print(stack)
            print('Message : ')
            print(msg)

    class Conf:
        """PipeBase's conf / constants section."""
        DB_NAME_TO_FILE = {
            'Greengenes_97': 'Greengenes_97.tgz',
            'Greengenes_99': 'Greengenes_99.tgz',
            # Mothur_format_Gg_13_8_99.taxonomy.tgz
            'Greengenes': 'Mothur_format_Gg_13_8_99.taxonomy.tgz',
            'ITS_97': 'ITS_97.tgz',
            'ITS_99': 'ITS_99.tgz',
            'SILVA_97': 'SILVA_97.tgz',
            'SILVA_99': 'SILVA_99.tgz',
            'SILVA': 'SILVA_99.tgz',
            'silva_DB_fname': 'silva.nr_v128.tgz',
        }
        REF_DB_BUCKET_NAME = tfvars.REFDB_BUCKET_ID
