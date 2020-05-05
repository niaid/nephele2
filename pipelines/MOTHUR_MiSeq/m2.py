#!/usr/bin/env python3

"""
Runs the mothur paired-end pipeline using mothur's batch commands.

.. highlight:: bash

"""

import subprocess
import shlex
import os
import re
import argparse
import traceback
import sys
from biom.cli.table_summarizer import summarize_table
import rpy2.rinterface
from rpy2.robjects.packages import importr
from nephele2 import config, tfvars
from nephele2.pipelines import pipeline_error
from nephele2.pipelines import r_pipeline_error
from nephele2.pipelines.pipebase import PipeBase
from nephele2.pipelines import file_utils
from nephele2.pipelines.QIIME import pick_otu
from nephele2.pipelines.QIIME import join_read
from nephele2.pipelines.QIIME import picrust
from nephele2.pipelines.MOTHUR_MiSeq import moth_utils
# useful constants go here
KBYTE = 1024

#: in kilobytes per du man page
MAX_MOTHUR_INPUT = 10000000

HUMAN_READABLE_SIZE = round(MAX_MOTHUR_INPUT/KBYTE**2, 2)

INPUTS_TOO_BIG_MSG = (
    'The inputs for this job are greater than {}GiB. mothur will not be able '
    'to analyse these data due to memory limitations. Please consider another '
    'pipeline.').format(HUMAN_READABLE_SIZE)

#: filename for biom summary
OTU_SUMMARY_TABLE = 'otu_summary_table.txt'

CTG_FNAME = ('combo.trim.contigs.renamed.good.unique.good.filter.unique'
             '.precluster')

#: output files that will be sent to the user
FILES_TO_KEEP = ['logfile.txt',
                 OTU_SUMMARY_TABLE,
                 '{}.denovo.vsearch.pick.pick.count_table'.format(CTG_FNAME),
                 '{}.denovo.vsearch.pick.pick.count.summary'.format(CTG_FNAME),
                 '{}.pick.fasta'.format(CTG_FNAME),
                 '{}.pick.nr_v128.wang.pick.tax.summary'.format(CTG_FNAME),
                 '{}.pick.opti_mcc.0.03.biom'.format(CTG_FNAME),
                 '{}.pick.opti_mcc.0.03.cons.tax.summary'.format(CTG_FNAME),
                 '{}.pick.opti_mcc.0.03.cons.taxonomy'.format(CTG_FNAME),
                 '{}.pick.opti_mcc.list'.format(CTG_FNAME),
                 '{}.pick.opti_mcc.shared'.format(CTG_FNAME),
                 '{}.pick.nr_v128.wang.pick.taxonomy'.format(CTG_FNAME),
                 '{}.pick.phylip.dist'.format(CTG_FNAME),
                 '{}.pick.phylip.tre'.format(CTG_FNAME)]

#: dict mapping reference database name to *ref_db* reference fasta and
#: *tax_db* reference taxonomy.  **Add new reference dbs here.**
DBS = {'homd': {'ref_db': 'HOMD_16S_rRNA_RefSeq_V15.1.aligned.fasta',
                'tax_db': 'HOMD_16S_rRNA_RefSeq_V15.1.mothur.taxonomy'},
       'sv99': {'ref_db': 'silva.nr_v128.align',
                'tax_db': 'silva.nr_v128.tax'}}

#: output directories sent to the user
DIRS_TO_KEEP = ['otus_picrust', 'PICRUSt_data', 'graphs']

UNKNOWN_ERR_MSG = 'An unknown error arose during execution of the pipeline'


def link_samples(samples, inputs_dir, outputs_dir):
    """
    Create soft links for all of the sample files in the output directory.
    Uses `PipeBase.create_link`.

    Args:
        samples (list): a list of namedtuples created from the mapping file
        inputs_dir (str): path to inputs directory
        outputs_dir (str): path to outputs directory

    """
    for sample in samples:
        PipeBase.create_link(inputs_dir+sample.ForwardFastqFile,
                             outputs_dir+sample.ForwardFastqFile)
        PipeBase.create_link(inputs_dir+sample.ReverseFastqFile,
                             outputs_dir+sample.ReverseFastqFile)

class MothurNeph(PipeBase):

    def __init__(self, args):
        """housekeeping stuff that always happens to init pipeline."""
        super().__init__(args)
        #: maximum allowed size of distance matrix in memory - currently 60% of
        #: available mem
        self.mem_size = os.sysconf('SC_PAGE_SIZE') \
            * os.sysconf('SC_PHYS_PAGES') * 0.6
        self.db_name = self.args.ref_db
        self.tax_db_name = DBS[self.db_name]['tax_db']
        self.ref_db_name = DBS[self.db_name]['ref_db']
        self.ref_db = self.outputs_dir + self.ref_db_name
        self.tax_db = self.outputs_dir + self.tax_db_name

    @staticmethod
    def mothurize(dname, cmnd, mothur_log):
        """
        Adds the generic mothur commands that all of the mothur batch scripts
        need to run.

        Args:
            dname (str): path for the outputs directory
            cmnd (str): the batch command we want mothur to run
            mothur_log (str): mothur log filename

        Returns:
            str: the full mothur batch command
        """
        return 'mothur '\
            '"#set.logfile(name={mothur_log}, append=T); '\
            ' set.dir(input={dname}, output={dname}); '\
            ' {cmnd} '\
            ' get.current();" '\
            .format(dname=dname, cmnd=cmnd, mothur_log=mothur_log)

    @staticmethod
    def summary_screen_unique_count(fasta, group, max_ambig, max_length):
        """
        Generates the mothur batch command to run::

            screen.seqs(maxambig=0, maxlength=262, group=x.contigs.groups)
            unique.seqs()
            count.seqs(group=x.contigs.good.groups)

        Args:
            fasta (str): the fasta file
            group (str): the group file
            max_ambig (int): maximum number of ambigous bases allowed,
                sequences with more than this will be removed
            max_length (int): maximum length of any sequence, all sequences
                longer than this number will be removed from the set

        Returns:
            str: the full mothur batch command


        :Outputs:

            * ``x.trim.contigs.good.names`` - from unique.seqs
            * ``x.trim.contigs.good.unique.fasta`` - from unique.seqs
            * ``x.trim.contigs.good.count_table`` - count.seqs
        """
        cmnd = "summary.seqs(fasta={fasta}); "\
               "screen.seqs(maxambig={max_ambig}, maxlength={max_length}, group={group}); "\
               "unique.seqs(); "\
               "count.seqs(group=current); ".format(fasta=fasta,
                                                    group=group,
                                                    max_ambig=max_ambig,
                                                    max_length=max_length)
        return cmnd

    @staticmethod
    def make_contigs(combo_fname, maxee=None):
        """
        Generates the mothur command to run::

            make.contigs(file=x.txt, maxee=maxee)
            rename.seqs(fasta=current, group=current)

        Args:
            combo_fname (str): name of the input file
            maxee (int): maximum number of errors allowed for contigs.
                  if ``None``, arg is not passed.  (I think then mothur uses `default 10000`_ -
                  functionally no filtering on quality/errors).

        .. _default 10000: https://github.com/mothur/mothur/blob/v1.40.5/source/commands/makecontigscommand.cpp#L110

        Returns:
            str: the full mothur batch command

        :Outputs:

            * ``x.trim.contigs.fasta``
            * ``x.trim.contigs.qual``
            * ``x.scrap.contigs.fasta``
            * ``x.scrap.contigs.qual``
            * ``x.contigs.report``
            * ``x.contigs.groups``

        """
        # updated on February 13, 2020, NPHL-1929
        maxee = '' if maxee is None else ', maxee=' + str(maxee)
        cmnd = ("make.contigs(file={combo_fname}{maxee}); "
                "rename.seqs(fasta=current, group=current); ")\
            .format(combo_fname=combo_fname, maxee=maxee)
        return cmnd

    @staticmethod
    def align_summary(seqs, ref):
        """Generates::

             align.seqs(fasta=x.trim.contigs.good.unique.fasta_rand_extract_p_at_0.1.fa,
                reference= silva.seed_v128.align, flip=T)

        Args:
            dname (str): path to the outputs directory
            seqs (str): the path to the fasta file
            ref (str): path to the reference database
            summ_file (str): path to summary file

        Returns:
            str: the full mothur batch command

        :Outputs:

            * ``x.trim.contigs.good.unique.fasta_rand_extract_p_at_0.1.align``
            * ``x.trim.contigs.good.unique.fasta_rand_extract_p_at_0.1.align.report``
        """
        cmnd = 'align.seqs(fasta={seqs}, reference={ref}, flip=T);'\
               'summary.seqs(); '.format(seqs=seqs, ref=ref)
        return cmnd

    @staticmethod
    def pcr_trim(start, end, ref):
        """Generates the command to run::

            pcr.seqs(fasta=silva.nr_v128.align, start={start}, end={end}, keepdots=F);


        Args:
            start (int): a starting position to trim to
            end (int): an ending position to trim from
            ref (str): path to the reference database file that we're trimming

        Returns:
            str: the full mothur batch command


        :Outputs:

            * ``dbname.pcr.align``
        """
        cmd = 'pcr.seqs(fasta={ref}, start={start}, end={end}, keepdots=F); '\
              .format(ref=ref, start=start, end=end)
        return cmd

    @staticmethod
    def align_to_tre_and_biom(align_fp, counts, ref, tax_ref, trimmed_ref, opt,
                              crit, remove_lineage):
        """
        Generates the command to run the main mothur program from align.seqs()
        to summary.tax(). Does taxonomic assignment of pre-clustered sequences.
        Details of the expected input and output can be seen in `mothur_spec`
        and :ref:`README <pipelines.mothur:Analysis Steps and Commands>`.


        Args:
            align_fp (str): path to the input fasta file generated by
                the summary_screen_unique_count steps
            counts (str): path to the input count file generated by
                the summary_screen_unique_count steps
            ref (str): path to the full reference database
            ref (str): path to the taxonomy database
            trimmed_ref (str): path to the trimmed reference database
                generated by pcr.seqs()
            opt (str): value accepted by screen.seqs() to determine where to
                start trimming
            crit (int): defines the percentage used with the optimize parameter
                to determing where to trim
            remove_lineage (str): the taxons to remove

        Returns:
            str: the full mothur batch command
        """
        cmd = 'align.seqs(fasta={align_fp}, reference={trimmed_ref}, flip=T); '\
              'summary.seqs(); '\
              'screen.seqs(optimize={opt}, criteria={crit}, fasta=current, count={counts}); '\
              'filter.seqs(fasta=current, vertical=T, trump=.); '\
              'unique.seqs(fasta=current, count=current); '\
              'pre.cluster(diffs=2, fasta=current, count=current); '\
              'chimera.vsearch(fasta=current, count=current, dereplicate=t); '\
              'remove.seqs(fasta=current, accnos=current); '\
              'summary.seqs(count=current); '\
              'classify.seqs(count=current, reference={ref}, taxonomy={tax_ref}, cutoff=80, probs=f); '\
              'remove.lineage(count=current, taxonomy=current, taxon={remove_lineage}); ' \
              'summary.tax(taxonomy=current, count=current); '\
              .format(align_fp=align_fp,
                      trimmed_ref=trimmed_ref,
                      opt=opt, crit=crit,
                      counts=counts,
                      ref=ref,
                      tax_ref=tax_ref,
                      remove_lineage=remove_lineage)
        return cmd

    # NPHL-1930
    @staticmethod
    def first_split(fasta, counts, taxonomy):
        """
        Run first cluster.split to generate distance matrices for future clustering::

            cluster.split(fasta=combo.trim.contigs.renamed.good.unique.good.filter.unique.precluster.pick.fasta, count=combo.trim.contigs.renamed.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=combo.trim.contigs.renamed.good.unique.good.filter.unique.precluster.pick.nr_v128.wang.pick.taxonomy, splitmethod=classify, taxlevel=4, large=T, cluster=F, cutoff=0.03)

        Args:
          fasta (str): fasta file with sequences to be clustered
          counts (str): count table
          taxonomy (str): preclustered taxonomy

        Returns:
            str: the full mothur batch command

        :Outputs:

            * file: ``*.precluster.pick.file`` - file with list of split dist
            and count files

        """
        cmd = ('cluster.split(fasta={fasta}, '
               'count={counts}, '
               'taxonomy={taxonomy}, '
               'splitmethod=classify, '
               'taxlevel=4, '
               'large=T, '
               'cluster=F, '
               'cutoff=0.03); ')\
            .format(fasta=fasta, counts=counts, taxonomy=taxonomy)
        return cmd

    @staticmethod
    def check_dist_matrix_size(pick_file, mem_size):
        """Check if distance matrix will fit in memory for clustering.

        Args:
          pick_file (str): *file* output of `first_split`
          mem_size (float): memory size limit; set to `MothurNeph.mem_size`

        Raises:
          `pipeline_error.FilesizeTooBig`: Distance matrix is too large.
        """
        msg = ('Distance matrix is too large. mothur will not be able to '
               'cluster OTUs. Visualizations will not be run either.')
        dist_patt = re.compile(r"([\w\./]+?\.dist)\s")
        with open(pick_file) as _file:
            for line in _file:
                for match in re.finditer(dist_patt, line):
                    x = match.group(1).strip()
                    if os.stat(x).st_size > mem_size:
                        raise pipeline_error.FilesizeTooBig(
                            msg=msg,
                            fname=os.path.basename(x),
                            fsize=os.stat(x).st_size/KBYTE**3)

    @staticmethod
    def cluster_split(pick_file, counts, taxonomy, tax_ref):
        """
        Runs second cluster.split which clusters the OTUs.
        Then runs make.shared, classify.otu, and make.biom::

            cluster.split(file=combo.trim.contigs.renamed.good.unique.good.filter.unique.precluster.pick.file, cutoff=0.03, runsensspec=f);
            make.shared(list=current, count=combo.trim.contigs.renamed.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=0.03);
            classify.otu(list=current, count=current, taxonomy=combo.trim.contigs.renamed.good.unique.good.filter.unique.precluster.pick.nr_v128.wang.pick.taxonomy, label=0.03);
            count.groups(count=current);
            make.biom(shared=current, constaxonomy=current, reftaxonomy=silva.nr_v128.tax)


        Args:
          pick_file (str): file output of `first_split`
          counts (str): count table
          taxonomy (str): preclustered taxonomy
          tax_ref (str): reference taxonomy

        Returns:
            str: the full mothur batch command

        :Outputs:

            * list:  ``*.opti_mcc.list`` - clustered OTUs
            * dist: ``*.unique.precluster.pick.dist`` - merged distance matrix
            * shared:  ``*.opti_mcc.shared``
            * taxonomy: ``*.nr_v128.wang.pick.taxonomy``
            * constaxonomy: ``*.opti_mcc.0.03.cons.taxonomy``
            * biom: ``*.opti_mcc.0.03.biom``
            * count: ``*.pick.count_table``

        """
        cmd = 'cluster.split(file={pick_file}, cutoff=0.03, runsensspec=f); '\
              'make.shared(list=current, count={counts}, label=0.03); '\
              'classify.otu(list=current, count=current, taxonomy={taxonomy}, label=0.03); '\
              'count.groups(count=current); '\
              'make.biom(shared=current, constaxonomy=current, reftaxonomy={tax_ref}); '\
              .format(pick_file=pick_file, counts=counts, taxonomy=taxonomy, tax_ref=tax_ref)
        return cmd

    @staticmethod
    def clearcut(fasta):
        """Generates::

            dist.seqs(fasta=combo.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, output=lt)
            clearcut(phylip=combo.trim.contigs.good.unique.good.filter.unique.precluster.pick.phylip.dist)

        Args:
            fasta (str): fasta file from which to build phylogenetic tree

        Returns:
            str: the full mothur batch command

        :Outputs:

            * ``*.precluster.pick.phylip.dist``
            * ``*.precluster.pick.phylip.tre``

        """
        cmd = 'dist.seqs(fasta={fasta}, output=lt); '\
              'clearcut(phylip=current); '\
              .format(fasta=fasta)
        return cmd

    @classmethod
    def check_if_input_too_big(cls, inputs_dir):
        """
        We want to stop processing mothur pipes if inputs are larger than
        `MAX_MOTHUR_INPUT`.  Uses `PipeBase.get_dir_size`

        Args:
           inputs_dir (str): directory to check

        Raises:
          `pipeline_error.FilesizeTooBig`: The inputs for this job are greater
          than {}GiB.
        """
        input_size = PipeBase.get_dir_size(inputs_dir)
        fsize = input_size/KBYTE**2
        if input_size > MAX_MOTHUR_INPUT:
            raise pipeline_error.FilesizeTooBig(msg=INPUTS_TOO_BIG_MSG,
                                                fname=inputs_dir,
                                                fsize=fsize)

    def exec_mothur(self, cmd, dname=None, mothur_log=None):
        """Generate mothur command with `mothurize` and
        execute with :py:func:`subprocess.run`

        Args:
            cmd (str): mothur command to execute
            dname (str): output directory - default is `PipeBase.outputs_dir`
            mothur_log (str): filename for mothur log file - default is
            `PipeBase.log`

        Raises:
            mothurError
        """
        if dname is None:
            dname = self.outputs_dir
        if mothur_log is None:
            mothur_log = self.log.name
        try:
            mothur_cmd = self.mothurize(dname, cmd, mothur_log=mothur_log)
            self.log.info(mothur_cmd)
            cmd_args = shlex.split(mothur_cmd)
            return subprocess.run(cmd_args, check=True)
        except subprocess.CalledProcessError as cpe:
            raise mothurError(msg=str(cpe), cmd=str(cpe.cmd)) from cpe


class mothurError(pipeline_error.PipelineError):
    """Exception thrown when `exec_mothur` command fails"""
    def __init__(self, msg=None, cmd=None, job_id=None):
        super().__init__(msg, job_id=job_id)
        #: mothur command
        self.cmd = cmd


def main(args):
    """
    main function to run mothur pipeline
    """
    exit_status = 0
    try:
        pipe = MothurNeph(args)
        pipe.log.info("Providing access to databases")
        # b
        PipeBase.create_link(pipe.db_dir+pipe.ref_db_name, pipe.ref_db)
        PipeBase.create_link(pipe.db_dir+pipe.tax_db_name, pipe.tax_db)

        # check total size of input directory
        pipe.log.info('Checking input size')
        pipe.check_if_input_too_big(pipe.inputs_dir)

        # link samples, DBS, and make mothur mapping file
        pipe.log.info('Trying to load samples from mapping file.')
        pipe.samples = moth_utils.load_samples_from_map(pipe.map_fp)
        link_samples(pipe.samples, pipe.inputs_dir, pipe.outputs_dir)
        pipe.log.info("Generating mothur input file")
        moth_utils.gen_mothur_sample_file(pipe.samples,
                                          pipe.outputs_dir+'combo.txt')
        cmd = pipe.make_contigs('combo.txt', pipe.args.maxee)
        pipe.exec_mothur(cmd)

        current_files = file_utils.get_current_files(pipe.outputs_dir)
        fasta = current_files['fasta']
        groups = current_files['group']

        pipe.log.info("Successfully generated contigs")

        # run picrust analysis; have to use qiime to join paired ends
        if pipe.args.picrust:
            pe = join_read.join_read(pipe.log.info,
                                     min_overlap=10,
                                     perc_max_diff=25,
                                     in_d=pipe.inputs_dir,
                                     out_d=pipe.outputs_dir)
            pe.run(pipe.args.map_file.name)
            otu = pick_otu.pick_otu(pipe.log.info,
                                    phred_quality_threshold=19,
                                    max_bad_run_length=3,
                                    sequence_max_n=0,
                                    phred_offset=33,
                                    in_d=pipe.inputs_dir,
                                    out_d=pipe.outputs_dir)
            otu.run_split(pe.get_output())
            pt = picrust.picrust(pipe.log.info,
                                 pipe.inputs_dir,
                                 pipe.outputs_dir)
            pt.run((otu.get_output())["seqs.fna"])

        # maxlength is a very sens param. If users spec it, use that number,
        # else generate it.
        if pipe.args.maxlength and pipe.args.maxlength > 0:
            max_len_used = pipe.args.maxlength
            pipe.log.info("Using specified maxlength: {0}".format(
                pipe.args.maxlength))
        else:
            pipe.log.info("Calculating max length")
            max_len_used = moth_utils.calc_maxlength(fasta)
            pipe.log.info("Maxlength calculated at: {}".format(max_len_used))
            pipe.log.info("Using estimated maxlength: {}".format(max_len_used))

        # might be able to simplify this:
        max_ambig = 0               # this is just a setting that Mariam picked
        cmd = pipe.summary_screen_unique_count(fasta,
                                               groups,
                                               max_ambig,
                                               max_len_used)
        pipe.exec_mothur(cmd)
        current_files = file_utils.get_current_files(pipe.outputs_dir)

        good_uniq_seqs = current_files['fasta']
        good_counts = current_files['count']

        pipe.log.info("Calculating median start and end...")
        # randomly extract some fraction of the good_uniq_seqs
        rand_extract_fa = moth_utils.rand_extract(0.1, good_uniq_seqs)
        # mothur chops off the tail of the file and adds .summary, emulate and
        # .summ  swap .ext for .summary
        rand_extract_fa_smmry_fp = \
            rand_extract_fa.rsplit('.', 1)[0] + '.summary.summ'

        # capture the contents of the mothur log, it contains median info
        cmd = pipe.align_summary(rand_extract_fa, pipe.ref_db)
        pipe.exec_mothur(cmd, mothur_log=rand_extract_fa_smmry_fp)

        # append output of align_summary to log file
        file_utils.cat_file(rand_extract_fa_smmry_fp, pipe.log.name)

        median_start, median_end = moth_utils.get_median_start_end(
            rand_extract_fa_smmry_fp)
        pipe.log.info('Median start: ' + median_start)
        pipe.log.info('Median end: ' + median_end)

        cmd = pipe.pcr_trim(median_start, median_end, pipe.ref_db)
        pipe.exec_mothur(cmd)
        current_files = file_utils.get_current_files(pipe.outputs_dir)
        trimmed_ref_db = current_files['fasta']

        pipe.log.info(
            "Pre-clustering and classifying pre-clustered sequences.")
        cmd = pipe.align_to_tre_and_biom(good_uniq_seqs,
                                         good_counts,
                                         pipe.ref_db,
                                         pipe.tax_db,
                                         trimmed_ref_db,
                                         pipe.args.optimize,
                                         pipe.args.criteria,
                                         pipe.args.remove_lineage)
        pipe.exec_mothur(cmd)

        pipe.log.info("Creating the split distance matrices")
        current_files = file_utils.get_current_files(pipe.outputs_dir)
        cmd = pipe.first_split(current_files["fasta"], current_files["count"],
                               current_files["taxonomy"])
        pipe.exec_mothur(cmd)
        current_files = file_utils.get_current_files(pipe.outputs_dir)

        # check the size of distance matrix; NPHL-1930
        pipe.log.info('Checking size of distance matrices')
        pipe.log.info('Max mem size: %f GiB', pipe.mem_size/KBYTE**3)
        pipe.check_dist_matrix_size(current_files['file'], pipe.mem_size)

        pipe.log.info(
            "Clustering and classifing otus, make shared and biom files")
        # needed for dist.seqs in pipe.clearcut below
        final_fasta = current_files['fasta']
        cmd = pipe.cluster_split(current_files['file'], current_files['count'],
                                 current_files['taxonomy'], pipe.tax_db)
        pipe.exec_mothur(cmd)

        # make phylip file
        current_files = file_utils.get_current_files(pipe.outputs_dir)
        biomfile = current_files['biom']
        cmd = pipe.clearcut(final_fasta)
        pipe.exec_mothur(cmd)

        # Summarize biom file
        pipe.log.info('Summarizing biom file to %s.', OTU_SUMMARY_TABLE)
        current_files = file_utils.get_current_files(pipe.outputs_dir)
        pipe.ensure_file_exists(biomfile)
        summarize_table.callback(biomfile,
                                 pipe.outputs_dir + OTU_SUMMARY_TABLE,
                                 False, False)

        # Visualizations
        # get sampling depth
        sampling_depth = pipe.get_depth(biomfile, args.sampling_depth)

        # if sampling depth is positive, try to make graphs
        if sampling_depth > 0:
            try:
                # if running in Nephele reinstall R package
                if pipe.args.job_id:
                    r_mod_name = 'datavis16s'
                    pipe.log.info('Loading R module: {r_mod_name}.'.format(
                        r_mod_name=r_mod_name))
                    pipe.load_R_mod(config.PIPELINES_LOC_ON_WRKR + r_mod_name)

                pipe.log.info('Running data visualization pipeline.')
                # import R library
                datavis16s = importr('datavis16s')
                # run the pipeline
                # note - trygraphwrapper returns an IntVector, take zeroth elt
                # as exit val.
                exit_status = datavis16s.trygraphwrapper(
                    biomfile,
                    pipe.outputs_dir,
                    pipe.map_fp,
                    str(pipe.log.name),
                    FUN='allgraphs',
                    sampdepth=sampling_depth)[0]
            # catch rpy2 exception
            except rpy2.rinterface.RRuntimeError as rpy2_err:
                raise r_pipeline_error.RPipelineError(
                    rpy2_err, args.job_id) from None
        else:
            pipe.log.info('Visualization pipeline will not be run, as there '
                          'are not enough samples.')

    # catch big input data and dist matrix
    except pipeline_error.FilesizeTooBig as ftb:
        if ftb.fname is not None and ftb.fsize is not None:
            pipe.log.warning('%s size: %.6g GiB', ftb.fname, ftb.fsize)
        pipe.log.warning(ftb.msg)
        msg = ftb.msg + (' See this <a href="{}/faq/#collapsemothur">'
                         'FAQ</a>.'.format(tfvars.SERVER_ADDR))
        pipe.log_to_db(job_id=pipe.job_id,
                       stack=traceback.format_exc(),
                       msg=msg)
        exit_status = 1

    # catch errors in linking samples and dbs
    except pipeline_error.NecessaryFileNotFound as fnf:
        pipe.log.error(fnf.msg)
        sys.stderr.write(traceback.format_exc())
        pipe.log_to_db(job_id=pipe.job_id,
                       stack=traceback.format_exc(), msg=fnf.msg)
        exit_status = 1

    # catch errors in mothur commands
    except mothurError as me:
        pipe.log.error('mothur Error:')
        pipe.log.error(me.msg)
        msg = ('mothur encountered an error while running. Please search '
               'logfile.txt for "ERROR" for more information.')
        pipe.log_to_db(job_id=pipe.job_id, stack=me.msg, msg=msg)
        exit_status = 1

    # catch R pipeline error
    except r_pipeline_error.RPipelineError as rerr:
        pipe.log.error('R Pipeline Error:')
        sys.stderr.write(traceback.format_exc())
        pipe.log.error(traceback.format_exc())
        pipe.log.error(rerr)
        pipe.log_to_db(job_id=pipe.job_id, stack=rerr.stack, msg=rerr.msg)
        exit_status = 1

    # catch general but still unknown pipeline error
    except pipeline_error.PipelineError as piperr:
        sys.stderr.write(traceback.format_exc())
        pipe.log.error(traceback.format_exc())
        pipe.log.error(piperr.msg)
        pipe.log_to_db(job_id=pipe.job_id,
                       stack=traceback.format_exc(), msg=piperr.msg)
        exit_status = 1

    # catch big scary error
    except Exception as err:
        sys.stderr.write(traceback.format_exc())
        pipe.log.error(traceback.format_exc())
        pipe.log.error(err)
        pipe.log_to_db(job_id=pipe.job_id, stack=traceback.format_exc(),
                       msg=UNKNOWN_ERR_MSG)
        exit_status = 1

    # on exit, clean up
    finally:
        try:
            if not args.keep:
                pipe.log.info('Removing intermediate files and dirs...')
                files_to_keep = FILES_TO_KEEP + [os.path.basename(pipe.map_fp)]
                file_utils.remove_intermediate_files(pipe.outputs_dir,
                                                     files_to_keep)
                file_utils.remove_intermediate_dirs(pipe.outputs_dir,
                                                    DIRS_TO_KEEP)
            pipe.log.info("Pipeline exiting. %d", exit_status)
        except pipeline_error.PipelineError as p_err:
            sys.stderr.write(traceback.format_exc())
            pipe.log.error(traceback.format_exc())
            pipe.log.error(p_err.msg)
            if exit_status == 0:
                pipe.log_to_db(job_id=pipe.job_id,
                               stack=traceback.format_exc(),
                               msg=p_err.msg)
            exit_status = 1
        except Exception as e:
            sys.stderr.write(traceback.format_exc())
            pipe.log.error(traceback.format_exc())
            pipe.log.error(e)
            pipe.log_to_db(job_id=pipe.job_id,
                           stack=traceback.format_exc(),
                           msg=UNKNOWN_ERR_MSG)
            exit_status = 1

    sys.exit(exit_status)


if __name__ == '__main__':
    PARSER = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    PARSER.add_argument("--job_id", type=str)
    PARSER.add_argument("--inputs_dir",
                        type=str,
                        help="for running outside of Nephele")
    PARSER.add_argument("--outputs_dir",
                        type=str,
                        help="for running outside of Nephele")
    PARSER.add_argument('--map_file',
                        required=True,
                        type=argparse.FileType('r'),
                        help='Name of the mapping file')

    PARSER.add_argument('--data_type', choices=['PE'], default='PE')
    PARSER.add_argument('--maxlength', type=int)
    PARSER.add_argument('--picrust', action='store_true')
    #: choices for optimize parameter in `align_to_tre_and_biom`
    PARSER.add_argument('--optimize',
                        choices=['start-end', 'start', 'end'],
                        default='start-end',
                        help='choices for optimize parameter in \
                        align_to_tre_and_biom')
    PARSER.add_argument('--criteria', type=str, default=90, help="%(type)s")
    PARSER.add_argument('--remove_lineage',
                        default='unknown-Chloroplast-Mitochondria-Eukaryote',
                        type=str, help="%(type)s")
    PARSER.add_argument('--sampling_depth',
                        type=int,
                        help="sampling depth for downstream analysis. "
                        "default %(default)s")
    PARSER.add_argument('--ref_db',
                        type=str,
                        default="sv99",
                        choices=DBS.keys(),
                        help='reference database')
    PARSER.add_argument('--keep',
                        action='store_true',
                        help="keep intermediate files")
    PARSER.add_argument('--maxee',
                        type=int,
                        help="%(type)s. maxee param for make.contigs")
    ARGS = PARSER.parse_args()
    main(ARGS)
