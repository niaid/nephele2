#!/usr/bin/env python
# -*- coding: utf-8 -*-


import random
import csv
from collections import namedtuple, defaultdict
from nephele2.pipelines import pipeline_error

Sample = namedtuple('Sample', ['SampleID',
                               'ForwardFastqFile',
                               'ReverseFastqFile',
                               'TreatmentGroup',
                               'Description'])


def load_samples_from_map(fname):
    """
    Reads the map file and creates a list of named tuples, one for each
    row.

    Args:
        fname (str): file path for the map file (must be a valid
           paired-end format)

    Returns:
        m2.Sample: a list of namedtuples representing the rows in the map
            file
    """
    samples = list()
    reader = csv.DictReader(open(fname), delimiter='\t')
    for row in reader:
        if row['#SampleID'] == '':
            continue
        sample = Sample(row['#SampleID'],
                        row['ForwardFastqFile'],
                        row['ReverseFastqFile'],
                        row['TreatmentGroup'],
                        row['Description'])
        samples.append(sample)
    if not samples:
        raise RuntimeError('No Samples in Mapping file')
    return samples


def gen_mothur_sample_file(samples, fname):
    """
    Generates the input sample file required by mothur. We home roll ours
    instead of using the mothur make.file() command because we only want to
    use the samples that are listed in the map file, not necessarily all of
    the samples found in the input directory.

    Args:
        samples (list): a list of namedtuples created from the mapping
        file fname (str): the name of the file we're creating

    Raises:
        `pipeline_error.PipelineError`: any of a number of exceptions that
        can be encountered while writing a file or reading an array.
    """
    with open(fname, 'w') as _file:
        for sample in samples:
            line = [sample.SampleID,
                    sample.ForwardFastqFile,
                    sample.ReverseFastqFile]
            print("\t".join(line), file=_file)


#  def link_databases(dname, db_dir, db="sv99"):
#      """
#      Create soft links to the databases in the output directory. We do this
#      to avoid permission issues when running the pcr trim of the database.
#      When mothur runs this pcr command, it tries to create outputs in the
#      location of the *dbs*, not the *output dir*.  The DB dir is owned by
#      root, and it can't (and shouldn't) be able to do this.
#      Uses `PipeBase.create_link`
#
#      Args:
#          dname (str): the path to the outputs directory
#          db_dir (str): the path to the database directory
#          db (str): reference db name from `MothurNeph.DBS`
#
#      Returns:
#          tuple: tuple containing:
#              * ref_db_dest (str): the path for the soft link of the reference
#              * tax_db_dest (str): the path for the soft link of the taxonomy db
#
#      """
#      ref_db = DBS[db]['ref_db']
#      tax_db = DBS[db]['tax_db']
#      ref_db_dest = dname + ref_db
#      tax_db_dest = dname + tax_db
#      PipeBase.create_link(db_dir+ref_db, ref_db_dest)
#      PipeBase.create_link(db_dir+tax_db, tax_db_dest)
#      return ref_db_dest, tax_db_dest
#


def rand_extract(p_value, in_fasta):
    """ Conrad Shyu, 10/17/2017
    Creates a file of randomly extracted sequences from the fasta file.

    Args:
        p_value (float): the cutoff value used to decide if we'll add a
        sequence to the file in_fasta (str): the input fasta file that we'll
        take sequences from

    Returns:
        str: the path of the output file that we created which contains a
        random selection of sequences from the input file

    Raises:
        `pipeline_error.PipelineError`: any exception that occurs while reading
        or writing the files
    """
    try:
        out_fasta = '{in_fasta}_rand_extract_p_at_{p}.fa'\
            .format(in_fasta=in_fasta, p=str(p_value))
        out = open(out_fasta, 'w')
        with open(in_fasta, 'r') as file:
            for entry in file:
                seq = next(file)
                if random.random() < p_value:
                    out.write("%s\n%s\n" % (entry.strip(), seq.strip()))
        return out_fasta
    except Exception as e:
        raise pipeline_error.PipelineError(
            msg='Error occurred when extracting random sequences - function rand_extract.') from e


def calc_maxlength(in_file):
    """updated by Conrad Shyu, 10/17/2017
    Calculates the maximum length of the samples within the fasta file.

    Args:
        in_file (str): the path of the fasta file

    Returns:
        int: the maximum sequence length in the fasta file

    Raises:
        `pipeline_error.PipelineError`: any exception encountered while
        performing this calculation
    """
    try:
        m = defaultdict(int)
        # accumulate the read lengths
        with open(in_file, "r") as file:
            for a in file:
                s = next(file)              # skip the annotation
                m[len(s.strip())] += 1      # construct the histogram
        # return the maximum length in the array
        return int(sorted(m.items(), key=lambda x: x[1], reverse=True)[0][0] + 10)
    except Exception as e:
        raise pipeline_error.PipelineError(
            msg='Error occurred when calculating maximum length of samples within fasta file.') from e


def get_median_start_end(fname):
    """
    Gets the median start and end from the mothur summary file.

    Args:
        fname (str): name of the file with the mothur summary output

    Returns:
        tuple: containing:

            * median start (int): the median start
            * median end (int): the median end

    Raises:
        `pipeline_error.PipelineError`: an error occurred while reading the
        file
    """
    try:
        with open(fname) as f:
            for line in f:
                if line.startswith('Median:'):
                    a = line.split("\t")
                    if len(a) > 3:
                        return a[1], a[2]
                    else:
                        raise pipeline_error.SummarySeqsError()
    except Exception as e:
        raise pipeline_error.PipelineError(msg='Error occurred in getting median start and end.') from e


