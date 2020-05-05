# -*- coding: utf-8 -*-

"""Pipeline utils class
=======================
putting some funcs in here so they can be made available to other pieces of nephele
"""

import os
import csv
import sh
import errno
from nephele2 import config


class PipeUtils:
    """things that are useful for pipelines but might ALSO be useful for other
    modules"""

    @staticmethod
    def decompr_files(dname, fnames):
        """Runs gunzip on a list of files.

        Args:
            dname (str): The Dirname containing files to gunzip.
            fnames (list): list of filenames in dir.

        Examples:
            Pipe.decompr_files('user/inputs/', ['a.gz', 'b.gz'])
        """
        try:
            for fname in fnames:
                if os.path.exists(dname + fname):
                    sh.gunzip(dname + fname)
                elif os.path.exists(dname + os.path.splitext(fname)[0]):
                    return
                else:
                    raise FileNotFoundError(
                        errno.ENOENT, os.strerror(errno.ENOENT), dname + fname)
        except BaseException:
            raise

    @staticmethod
    def check_map_file_for_gz(map_fp):
        """Examines a map file for entries ending in .gz. Only examines col names in cols_to_check.

        Args:
            map_fp (str): fully pathed map file
            cols_to_check (list): 'col_a', 'col_b' etc."""
        cols_to_check = config.GZ_COL_NAMES
        gz_files = list()
        try:
            with open(map_fp) as csvfile:
                reader = csv.DictReader(csvfile, delimiter='\t')
                for row in reader:
                    for col in cols_to_check:
                        if col in row.keys():
                            if row[col].endswith('.gz'):
                                gz_files.append(row[col])
            return gz_files
        except BaseException:
            raise

    class Error(Exception):
        pass


"""
from: https://codereview.stackexchange.com/questions/32897/efficient-parsing-of-fastq
"""


class Line(str):
    """A line of text with associated filename and line number."""

    def error(self, message):
        """Return an error relating to this line."""
        return Error("{0}({1}): {2}\n{3}"
                     .format(self.filename, self.lineno, message, self))


class Lines():
    """Lines(filename, iterator) wraps 'iterator' so that it yields Line
    objects, with line numbers starting from 1. 'filename' is used in
    error messages.

    """

    def __init__(self, filename, iterator):
        self.filename = filename
        self.lines = enumerate(iterator, start=1)

    def __iter__(self):
        return self

    def __next__(self):
        lineno, s = next(self.lines)
        line = Line(s)
        line.filename = self.filename
        line.lineno = lineno
        return line


def compare_read_names(fname1, fname2):
    found = list()
    try:
        with open(fname1, 'r') as f1_in:
            seq_ids = read_fastq(fname1, f1_in)
            for seq_id in seq_ids:
                found.append(seq_id)
        with open(fname2, 'r') as f2_in:
            seq_ids = read_fastq(fname2, f2_in)
            for seq_id in seq_ids:
                if seq_id not in found:
                    raise ValueError('Bad seq file pair : ' + seq_id)
    except BaseException:
        raise


def read_fastq(filename, iterator):
    """Read FASTQ data from 'iterator' (which may be a file object or any
    other iterator that yields strings) and generate tuples (sequence
    name, sequence data, quality data). 'filename' is used in error
    messages.

    """
    # This implementation follows the FASTQ specification given here:
    # <http://nar.oxfordjournals.org/content/38/6/1767.full>
    import re
    at_seqname_re = re.compile(r'@(.+)$')
    sequence_re = re.compile(r'[!-*,-~]*$')
    plus_seqname_re = re.compile(r'\+(.*)$')
    quality_re = re.compile(r'[!-~]*$')

    lines = Lines(filename, iterator)
    for line in lines:
        # First line of block is @<seqname>.
        m = at_seqname_re.match(line)
        if not m:
            raise line.error("Expected @<seqname> but found:")
        seqname = m.group(1)
        try:
            # One or more lines of sequence data.
            sequence = []
            for line in lines:
                m = sequence_re.match(line)
                if not m:
                    break
                sequence.append(m.group(0))
            if not sequence:
                raise line.error("Expected <sequence> but found:")

            # The line following the sequence data consists of a plus
            # sign and an optional sequence name (if supplied, it must
            # match the sequence name from the start of the block).
            m = plus_seqname_re.match(line)
            if not m:
                raise line.error("Expected +[<seqname>] but found:")
            if m.group(1) not in ['', seqname]:
                raise line.error("Expected +{} but found:".format(seqname))

            # One or more lines of quality data, containing the same
            # number of characters as the sequence data.
            quality = []
            n = sum(map(len, sequence))
            while n > 0:
                line = next(lines)
                m = quality_re.match(line)
                if not m:
                    raise line.error("Expected <quality> but found:")
                n -= len(m.group(0))
                if n < 0:
                    raise line.error("<quality> is longer than <sequence>:")
                quality.append(m.group(0))

#            yield seqname, ''.join(sequence), ''.join(quality)
            # only want the seqname here, everything to first space.
            yield seqname.split(' ')[0]

        except StopIteration:
            raise line.error("End of input before sequence was complete:")
