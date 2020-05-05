#!/usr/bin/python3

import os
import csv
import string

from nephele2.pipelines.pipebase import PipeBase
import nephele2.pipelines.pipeline_error

class join_read(PipeBase):
    """
    The script takes forward and reverse Illumina reads and joins them usign the method chosen,
    fastq-join or SeqPrep.
    """
    def __init__(self, log_info, min_overlap=10, perc_max_diff=25, in_d="inputs", out_d="outputs"):
        self.log_info = log_info
        self.pe_join_method = "fastq-join"
        self.min_overlap = min_overlap
        self.perc_max_diff = perc_max_diff
        self.in_dir = in_d
        self.out_dir = out_d
        self.output = {}
        self.cmds = []

    def get_cmds(self):
        return(self.cmds)

    def get_output(self):
        return(self.output)

    def run(self, map_fp):
        """ run the qiime command """
        with open(map_fp) as csvfile:
            reader = csv.DictReader(csvfile, delimiter='\t')
            for row in reader:
                self.output[row['#SampleID']] = os.path.join(self.out_dir, "join_pair", row['#SampleID'], "fastqjoin.join.fastq")
                cmd = "join_paired_ends.py -f \"%s\" -r \"%s\" -o \"%s\" -m %s -j %d -p %d" % (
                    os.path.join(self.in_dir, row['ForwardFastqFile']),         # forward sequence file
                    os.path.join(self.in_dir, row['ReverseFastqFile']),         # reverse sequence file
                    os.path.join(self.out_dir, "join_pair", row['#SampleID']),  # sample id as output directory
                    self.pe_join_method,                                        # fastq-join C++ program
                    self.min_overlap,                                           # minimum allowed overlap in bases
                    self.perc_max_diff)                                         # maximum allowed % differences within overlap
                self.log_info(cmd)
                self.cmds.append(cmd)
        self.exec_cmnd_parallel(self.cmds)  # run the command
        return(True)
