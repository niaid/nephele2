#!/usr/bin/env python3

import os
import csv
import sys
import sh
import biom
import numpy
import argparse
import traceback
from rpy2.robjects.vectors import IntVector
from rpy2.robjects.packages import importr
import rpy2.rinterface
from collections import namedtuple
from nephele2 import config
from nephele2.pipelines import pipeline_error
from nephele2.pipelines.pipebase import PipeBase
from biom.cli.table_summarizer import summarize_table

import nephele2.pipelines.QIIME.pick_otu as pick_otu
import nephele2.pipelines.QIIME.diversity as diversity
import nephele2.pipelines.QIIME.picrust as picrust

class Pipeline(PipeBase):
    def check_file(self, fx, m):
        try:
            self.ensure_file_exists(fx)
        except pipeline_error.NecessaryFileNotFound as fn:
            m = "%s \"%s\" does not exist." % (m, fx)
            self.log.error(m)
            self.log_to_db(job_id = self.job_id, stack = "", msg = m)
            exit(1)

def main(args):
    pipe = Pipeline(args)
    db_tree = {
        "homd": "/mnt/EFS/dbs/homd/HOMD_16S_rRNA_RefSeq_V15.11.tre",
        "sv97": "/mnt/EFS/dbs/SILVA_97/97_otus.tre",
        "sv99": "/mnt/EFS/dbs/SILVA_99/99_otus.tre",
        "gg97": "/mnt/EFS/dbs/Greengenes_97/97_otus.tree",
        "gg99": "/mnt/EFS/dbs/Greengenes_99/99_otus.tree"}
    fastq = {}

    # process the mapping file
    with open(pipe.args.map_file.name) as f:
        reader = csv.DictReader(f, delimiter = "\t")
        for r in reader:
            fastq[r["#SampleID"]] = os.path.join(pipe.inputs_dir, r['ForwardFastqFile'])
            pipe.check_file(fastq[r["#SampleID"]], "Sample file")

    # run otu picking; gg99, gg97, sv99, or sv97
    otu = pick_otu.pick_otu(pipe.log.info, pipe.args.phred_quality, pipe.args.max_bad_run, pipe.args.max_n,
        pipe.args.phred_offset, pipe.inputs_dir, pipe.outputs_dir)
    otu.run_split(fastq)

    otu_tree = ''
    otu_biom = ''
    # specify otu picking strategy and run core diversity
    if pipe.args.otu_strategy == "de_novo":
        otu.run_denovo(pipe.args.ref_db)
        o = otu.get_output()
        otu_tree = o["rep_set.tre"]
        otu_biom = o["otu_table.biom"]
    elif pipe.args.otu_strategy == "closed":
        otu.run_closed(pipe.args.ref_db)
        o = otu.get_output()
        otu_tree = db_tree[pipe.args.ref_db]
        otu_biom = o["otu_table.biom"]
    else:
        otu.run_open(pipe.args.ref_db)
        o = otu.get_output()
        otu_tree = o["rep_set.tre"]
        otu_biom = o["otu_table_mc2_w_tax_no_pynast_failures.biom"]

    # make sure the biom file is not empty
    if not os.path.getsize(otu_biom) > 0:
        pipe.log.info("BIOM file is empty.")
        exit(0)

    pipe.check_file(otu_biom, "OTU table")
    summarize_table.callback(otu_biom, os.path.join(pipe.outputs_dir, "otu_summary_table.txt"), False, False)
    pipe.exec_cmnd("biom convert -i %s -o %s --to-tsv --header-key taxonomy" % (
        otu_biom, os.path.join(pipe.outputs_dir, "otu_table.txt")))
    pipe.exec_cmnd("biom convert -i %s -o %s --to-json --table-type 'OTU table' --process-obs-metadata sc_separated" % (
        os.path.join(pipe.outputs_dir, "otu_table.txt"), os.path.join(pipe.outputs_dir, "otu_table.v1.biom")))
    depth = pipe.get_depth(otu_biom, pipe.args.sampling_depth)

    # no sample has the number of reads greater than 10,000
    if depth < 0:
        pipe.log.info("Visualization pipeline will not be run, as there are not enough samples.")
        exit(0)

    try:
        core = diversity.diversity(pipe.log.info, pipe.args.map_file.name, pipe.inputs_dir, pipe.outputs_dir)
        core.run(depth, otu_biom, otu_tree)
    except Exception:
        m = "QIIME core diversity pipeline failed with unknown errors."
        pipe.log.error(m)
        pipe.log_to_db(job_id = pipe.job_id, stack = traceback.format_exc(), msg = m)
        exit(0)

    # run picrust analysis
    if pipe.args.picrust is True:
        pt = picrust.picrust(pipe.log.info, pipe.inputs_dir, pipe.outputs_dir)
        pipe.check_file(o["seqs.fna"], "Sequence file")
        pt.run(o["seqs.fna"])

    # run additional visualization
    if pipe.args.job_id is not None:
        r_mod_name = 'datavis16s'
        pipe.log.info('Loading R module: {r_mod_name}.'.format(r_mod_name=r_mod_name))
        pipe.load_R_mod(config.PIPELINES_LOC_ON_WRKR + r_mod_name)
    vis = importr("datavis16s")
    vis.trygraphwrapper(otu_biom, pipe.outputs_dir, pipe.args.map_file.name,
        logfilename = os.path.join(pipe.outputs_dir, "logfile.txt"),
        FUN = "allgraphs",
        sampdepth = depth)

if __name__ == "__main__":
    PARSER = argparse.ArgumentParser()
    MUTEX = PARSER.add_mutually_exclusive_group(required = True)
    MUTEX.add_argument("--job_id", type = str)
    MUTEX.add_argument("--inputs_dir", type = str)

    PARSER.add_argument("--data_type", type = str, default = "SE")
    PARSER.add_argument("--map_file", type = argparse.FileType('r'), required = True)
    PARSER.add_argument("--phred_quality", type = int, default = 19)
    PARSER.add_argument("--phred_offset", type = int, default = 33)
    PARSER.add_argument("--max_bad_run", type = int, default = 3)
    PARSER.add_argument("--max_n", type = int, default = 0)
    PARSER.add_argument("--sampling_depth", type=int)
    PARSER.add_argument("--otu_strategy", type = str, default = "open", required = True,
        choices = ["open", "closed", "de_novo"], help = "OTU picking strategy: open, closed, or de_novo")
    PARSER.add_argument('--ref_db', type = str, default = "gg99", required = True,
        choices = ["homd", "sv97", "sv99", "gg97", "gg99"], help = 'reference database; homd, gg99, gg97, sv99, or sv97')
    PARSER.add_argument("--picrust", action = "store_true", default = False)

    main(PARSER.parse_args())
