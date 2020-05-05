#!/usr/bin/python3

import os
import random
import string

from nephele2.pipelines.pipebase import PipeBase
import nephele2.pipelines.pipeline_error

class pick_otu(PipeBase):
    """
    the script join the contigs and perform otu picking, open reference, closed reference, and de novo
    """
    ref_db = {
        "homd": "/mnt/EFS/dbs/homd/HOMD_16S_rRNA_RefSeq_V15.1.p9.fasta",
        "gg99": "/mnt/EFS/dbs/Greengenes_99/99_otus.fasta",
        "gg97": "/mnt/EFS/dbs/Greengenes_97/97_otus.fasta",
        "sv99": "/mnt/EFS/dbs/SILVA_99/otus_16S.fasta",
        "sv97": "/mnt/EFS/dbs/SILVA_97/otus_16S.fasta",
        "its99": "/mnt/EFS/dbs/ITS_99/sh_refs_qiime_ver8_99_02.02.2019.fasta",
        "its97": "/mnt/EFS/dbs/ITS_97/sh_refs_qiime_ver8_97_02.02.2019.fasta"}

    barcode_type = "not-barcoded"   # type of barcode used; golay_12 or not-barcoded
    otu_dir = "otus"
    lib_dir = "split_lib_out"
    join_fasta = "seqs.fna"

    def __init__(self, log_info, phred_quality_threshold=19, max_bad_run_length=3, sequence_max_n=0, phred_offset=33,
        in_d = "inputs", out_d = "outputs"):
        self.log_info = log_info
        self.in_dir = in_d
        self.out_dir = out_d
        self.phred_quality_threshold = phred_quality_threshold
        self.max_bad_run_length = max_bad_run_length
        self.sequence_max_n = sequence_max_n
        self.phred_offset = phred_offset
        self.output = {}
        self.cmds = []

    def param_file(self, p):
        """
        generate parameter files
        written by Conrad Shyu, 3/26/2018
        """
        pf = os.path.join(self.out_dir, "%s.txt" % "".join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8)))
        with open(pf, "w") as f:
            f.write("\n".join(p + [""]))
        return(pf)

    def get_cmds(self):
        return(self.cmds)

    def get_output(self):
        self.output.update(self.scan_dir(os.path.join(self.out_dir, pick_otu.otu_dir)))
        return(self.output)

    def run_split(self, fastq):
        """
        This script performs demultiplexing of Fastq sequence data where barcodes and sequences are contained
        in two separate fastq files (common on Illumina runs).
        """
        cmd = "split_libraries_fastq.py -o \"%s\" -i %s -q %d -n %d -r %d --phred_offset %d --barcode_type=\"%s\" --sample_ids=%s" % (
            os.path.join(self.out_dir, pick_otu.lib_dir),   # output directory
            ",".join(["\"%s\"" % fastq[i] for i in fastq]), # list of input files
            self.phred_quality_threshold,                   # phred quality threshold
            self.sequence_max_n,                            # maximum number of degenerate bases
            self.max_bad_run_length,                        # maximum number of consecutive low quality base calls allowed
            self.phred_offset,                              # ascii offset to use when decoding phred scores (either 33 or 64)
            pick_otu.barcode_type,                          # golay_12 or not-barcoded
            ",".join(["\"%s\"" % i for i in fastq]))        # list of sample ids
        self.cmds.append(cmd)
        self.log_info(cmd)
        self.exec_cmnd(cmd)
        self.output[pick_otu.join_fasta] = os.path.join(self.out_dir, pick_otu.lib_dir, pick_otu.join_fasta)
        return(True if os.path.isfile(self.output[pick_otu.join_fasta]) else False)

    def run_open(self, otus_p = "gg99", user_p = []):
        """
        In an open-reference OTU picking process, reads are clustered against a reference sequence collection
        and any reads which do not hit the reference sequence collection are subsequently clustered de novo.
        pick_open_reference_otus.py is the primary interface for open-reference OTU picking in QIIME, and
        includes taxonomy assignment, sequence alignment, and tree-building steps.
        """
        params = {
            "its99": [
                "pick_otus:enable_rev_strand_match\tTrue",
                "pick_otus:sortmerna_db\t/mnt/EFS/dbs/ITS_99/sh_refs_qiime_ver8_99_02.02.2019",
                "pick_otus.py:pick_otus_reference_seqs_fp\t/mnt/EFS/dbs/ITS_99/sh_refs_qiime_ver8_99_02.02.2019.fasta",
                "assign_taxonomy:id_to_taxonomy_fp\t/mnt/EFS/dbs/ITS_99/sh_taxonomy_qiime_ver8_99_02.02.2019.txt",
                "assign_taxonomy:reference_seqs_fp\t/mnt/EFS/dbs/ITS_99/sh_refs_qiime_ver8_99_02.02.2019.fasta",
                "assign_taxonomy:assignment_method\tsortmerna",
                "assign_taxonomy:sortmerna_db\t/mnt/EFS/dbs/ITS_99/sh_refs_qiime_ver8_99_02.02.2019"],
            "its97": [
                "pick_otus:enable_rev_strand_match\tTrue",
                "pick_otus:sortmerna_db\t/mnt/EFS/dbs/ITS_97/sh_refs_qiime_ver8_97_02.02.2019",
                "pick_otus.py:pick_otus_reference_seqs_fp\t/mnt/EFS/dbs/ITS_97/sh_refs_qiime_ver8_97_02.02.2019.fasta",
                "assign_taxonomy:id_to_taxonomy_fp\t/mnt/EFS/dbs/ITS_97/sh_taxonomy_qiime_ver8_97_02.02.2019.txt",
                "assign_taxonomy:reference_seqs_fp\t/mnt/EFS/dbs/ITS_97/sh_refs_qiime_ver8_97_02.02.2019.fasta",
                "assign_taxonomy:assignment_method\tsortmerna",
                "assign_taxonomy:sortmerna_db\t/mnt/EFS/dbs/ITS_97/sh_refs_qiime_ver8_97_02.02.2019"],
            "homd": [
                "pick_otus:enable_rev_strand_match\tTrue",
                "pick_otus:sortmerna_db\t/mnt/EFS/dbs/homd/HOMD_16S_rRNA_RefSeq_V15",
                "assign_taxonomy:id_to_taxonomy_fp\t/mnt/EFS/dbs/homd/HOMD_16S_rRNA_RefSeq_V15.1.qiime.taxonomy",
                "assign_taxonomy:reference_seqs_fp\t/mnt/EFS/dbs/homd/HOMD_16S_rRNA_RefSeq_V15.1.p9.fasta",
                "assign_taxonomy:assignment_method\tsortmerna",
                "assign_taxonomy:sortmerna_db\t/mnt/EFS/dbs/homd/HOMD_16S_rRNA_RefSeq_V15"],
            "gg99": [
                "pick_otus:enable_rev_strand_match\tTrue",
                "pick_otus:sortmerna_db\t/mnt/EFS/dbs/Greengenes_99/99_otus",
                "assign_taxonomy:id_to_taxonomy_fp\t/mnt/EFS/dbs/Greengenes_99/99_otu_taxonomy.txt",
                "assign_taxonomy:reference_seqs_fp\t/mnt/EFS/dbs/Greengenes_99/99_otus.fasta",
                "assign_taxonomy:assignment_method\tsortmerna",
                "assign_taxonomy:sortmerna_db\t/mnt/EFS/dbs/Greengenes_99/99_otus"],
            "gg97": [
                "pick_otus:enable_rev_strand_match\tTrue",
                "pick_otus:sortmerna_db\t/mnt/EFS/dbs/Greengenes_97/97_otus",
                "assign_taxonomy:id_to_taxonomy_fp\t/mnt/EFS/dbs/Greengenes_97/97_otu_taxonomy.txt",
                "assign_taxonomy:reference_seqs_fp\t/mnt/EFS/dbs/Greengenes_97/97_otus.fasta",
                "assign_taxonomy:assignment_method\tsortmerna",
                "assign_taxonomy:sortmerna_db\t/mnt/EFS/dbs/Greengenes_97/97_otus"],
            "sv99": [
                "pick_otus:enable_rev_strand_match\tTrue",
                "pick_otus:sortmerna_db\t/mnt/EFS/dbs/SILVA_99/otus_16S",
                "assign_taxonomy:id_to_taxonomy_fp\t/mnt/EFS/dbs/SILVA_99/majority_taxonomy_7_levels.txt",
                "assign_taxonomy:reference_seqs_fp\t/mnt/EFS/dbs/SILVA_99/otus_16S.fasta",
                "assign_taxonomy:assignment_method\tsortmerna",
                "assign_taxonomy:sortmerna_db\t/mnt/EFS/dbs/SILVA_99/otus_16S"],
            "sv97": [
                "pick_otus:enable_rev_strand_match\tTrue",
                "pick_otus:sortmerna_db\t/mnt/EFS/dbs/SILVA_97/otus_16S",
                "assign_taxonomy:id_to_taxonomy_fp\t/mnt/EFS/dbs/SILVA_97/majority_taxonomy_7_levels.txt",
                "assign_taxonomy:reference_seqs_fp\t/mnt/EFS/dbs/SILVA_97/otus_16S.fasta",
                "assign_taxonomy:assignment_method\tsortmerna",
                "assign_taxonomy:sortmerna_db\t/mnt/EFS/dbs/SILVA_97/otus_16S"]}

        cmd = "pick_open_reference_otus.py -i \"%s\" -o \"%s\" -m sortmerna_sumaclust -r \"%s\" -p \"%s\" -f %s" % (
            self.output[pick_otu.join_fasta],               # merged sequences in fasta format
            os.path.join(self.out_dir, pick_otu.otu_dir),   # output directory (absolute path)
            pick_otu.ref_db[otus_p],                        # reference database
            self.param_file(params[otus_p] + user_p),       # otus parameters
            "--suppress_align_and_tree" if (otus_p == "its97") or (otus_p == "its99") else "")
        self.cmds.append(cmd)
        self.log_info(cmd)
        self.exec_cmnd(cmd)
        return(True)

    def run_closed(self, otus_p = "gg99", user_p = []):
        """
        In a closed-reference OTU picking process, reads are clustered against a reference sequence collection
        and any reads which do not hit a sequence in the reference sequence collection are excluded from
        downstream analyses. pick_closed_reference_otus.py is the primary interface for closed-reference OTU
        picking in QIIME. If the user provides taxonomic assignments for sequences in the reference database,
        those are assigned to OTUs.
        """
        params = {
            "its99": [
                "pick_otus:enable_rev_strand_match\tTrue",
                "pick_otus:otu_picking_method\tsortmerna",
                "pick_otus:sortmerna_db\t/mnt/EFS/dbs/ITS_99/sh_refs_qiime_ver8_99_02.02.2019",
                "pick_otus.py:pick_otus_reference_seqs_fp\t/mnt/EFS/dbs/ITS_99/sh_refs_qiime_ver8_99_02.02.2019.fasta",
                "assign_taxonomy:id_to_taxonomy_fp\t/mnt/EFS/dbs/ITS_99/sh_taxonomy_qiime_ver8_99_02.02.2019.txt",
                "assign_taxonomy:reference_seqs_fp\t/mnt/EFS/dbs/ITS_99/sh_refs_qiime_ver8_99_02.02.2019.fasta",
                "assign_taxonomy:assignment_method\tsortmerna",
                "assign_taxonomy:sortmerna_db\t/mnt/EFS/dbs/ITS_99/sh_refs_qiime_ver8_99_02.02.2019.fasta"],
            "its97": [
                "pick_otus:enable_rev_strand_match\tTrue",
                "pick_otus:otu_picking_method\tsortmerna",
                "pick_otus:sortmerna_db\t/mnt/EFS/dbs/ITS_97/sh_refs_qiime_ver8_97_02.02.2019",
                "pick_otus.py:pick_otus_reference_seqs_fp\t/mnt/EFS/dbs/ITS_97/sh_refs_qiime_ver8_97_02.02.2019.fasta",
                "assign_taxonomy:id_to_taxonomy_fp\t/mnt/EFS/dbs/ITS_97/sh_taxonomy_qiime_ver8_97_02.02.2019.txt",
                "assign_taxonomy:reference_seqs_fp\t/mnt/EFS/dbs/ITS_97/sh_refs_qiime_ver8_97_02.02.2019.fasta",
                "assign_taxonomy:assignment_method\tsortmerna",
                "assign_taxonomy:sortmerna_db\t/mnt/EFS/dbs/ITS_97/sh_refs_qiime_ver8_97_02.02.2019"],
            "homd": [
                "pick_otus:enable_rev_strand_match\tTrue",
                "pick_otus:otu_picking_method\tsortmerna",
                "pick_otus:sortmerna_db\t/mnt/EFS/dbs/homd/HOMD_16S_rRNA_RefSeq_V15",
                "pick_otus:refseqs_fp\t/mnt/EFS/dbs/homd/HOMD_16S_rRNA_RefSeq_V15.1.p9.fasta",
                "make_otu_table:taxonomy\t/mnt/EFS/dbs/homd/HOMD_16S_rRNA_RefSeq_V15.1.qiime.taxonomy",
                "assign_taxonomy:id_to_taxonomy_fp\t/mnt/EFS/dbs/homd/HOMD_16S_rRNA_RefSeq_V15.1.qiime.taxonomy",
                "assign_taxonomy:reference_seqs_fp\t/mnt/EFS/dbs/homd/HOMD_16S_rRNA_RefSeq_V15.1.p9.fasta",
                "assign_taxonomy:assignment_method\tsortmerna",
                "assign_taxonomy:sortmerna_db\t/mnt/EFS/dbs/homd/HOMD_16S_rRNA_RefSeq_V15"],
            "gg99": [
                "pick_otus:enable_rev_strand_match\tTrue",
                "pick_otus:otu_picking_method\tsortmerna",
                "pick_otus:sortmerna_db\t/mnt/EFS/dbs/Greengenes_99/99_otus",
                "pick_otus:refseqs_fp\t/mnt/EFS/dbs/Greengenes_99/99_otus.fasta",
                "make_otu_table:taxonomy\t/mnt/EFS/dbs/Greengenes_99/99_otu_taxonomy.txt",
                "assign_taxonomy:id_to_taxonomy_fp\t/mnt/EFS/dbs/Greengenes_99/99_otu_taxonomy.txt",
                "assign_taxonomy:reference_seqs_fp\t/mnt/EFS/dbs/Greengenes_99/99_otus.fasta",
                "assign_taxonomy:assignment_method\tsortmerna",
                "assign_taxonomy:sortmerna_db\t/mnt/EFS/dbs/Greengenes_99/99_otus"],
            "gg97": [
                "pick_otus:enable_rev_strand_match\tTrue",
                "pick_otus:otu_picking_method\tsortmerna",
                "pick_otus:sortmerna_db\t/mnt/EFS/dbs/Greengenes_97/97_otus",
                "pick_otus:refseqs_fp\t/mnt/EFS/dbs/Greengenes_97/97_otus.fasta",
                "make_otu_table:taxonomy\t/mnt/EFS/dbs/Greengenes_97/97_otu_taxonomy.txt",
                "assign_taxonomy:id_to_taxonomy_fp\t/mnt/EFS/dbs/Greengenes_97/97_otu_taxonomy.txt",
                "assign_taxonomy:reference_seqs_fp\t/mnt/EFS/dbs/Greengenes_97/97_otus.fasta",
                "assign_taxonomy:assignment_method\tsortmerna",
                "assign_taxonomy:sortmerna_db\t/mnt/EFS/dbs/Greengenes_97/97_otus"],
            "sv99": [
                "pick_otus:enable_rev_strand_match\tTrue",
                "pick_otus:otu_picking_method\tsortmerna",
                "pick_otus:sortmerna_db\t/mnt/EFS/dbs/SILVA_99/otus_16S",
                "pick_otus:refseqs_fp\t/mnt/EFS/dbs/SILVA_99/otus_16S.fasta",
                "make_otu_table:taxonomy\t/mnt/EFS/dbs/SILVA_99/majority_taxonomy_7_levels.txt",
                "assign_taxonomy:id_to_taxonomy_fp\t/mnt/EFS/dbs/SILVA_99/majority_taxonomy_7_levels.txt",
                "assign_taxonomy:reference_seqs_fp\t/mnt/EFS/dbs/SILVA_99/otus_16S.fasta",
                "assign_taxonomy:assignment_method\tsortmerna",
                "assign_taxonomy:sortmerna_db\t/mnt/EFS/dbs/SILVA_99/otus_16S"],
            "sv97": [
                "pick_otus:enable_rev_strand_match\tTrue",
                "pick_otus:otu_picking_method\tsortmerna",
                "pick_otus:sortmerna_db\t/mnt/EFS/dbs/SILVA_97/otus_16S",
                "pick_otus:refseqs_fp\t/mnt/EFS/dbs/SILVA_97/otus_16S.fasta",
                "make_otu_table:taxonomy\t/mnt/EFS/dbs/SILVA_97/majority_taxonomy_7_levels.txt",
                "assign_taxonomy:id_to_taxonomy_fp\t/mnt/EFS/dbs/SILVA_97/majority_taxonomy_7_levels.txt",
                "assign_taxonomy:reference_seqs_fp\t/mnt/EFS/dbs/SILVA_97/otus_16S.fasta",
                "assign_taxonomy:assignment_method\tsortmerna",
                "assign_taxonomy:sortmerna_db\t/mnt/EFS/dbs/SILVA_97/otus_16S"]}

        cmd = "pick_closed_reference_otus.py -i \"%s\" -o \"%s\" -r \"%s\" -p \"%s\" -f" % (
            self.output[pick_otu.join_fasta],               # merged sequences in fasta format
            os.path.join(self.out_dir, pick_otu.otu_dir),   # output directory (absolute path)
            pick_otu.ref_db[otus_p],                        # reference database
            self.param_file(params[otus_p] + user_p))       # otus parameters
        self.cmds.append(cmd)
        self.log_info(cmd)
        self.exec_cmnd(cmd)
        return(True)

    def run_denovo(self, otus_p = "gg99", user_p = []):
        """
        In a de novo OTU picking process, reads are clustered against one another without any external
        reference sequence collection. pick_de_novo_otus.py is the primary interface for de novo OTU picking
        in QIIME, and includes taxonomy assignment, sequence alignment, and tree-building steps. A benefit of
        de novo OTU picking is that all reads are clustered. A drawback is that there is no existing support
        for running this in parallel in QIIME, so it can be too slow to apply to large datasets (e.g., more
        than 10 million reads).
        """
        params = {
            "its99": [
                "pick_otus:enable_rev_strand_match\tTrue",
                "pick_otus:otu_picking_method\tsortmerna",
                "pick_otus:sortmerna_db\t/mnt/EFS/dbs/ITS_99/sh_refs_qiime_ver8_99_02.02.2019",
                "pick_otus.py:pick_otus_reference_seqs_fp\t/mnt/EFS/dbs/ITS_99/sh_refs_qiime_ver8_99_02.02.2019.fasta",
                "assign_taxonomy:id_to_taxonomy_fp\t/mnt/EFS/dbs/ITS_99/sh_taxonomy_qiime_ver8_99_02.02.2019.txt",
                "assign_taxonomy:reference_seqs_fp\t/mnt/EFS/dbs/ITS_99/sh_refs_qiime_ver8_99_02.02.2019.fasta",
                "assign_taxonomy:assignment_method\tsortmerna",
                "assign_taxonomy:sortmerna_db\t/mnt/EFS/dbs/ITS_99/sh_refs_qiime_ver8_99_02.02.2019"],
            "its97": [
                "pick_otus:enable_rev_strand_match\tTrue",
                "pick_otus:otu_picking_method\tsortmerna",
                "pick_otus:sortmerna_db\t/mnt/EFS/dbs/ITS_97/sh_refs_qiime_ver8_97_02.02.2019",
                "pick_otus.py:pick_otus_reference_seqs_fp\t/mnt/EFS/dbs/ITS_97/sh_refs_qiime_ver8_97_02.02.2019.fasta",
                "assign_taxonomy:id_to_taxonomy_fp\t/mnt/EFS/dbs/ITS_97/sh_taxonomy_qiime_ver8_97_02.02.2019.txt",
                "assign_taxonomy:reference_seqs_fp\t/mnt/EFS/dbs/ITS_97/sh_refs_qiime_ver8_97_02.02.2019.fasta",
                "assign_taxonomy:assignment_method\tsortmerna",
                "assign_taxonomy:sortmerna_db\t/mnt/EFS/dbs/ITS_97/sh_refs_qiime_ver8_97_02.02.2019"],
            "homd": [
                "pick_otus:enable_rev_strand_match\tTrue",
                "pick_otus:otu_picking_method\tsumaclust",
                "pick_otus:sortmerna_db\t/mnt/EFS/dbs/homd/HOMD_16S_rRNA_RefSeq_V15",
                "pick_otus:refseqs_fp\t/mnt/EFS/dbs/homd/HOMD_16S_rRNA_RefSeq_V15.1.p9.fasta",
                "assign_taxonomy:id_to_taxonomy_fp\t/mnt/EFS/dbs/homd/HOMD_16S_rRNA_RefSeq_V15.1.qiime.taxonomy",
                "assign_taxonomy:reference_seqs_fp\t/mnt/EFS/dbs/homd/HOMD_16S_rRNA_RefSeq_V15.1.p9.fasta",
                "assign_taxonomy:assignment_method\tsortmerna",
                "assign_taxonomy:sortmerna_db\t/mnt/EFS/dbs/homd/HOMD_16S_rRNA_RefSeq_V15"],
            "gg99": [
                "pick_otus:enable_rev_strand_match\tTrue",
                "pick_otus:otu_picking_method\tsumaclust",
                "pick_otus:sortmerna_db\t/mnt/EFS/dbs/Greengenes_99/99_otus",
                "pick_otus:refseqs_fp\t/mnt/EFS/dbs/Greengenes_99/99_otus.fasta",
                "assign_taxonomy:id_to_taxonomy_fp\t/mnt/EFS/dbs/Greengenes_99/99_otu_taxonomy.txt",
                "assign_taxonomy:reference_seqs_fp\t/mnt/EFS/dbs/Greengenes_99/99_otus.fasta",
                "assign_taxonomy:assignment_method\tsortmerna",
                "assign_taxonomy:sortmerna_db\t/mnt/EFS/dbs/Greengenes_99/99_otus"],
            "gg97": [
                "pick_otus:enable_rev_strand_match\tTrue",
                "pick_otus:otu_picking_method\tsumaclust",
                "pick_otus:sortmerna_db\t/mnt/EFS/dbs/Greengenes_97/97_otus",
                "pick_otus:refseqs_fp\t/mnt/EFS/dbs/Greengenes_97/97_otus.fasta",
                "assign_taxonomy:id_to_taxonomy_fp\t/mnt/EFS/dbs/Greengenes_97/97_otu_taxonomy.txt",
                "assign_taxonomy:reference_seqs_fp\t/mnt/EFS/dbs/Greengenes_97/97_otus.fasta",
                "assign_taxonomy:assignment_method\tsortmerna",
                "assign_taxonomy:sortmerna_db\t/mnt/EFS/dbs/Greengenes_97/97_otus"],
            "sv99": [
                "pick_otus:enable_rev_strand_match\tTrue",
                "pick_otus:otu_picking_method\tsumaclust",
                "pick_otus:sortmerna_db\t/mnt/EFS/dbs/SILVA_99/otus_16S",
                "pick_otus:refseqs_fp\t/mnt/EFS/dbs/SILVA_99/otus_16S.fasta",
                "assign_taxonomy:id_to_taxonomy_fp\t/mnt/EFS/dbs/SILVA_99/majority_taxonomy_7_levels.txt",
                "assign_taxonomy:reference_seqs_fp\t/mnt/EFS/dbs/SILVA_99/otus_16S.fasta",
                "assign_taxonomy:assignment_method\tsortmerna",
                "assign_taxonomy:sortmerna_db\t/mnt/EFS/dbs/SILVA_99/otus_16S"],
            "sv97": [
                "pick_otus:enable_rev_strand_match\tTrue",
                "pick_otus:otu_picking_method\tsumaclust",
                "pick_otus:sortmerna_db\t/mnt/EFS/dbs/SILVA_97/otus_16S",
                "pick_otus:refseqs_fp\t/mnt/EFS/dbs/SILVA_97/otus_16S.fasta",
                "assign_taxonomy:id_to_taxonomy_fp\t/mnt/EFS/dbs/SILVA_97/majority_taxonomy_7_levels.txt",
                "assign_taxonomy:reference_seqs_fp\t/mnt/EFS/dbs/SILVA_97/otus_16S.fasta",
                "assign_taxonomy:assignment_method\tsortmerna",
                "assign_taxonomy:sortmerna_db\t/mnt/EFS/dbs/SILVA_97/otus_16S"]}

        cmd = "pick_de_novo_otus.py -i \"%s\" -o \"%s\" -p \"%s\" -f" % (
            self.output[pick_otu.join_fasta],               # merged sequences in fasta format
            os.path.join(self.out_dir, pick_otu.otu_dir),   # output directory (absolute path)
            self.param_file(params[otus_p] + user_p))       # otus parameters
        self.cmds.append(cmd)
        self.log_info(cmd)
        self.exec_cmnd(cmd)
        return(True)
