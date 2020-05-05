#!/usr/bin/python3

import os
import random
import string

from nephele2.pipelines.pipebase import PipeBase
import nephele2.pipelines.pipeline_error

class picrust(PipeBase):
    """
    PICRUSt (pronounced pie crust) is a bioinformatics software package designed
    to predict metagenome functional content from marker gene (e.g., 16S rRNA)
    surveys and full genomes.
    """
    gg99_db = "/mnt/EFS/dbs/Greengenes_99/99_otus.fasta"
    otus_dir = "otus_picrust"       # picrust otu output directory
    data_dir = "PICRUSt_data"       # picrust output directory
    l2_param = [
        "summarize_taxa:md_identifier\t\"KEGG_Pathways\"",
        "summarize_taxa:absolute_abundance\tTrue",
        "summarize_taxa:level\t2"]
    l3_param = [
        "summarize_taxa:md_identifier\t\"KEGG_Pathways\"",
        "summarize_taxa:absolute_abundance\tTrue",
        "summarize_taxa:level\t3"]
    gg_param = [
        "pick_otus:enable_rev_strand_match\tTrue",
        "pick_otus:otu_picking_method\tsortmerna",
        "pick_otus:sortmerna_db\t/mnt/EFS/dbs/Greengenes_99/99_otus",
        "pick_otus:refseqs_fp\t/mnt/EFS/dbs/Greengenes_99/99_otus.fasta",
        "alpha_diversity:metrics\tobserved_species,chao1,PD_whole_tree,shannon",
        "make_distance_boxplots:num_permutations\t0",
        "summarize_taxa:level\t2,3,4,5,6,7",
        "summarize_taxa:absolute_abundance\tTrue",
        "make_otu_table:taxonomy\t/mnt/EFS/dbs/Greengenes_99/99_otu_taxonomy.txt",
        "assign_taxonomy:id_to_taxonomy_fp\t/mnt/EFS/dbs/Greengenes_99/99_otu_taxonomy.txt",
        "assign_taxonomy:reference_seqs_fp\t/mnt/EFS/dbs/Greengenes_99/99_otus.fasta",
        "assign_taxonomy:assignment_method\tsortmerna",
        "assign_taxonomy:sortmerna_db\t/mnt/EFS/dbs/Greengenes_99/99_otus"]

    def __init__(self, log_info, in_d = "inputs", out_d = "outputs"):
        self.log_info = log_info
        self.in_dir = in_d
        self.out_dir = out_d
        self.cmds = []

    def param_file(self, p):
        pf = os.path.join(self.out_dir, "%s.txt" % "".join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8)))
        with open(pf, "w") as f:
            f.write("\n".join(p + [""]))
        return(pf)

    def get_cmds(self):
        return(self.cmds)

    def get_output(self):
        self.output.update(self.scan_dir(os.join.path(self.out_dir, picrust.otus_dir)))
        return(self.output)

    def run(self, fasta):
        """
        pick_closed_reference_otus.py
            -i split_lib_out/seqs.fna
            --output_dir=otus_picrust
            --reference_fp=Greengenes_99/99_otus.fasta
            --taxonomy_fp=Greengenes_99/99_otu_taxonomy.txt
            --parameter_fp=picrust_params.txt
            --force

        normalize_by_copy_number.py
            -i otus_picrust/otu_table.biom
            -o PICRUSt_data/normalized_otus.biom

        predict_metagenomes.py
            -i PICRUSt_data/normalized_otus.biom
            -o PICRUSt_data/metagenome_predictions.biom

        categorize_by_function.py
            -i PICRUSt_data/metagenome_predictions.biom
            -c "KEGG_Pathways"
            -l 2
            -o PICRUSt_data/predicted_metagenomes.L2.biom

        categorize_by_function.py
            -i PICRUSt_data/metagenome_predictions.biom
            -c "KEGG_Pathways"
            -l 3
            -o PICRUSt_data/predicted_metagenomes.L3.biom

        summarize_taxa_through_plots.py
            -i PICRUSt_data/predicted_metagenomes.L2.biom
            -o PICRUSt_data/picrust_at_lvl2
            -p qiime_params_picrust2.txt
        """
        # run closed reference otu picking with greengenes 99 database
        cmd = "pick_closed_reference_otus.py -i \"%s\" -o \"%s\" -r \"%s\" -p \"%s\" -f" % (
            fasta, os.path.join(self.out_dir, picrust.otus_dir),
            picrust.gg99_db, self.param_file(picrust.gg_param))
        self.cmds.append(cmd)
        self.log_info(cmd)
        self.exec_cmnd(cmd)

        # normalize an OTU table by marker gene copy number
        cmd = "normalize_by_copy_number.py -i \"%s\" -o \"%s\"" % (
            os.path.join(self.out_dir, picrust.otus_dir, "otu_table.biom"),
            os.path.join(self.out_dir, picrust.data_dir, "normalized_otus.biom"))
        self.cmds.append(cmd)
        self.log_info(cmd)
        self.exec_cmnd(cmd)

        # produces the actual metagenome functional predictions for a given OTU table
        cmd = "predict_metagenomes.py -i \"%s\" -o \"%s\"" % (
            os.path.join(self.out_dir, picrust.data_dir, "normalized_otus.biom"),
            os.path.join(self.out_dir, picrust.data_dir, "metagenome_predictions.biom"))
        self.cmds.append(cmd)
        self.log_info(cmd)
        self.exec_cmnd(cmd)

        # collapse table data to a specified level in a hierarchy
        cmd = "categorize_by_function.py -i \"%s\" -c \"KEGG_Pathways\" -l 2 -o \"%s\"" % (
            os.path.join(self.out_dir, picrust.data_dir, "metagenome_predictions.biom"),
            os.path.join(self.out_dir, picrust.data_dir, "predicted_metagenomes.L2.biom"))
        self.cmds.append(cmd)
        self.log_info(cmd)
        self.exec_cmnd(cmd)

        cmd = "categorize_by_function.py -i \"%s\" -c \"KEGG_Pathways\" -l 3 -o \"%s\"" % (
            os.path.join(self.out_dir, picrust.data_dir, "metagenome_predictions.biom"),
            os.path.join(self.out_dir, picrust.data_dir, "predicted_metagenomes.L3.biom"))
        self.cmds.append(cmd)
        self.log_info(cmd)
        self.exec_cmnd(cmd)

        # summarize taxonomy and generate plots
        cmd = "summarize_taxa_through_plots.py -i \"%s\" -o \"%s\" -p \"%s\"" % (
            os.path.join(self.out_dir, picrust.data_dir, "predicted_metagenomes.L2.biom"),
            os.path.join(self.out_dir, picrust.data_dir, "picrust_at_lvl2"),
            self.param_file(picrust.l2_param))
        self.cmds.append(cmd)
        self.log_info(cmd)
        self.exec_cmnd(cmd)

        cmd = "summarize_taxa_through_plots.py -i \"%s\" -o \"%s\" -p \"%s\"" % (
            os.path.join(self.out_dir, picrust.data_dir, "predicted_metagenomes.L3.biom"),
            os.path.join(self.out_dir, picrust.data_dir, "picrust_at_lvl3"),
            self.param_file(picrust.l3_param))
        self.cmds.append(cmd)
        self.log_info(cmd)
        self.exec_cmnd(cmd)
        return(True)
 