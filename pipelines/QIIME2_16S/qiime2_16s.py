#!/usr/bin/env python3

import csv
import argparse
import traceback
import os

from qiime2 import Artifact, Metadata
from qiime2.plugins.vsearch.methods import join_pairs
from qiime2.plugins.quality_filter.methods import q_score
from qiime2.plugins.demux.visualizers import summarize
from qiime2.plugins.feature_table.visualizers import summarize as feature_table_summarize
from qiime2.plugins.vsearch.methods import dereplicate_sequences
from qiime2.plugins.vsearch.methods import uchime_denovo
from qiime2.plugins.feature_table.methods import filter_features
from qiime2.plugins.feature_table.methods import filter_seqs
from qiime2.plugins.vsearch.methods import cluster_features_de_novo
from qiime2.plugins.feature_table.visualizers import tabulate_seqs
from qiime2.plugins.vsearch.pipelines import cluster_features_open_reference
from qiime2.plugins.vsearch.methods import cluster_features_closed_reference
from qiime2.plugins.deblur.methods import denoise_16S
from qiime2.plugins.fragment_insertion.methods import sepp
from qiime2.plugins.phylogeny.pipelines import align_to_tree_mafft_fasttree
from qiime2.plugins.diversity.visualizers import alpha_rarefaction
from qiime2.plugins.feature_classifier.methods import classify_sklearn
from qiime2.plugins.feature_classifier.methods import classify_consensus_vsearch
from qiime2.plugins.metadata.visualizers import tabulate
from qiime2.plugins.feature_table.methods import filter_samples
from qiime2.plugins.taxa.visualizers import barplot
from biom.cli.metadata_adder import add_metadata
import wurlitzer
import biom

from nephele2.pipelines.pipebase import PipeBase
from nephele2.infra.utils.qiime2_utils import generate_manifest


DEFAULT_NUM_THREADS = 12
EFS_DB = "/mnt/EFS/dbs/"
DATA_TYPE_Q2_MAPPING = {
    "QIIME2_16S_PE": ("SampleData[PairedEndSequencesWithQuality]", "PairedEndFastqManifestPhred33V2"),
    "QIIME2_16S_SE": ("SampleData[SequencesWithQuality]", "SingleEndFastqManifestPhred33V2"),
    # "PJ": ("SampleData[JoinedSequencesWithQuality]", "SingleEndFastqManifestPhred33V2"),
}
OTU_REFERENCE_DCT = {
    "99": EFS_DB + "gg_13_8_otus/rep_set/99_otus.fasta",
    "97": EFS_DB + "gg_13_8_otus/rep_set/97_otus.fasta",
    "85": EFS_DB + "gg_13_8_otus/rep_set/85_otus.fasta",
}
SEPP_REFERENCE_DCT = {
    "silva": EFS_DB + "sepp-refs-silva-128.qza",
    "greengenes": EFS_DB + "sepp-refs-gg-13-8.qza",
}
SKLEARN_OPTIONS_DCT = {
    "silva_v4": EFS_DB + "classifiers/silva-138-99-515-806-nb-classifier.qza",
    "silva_full": EFS_DB + "classifiers/silva-138-99-nb-classifier.qza",
    "greengenes_v4": EFS_DB + "classifiers/gg-13-8-99-515-806-nb-classifier.qza",
    "greengenes_full": EFS_DB + "classifiers/gg-13-8-99-nb-classifier.qza",
}
VSEARCH_REF_READS = EFS_DB + "classifiers/silva-138-99-seqs.qza"
VSEARCH_REF_TAXONOMY = EFS_DB + "classifiers/silva-138-99-tax.qza"

class Qiime2_16S_Error(Exception):
    pass

class QIIME2_16S_Pipeline(PipeBase):
    def __init__(self, args):
        super().__init__(args)


def main(args):
    try:
        exit_status = 0
        pipe = QIIME2_16S_Pipeline(args)
        pipe.log.info("Starting QIIME2 16S pipeline...")

        with open(pipe.outputs_dir + "logfile_debug.txt", "w") as logfile_debug, \
            wurlitzer.pipes(stdout=logfile_debug if args.wurlitzer_stdout == "file" else None, stderr=logfile_debug if args.wurlitzer_stderr == "file" else None):
            # Sanity check
            if args.data_type == "QIIME2_16S_PE" and args.algorithms == "denoising":
                # Front-end will not provide option to choose PE + denoising. However, for sanity check, just put it here
                raise Qiime2_16S_Error("The pipeline does not support Paired End Denoising")

            """
            Determine the number os threads/cpus to use for sepp tree & vsearch classification
            Note: This number is selected based on testing with a specific dataset/AWS EC2 machine/software dependencies
            If you are trying to optimize this number, feel free to do so
            """
            try:
                system_cpu_count = os.cpu_count()
            except Exception as _:
                # If there is anything happens with this cpu_count func, keep continue, no reason to stop the pipeline
                system_cpu_count = DEFAULT_NUM_THREADS + 1

            # If machine has fewer cpus than default, using the max system cpu - 1
            NUM_THREADS = min(DEFAULT_NUM_THREADS, system_cpu_count - 1)

            pipe.log.info("Generating manifest input file...")
            fname = generate_manifest(pipe.map_fp, pipe.inputs_dir)

            # Part 1: Import data
            pipe.log.info("Importing data. Type: {}".format(args.data_type))
            q2_type = DATA_TYPE_Q2_MAPPING[args.data_type]
            input_data = Artifact.import_data(q2_type[0], fname, q2_type[1])
            if args.data_type == "QIIME2_16S_PE":
                joined_pairs = join_pairs(demultiplexed_seqs=input_data)
                input_data = joined_pairs.joined_sequences

            input_data_filtered, _ = q_score(demux=input_data, min_quality=args.min_quality)
            dereplicated_table, dereplicated_sequences = dereplicate_sequences(sequences=input_data_filtered)
            feature_table_visualization = feature_table_summarize(table=dereplicated_table)
            feature_table_visualization.visualization.save(pipe.outputs_dir+'table.qzv')

            # Part 2: Clustering or Denoising
            if args.algorithms == "clustering":
                # Chimera Removal prior to clustering (optional)
                if args.chimeras and args.otu_strategy in ["de_novo", "open"]:
                    pipe.log.info("Chimeras removal running...")
                    chimeras, nonchimeras, chimera_stats = uchime_denovo(sequences=dereplicated_sequences, table=dereplicated_table)
                    filtered_table = filter_features(table=dereplicated_table, metadata=chimeras.view(Metadata), exclude_ids=True)
                    filtered_sequences = filter_seqs(data=dereplicated_sequences, metadata=chimeras.view(Metadata), exclude_ids=True)
                    dereplicated_table = filtered_table.filtered_table
                    dereplicated_sequences = filtered_sequences.filtered_data

                if args.otu_strategy in ["open", "closed"]:
                    # Prepare reference and store in common Nephele reference directory: Import OTU reference sets for clustering
                    pipe.log.info("Prepare OTU reference sets for clustering...")
                    otu_rep_set = Artifact.import_data("FeatureData[Sequence]", OTU_REFERENCE_DCT[args.otu_reference])

                # OTU strategy
                pipe.log.info("Running {} OTU strategy".format(args.otu_strategy))
                if args.otu_strategy == "de_novo":
                    table, sequences = cluster_features_de_novo(sequences=dereplicated_sequences, table=dereplicated_table, perc_identity=args.perc_identity)
                elif args.otu_strategy == "open":
                    table, sequences, _ = cluster_features_open_reference(sequences=dereplicated_sequences, table=dereplicated_table, reference_sequences=otu_rep_set, perc_identity=args.perc_identity)
                elif args.otu_strategy == "closed":
                    table, sequences, _ = cluster_features_closed_reference(sequences=dereplicated_sequences, table=dereplicated_table, reference_sequences=otu_rep_set, perc_identity=args.perc_identity)
                else:
                    # Sanity check
                    raise Qiime2_16S_Error("Unknown OTU strategy")
            elif args.algorithms == "denoising":
                # Denoising with Deblur (only Single End)
                pipe.log.info("Denoising running...")
                table, sequences, _ = denoise_16S(demultiplexed_seqs=input_data_filtered, trim_length=args.trim_length, sample_stats=True)
            else:
                # Sanity check
                raise Qiime2_16S_Error("Unknown {} algorithm".format(str(args.algorithms)))

            pipe.log.info("Exporting biom file...")
            table.save(pipe.outputs_dir + "table.qza")
            table.export_data(output_dir=pipe.outputs_dir) # feature-table.biom
            sequences.save(pipe.outputs_dir + "rep-seqs.qza")
            sequences.export_data(output_dir=pipe.outputs_dir) # dna-sequences.fasta
            sequences_visualization = tabulate_seqs(data=sequences)
            sequences_visualization.visualization.save(pipe.outputs_dir + "rep-seqs.qzv")


            # Part 4: Create a phylogenetic tree to be used in downstream diversity steps
            pipe.log.info("Create a phylogenetic tree: {} with {} threads".format(args.phylogenetic, NUM_THREADS))
            if args.phylogenetic == "sepp":
                reference_database = Artifact.load(SEPP_REFERENCE_DCT[args.sepp_reference])
                tree, placements = sepp(representative_sequences=sequences, reference_database=reference_database, threads=NUM_THREADS, debug=False)
                tree_name = "insertion-tree.qza"
            elif args.phylogenetic == "mafft":
                alignment, masked_alignment, unrooted_tree, tree = align_to_tree_mafft_fasttree(sequences=sequences, n_threads=NUM_THREADS)
                tree_name = "rooted-tree.qza"
                unrooted_tree.export_data(output_dir=pipe.outputs_dir + "unrooted-tree") # unrooted-tree/tree.nwk
                unrooted_tree.save(pipe.outputs_dir + "unrooted-tree/unrooted-tree.qza")
            else:
                # Sanity check
                raise Qiime2_16S_Error("The pipeline does not support {} phylogenetic tree".format(args.phylogenetic))
            pipe.log.info("Saving phylogenetic tree...")
            tree.export_data(output_dir=pipe.outputs_dir + "rooted-tree") # rooted-tree/tree.nwk
            tree.save(pipe.outputs_dir + "rooted-tree/" + tree_name)


            # Part 5: Alpha rarefaction
            pipe.log.info("Running Alpha rarefaction...")
            biom_table = table.view(biom.Table)
            max_depth=int(max(biom_table.sum('sample')))
            pipe.log.info("Alpha rarefaction max_depth: {}".format(max_depth))
            alpha_rarefaction_visualization = alpha_rarefaction(table=table, phylogeny=tree, max_depth=max_depth, metadata=Metadata.load(fname))
            alpha_rarefaction_visualization.visualization.save(pipe.outputs_dir+'alpha-rarefaction.qzv')


            # Part 7: Taxonomy classification
            pipe.log.info("Taxonomy classification: {}...".format(args.taxonomy_methods))
            if args.taxonomy_methods == "sklearn":
                classifier = Artifact.load(SKLEARN_OPTIONS_DCT[args.sklearn_options])
                classification = classify_sklearn(reads=sequences, classifier=classifier, n_jobs=NUM_THREADS)
            elif args.taxonomy_methods == "vsearch":
                ref_reads = Artifact.load(VSEARCH_REF_READS)
                ref_taxonomy = Artifact.load(VSEARCH_REF_TAXONOMY)
                pipe.log.info("Running vsearch classification with {} threads".format(NUM_THREADS))
                classification = classify_consensus_vsearch(query=sequences, reference_reads=ref_reads, reference_taxonomy=ref_taxonomy, threads=NUM_THREADS)
            else:
                # Sanity check
                raise Qiime2_16S_Error("The pipeline does not support {} Taxonomy classification".format(args.taxonomy_methods))
            pipe.log.info("Exporting Taxonomy...")
            classification.classification.save(pipe.outputs_dir + "taxonomy.qza")
            classification.classification.export_data(output_dir=pipe.outputs_dir) # taxonomy.tsv

            # Add taxonomy and metadata to biom file
            pipe.log.info("Adding taxonomy and metadata into biom file...")
            biom_fp = pipe.outputs_dir + "feature-table.biom"
            add_metadata.callback(
                input_fp=biom_fp, output_fp=biom_fp, sample_metadata_fp=args.map_file.name,
                observation_metadata_fp=pipe.outputs_dir+"taxonomy.tsv", sc_separated="taxonomy", sc_pipe_separated=None, int_fields=None,
                float_fields=None, sample_header=None, observation_header="Feature ID,taxonomy", output_as_json=False
            )


            # Part 8: Make plots
            pipe.log.info("Making plots...")
            taxonomy_visualization = tabulate(input=classification.classification.view(Metadata))
            taxonomy_visualization.visualization.save(pipe.outputs_dir+'taxonomy.qzv')

            filtered_table = filter_samples(table=table, min_frequency=args.min_frequency)
            taxa_barplot_visualization = barplot(table=filtered_table.filtered_table, taxonomy=classification.classification, metadata=Metadata.load(fname))
            taxa_barplot_visualization.visualization.save(pipe.outputs_dir+'taxa_barplot.qzv')
    except Exception as exc:
        exit_status = 1
        pipe.log.error("Pipeline Error:")
        pipe.log.error(traceback.format_exc())
        stack = traceback.format_exc()
        msg = str(exc) + '\nPipeline error. Please see logfile.txt for more information.'
        pipe.log_to_db(job_id=args.job_id, stack=stack, msg=msg)
    finally:
        if exit_status == 1:
            pipe.log.error('The pipeline completed with errors.')
        else:
            pipe.log.info('The pipeline completed.')
        pipe.log.info(exit_status)
        exit(exit_status)

    exit(0)

if __name__ == "__main__":
    PARSER = argparse.ArgumentParser()
    MUTEX = PARSER.add_mutually_exclusive_group(required = True)
    MUTEX.add_argument("--job_id", type = str)
    MUTEX.add_argument("--inputs_dir", type = str)
    PARSER.add_argument("--data_type", type = str, default = "QIIME2_16S_PE", choices = ["QIIME2_16S_PE", "QIIME2_16S_SE"], help = "Data type: QIIME2_16S_PE, QIIME2_16S_SE")
    PARSER.add_argument("--map_file", type = argparse.FileType('r'), required = True)
    PARSER.add_argument("--algorithms", type = str, default = "clustering",
        choices = ["clustering", "denoising"], help = "Algorithm: clustering, denoising")

    known_args, _ = PARSER.parse_known_args()
    otu_strategy_default = "closed" if known_args.algorithms == "clustering" else None
    PARSER.add_argument("--otu_strategy", type = str, default=otu_strategy_default,
        choices = ["open", "closed", "de_novo"], help = "OTU picking strategy: open, closed, or de_novo")

    known_args, _ = PARSER.parse_known_args()
    otu_reference_default = "97" if known_args.otu_strategy in ["open", "closed"] else None
    PARSER.add_argument("--otu_reference", type = str, default=otu_reference_default,
        choices = ["99", "97", "85"], help = "OTU reference: 99, 97, or 85")

    known_args, _ = PARSER.parse_known_args()
    chimeras_default = False if known_args.otu_strategy in ["de_novo", "open"] else None
    PARSER.add_argument("--chimeras", action = "store_true", default=chimeras_default)
    PARSER.add_argument("--phylogenetic", type = str, default = "mafft",
        choices = ["sepp", "mafft"], help = "Phylogenetic Tree: sepp, mafft")

    known_args, _ = PARSER.parse_known_args()
    sepp_reference_default = "silva" if known_args.phylogenetic == "sepp" else None
    PARSER.add_argument("--sepp_reference", type = str, default = sepp_reference_default,
        choices = ["silva", "greengenes"], help = "sepp Reference: silva, greengenes")

    PARSER.add_argument("--taxonomy_methods", type = str, default = "sklearn",
        choices = ["sklearn", "vsearch"], help = "Taxonomy classification method: sklearn, vsearch")

    known_args, _ = PARSER.parse_known_args()
    sklearn_options_default = "silva_v4" if known_args.taxonomy_methods == "sklearn" else None
    PARSER.add_argument("--sklearn_options", type = str, default = sklearn_options_default,
        choices = ["silva_v4", "silva_full", "greengenes_v4", "greengenes_full"], help = "sklearn options: silva_v4, silva_full, greengenes_v4, greengenes_full")

    known_args, _ = PARSER.parse_known_args()
    perc_identity_default = 0.97 if known_args.algorithms == "clustering" else None
    PARSER.add_argument('--perc_identity', type=float, default=perc_identity_default)

    known_args, _ = PARSER.parse_known_args()
    trim_length_default = 200 if known_args.algorithms == "denoising" else None
    PARSER.add_argument('--trim_length', type=int, default=trim_length_default)

    # Always there
    PARSER.add_argument('--min_quality', type=int, default=20)
    PARSER.add_argument('--min_frequency', type=int, default=2000)

    PARSER.add_argument("--wurlitzer_stdout", type = str, default = "file", choices = ["file", "std"], help = "choose where to put wurlitzer stdout. Using std to see it from terminal")
    PARSER.add_argument("--wurlitzer_stderr", type = str, default = "file", choices = ["file", "std"], help = "choose where to put wurlitzer stderr. Using std to see it from terminal")

    main(PARSER.parse_args())
