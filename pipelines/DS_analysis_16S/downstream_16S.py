#!/usr/bin/env python3

"""

Pipeline for performing downstream analysis on a biom file.

| NEEDS TO RUN:
| PATH="/usr/local/bin/miniconda3/bin:$PATH"
| source activate qiime2-2018.6
| PYTHONPATH=$PYTHONPATH:/usr/local/lib/python3.5/dist-packages/

- You can't call this script directly due to qiime2 dependency; instead use downstream_16S.sh_ which contains ^^^
- but the arguments and usage are otherwise the same.

.. _downstream_16S.sh: https://github.com/niaid/nephele2/blob/master/pipelines/DS_analysis_16S/downstream_16S.sh

"""

import argparse
import traceback
import sys
import os

from nephele2.pipelines.pipebase import PipeBase

import qiime2.plugin
import qiime2.sdk
from qiime2.plugins.feature_table.visualizers import summarize
from qiime2.metadata import metadata
from qiime2.plugins.diversity.pipelines import core_metrics
from qiime2.plugins.diversity.visualizers import alpha_group_significance

import biom
from q2_types.feature_data._transformer import _biom_to_tsv_taxonomy_format
from qiime2.plugins.feature_table.methods import filter_samples
from qiime2.plugins.taxa.visualizers import barplot


def main(args):
    Q2_TYPE = 'FeatureTable[Frequency]'  # AKA type
    VIEW_TYPE = 'BIOMV100Format'  # AKA input_type AKA source_type
    try:
        exit_status = 0
        pipe = DSAnalysis16S(args)

        # capture stderr and stdout to log file
        origout = sys.stdout
        origerror = sys.stderr
        sys.stderr = open(pipe.log.name, 'a')
        sys.stdout = sys.stderr
        ####
        if (biom.util.is_hdf5_file(args.biom_fp.name)):
            VIEW_TYPE='BIOMV210Format'
        pipe.log.info(
            'Map file : {}, Biom file: {}, Input type : {}'.format(
                args.map_file.name,
                args.biom_fp.name,
                VIEW_TYPE))

        """
        qiime tools import \
        --input-path combo.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.0.03.biom \
        --type 'FeatureTable[Frequency]' \
        --source-format BIOMV100Format \
        --output-path table-qiime2.qza
        """
        pipe.log.info('Trying to create artifact.')

        artifact = qiime2.sdk.Artifact.import_data(
            Q2_TYPE, args.biom_fp.name, view_type=VIEW_TYPE)
        # pipe.log.info('Done.\nSaving artifact.')
        # artifact.save(pipe.outputs_dir + '/Q2_table')
        pipe.log.info('Done.\nLoading Metadata from map file.')
        sys.stdout.flush()

        """
        qiime tools import \
        --input-path combo.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.0.03.biom \
        --type 'FeatureTable[Frequency]' \
        --source-format BIOMV100Format \
        --output-path table-qiime2.qza
        """
        md = metadata.Metadata.load(args.map_file.name)
        pipe.log.info('Done.\n')
        pipe.log.info(
            'Filter out samples from biom file which are not in the mapping file.\n'
            'artifact = filter_samples(table=artifact, metadata=md).filtered_table')
        artifact = filter_samples(table=artifact, metadata=md).filtered_table
        pipe.log.info('Done.\n')

        pipe.log.info('Creating Summary')
        """
        qiime feature-table summarize \
        --i-table table-qiime2.qza \
        --o-visualization table-qiime2.qzv \
        --m-sample-metadata-file Mapping_file_corrected.txt.nogz
        """
        summary = summarize(table=artifact, sample_metadata=md)
        summary.visualization.save(pipe.outputs_dir + '/summary')
        sys.stdout.flush()
        pipe.log.info('Summary done.\n')


        """
        BARPLOT
        -------
        import taxonomy from biom file
        """
        pipe.log.info('Starting barplot steps.')
        pipe.log.info('Importing taxonomy.\n')
        if (VIEW_TYPE == 'BIOMV100Format'):
            pipe.log.info('biomtable = biom.load_table(args.biom_fp.name)')
            biomtable = biom.load_table(args.biom_fp.name)
            pipe.log.info('tf = _biom_to_tsv_taxonomy_format(biomtable)')
            tf = _biom_to_tsv_taxonomy_format(biomtable)
            TAX_VIEW_TYPE = None
        else:
            tf = args.biom_fp.name
            TAX_VIEW_TYPE = VIEW_TYPE
        pipe.log.info(
            'tax_tab = qiime2.sdk.Artifact.import_data(type="FeatureData[Taxonomy]",'
            'view=tf, view_type=TAX_VIEW_TYPE)')
        tax_tab = qiime2.sdk.Artifact.import_data(type="FeatureData[Taxonomy]",
                                                  view=tf, view_type=TAX_VIEW_TYPE)
        os.remove(str(tf))

        """
        filter feature table by sampling depth
        """
        pipe.log.info(
            'Filter samples out of table based on sampling depth.\n'
            'filtartifact = filter_samples(table=artifact, min_frequency=args.sampling_depth)')
        filtartifact = filter_samples(
            table=artifact, min_frequency=args.sampling_depth)

        """
        make barplot
        """
        pipe.log.info(
            'Creating barplot.\n'
            'bp = barplot(table=filtartifact.filtered_table, taxonomy=tax_tab, metadata=md)')
        bp = barplot(
            table=filtartifact.filtered_table,
            taxonomy=tax_tab,
            metadata=md)
        pipe.log.info(
            'Saving barplot visualization to barplot.qzv.\n'
            'bp.visualization.save(pipe.outputs_dir + "barplot.qzv")')
        bp.visualization.save(pipe.outputs_dir + "barplot.qzv")
        pipe.log.info('Barplot done.\n')


        """
        qiime diversity core-metrics
        --i-table table-qiime2.qza
        --p-sampling-depth 20000
        --m-metadata-file Mapping_file_corrected.txt.nogz
        --output-dir qiime2-core-metrics-results
        """
        pipe.log.info('Creating Core Metrics.')
        cm = core_metrics(artifact, args.sampling_depth, md)
        sys.stdout.flush()
        pipe.log.info('Saving.')
        for field in cm._fields:
            cm.__getattribute__(field).save(pipe.outputs_dir + '/'+field)
        for field in ['bray_curtis_distance_matrix', 'bray_curtis_pcoa_results',
                      'evenness_vector', 'jaccard_distance_matrix',
                      'jaccard_pcoa_results', 'observed_otus_vector',
                      'rarefied_table', 'shannon_vector']:
            cm.__getattribute__(field).export_data(
                pipe.outputs_dir + '/'+field)
        pipe.log.info('Core metrics done.\n')

        if args.alpha_group_sig:
            """
            qiime diversity alpha-group-significance
            --i-alpha-diversity qiime2-core-metrics-results/evenness_vector.qza
            --m-metadata-file Mapping_file_corrected.txt.nogz
            --o-visualization qiime2-core-metrics-results/evenness-group-significance.qzv
            """
            pipe.log.info('Creating Alpha Group Significance for evenness index.\n'
                          'ags_even = alpha_group_significance(cm.evenness_vector, md)')
            ags_even = alpha_group_significance(cm.evenness_vector, md)
            sys.stdout.flush()
            pipe.log.info('Saving to alpha_significance_evenness.\n'
                          'ags_even.visualization.save(pipe.outputs_dir+"alpha_group_significance_evenness")')
            ags_even.visualization.save(pipe.outputs_dir+'alpha_group_significance_evenness')
            pipe.log.info('Done.\n')

            """
            qiime diversity alpha-group-significance
            --i-alpha-diversity qiime2-core-metrics-results/shannon_vector.qza
            --m-metadata-file Mapping_file_corrected.txt.nogz
            --o-visualization qiime2-core-metrics-results/shannon-group-significance.qzv
            """
            pipe.log.info('Creating Alpha Group Significance for Shannon index.\n'
                          'ags_shannon = alpha_group_significance(cm.shannon_vector, md)')
            ags_shannon = alpha_group_significance(cm.shannon_vector, md)
            sys.stdout.flush()
            pipe.log.info('Saving to alpha_significance_shannon.\n'
                          'ags_shannon.visualization.save(pipe.outputs_dir+"alpha_group_significance_shannon")')
            ags_shannon.visualization.save(pipe.outputs_dir+'alpha_group_significance_shannon')
            pipe.log.info('Done.\n')


        pipe.log.info('Pipeline done.')
    except ValueError as ve:
        exit_status = 1
        pipe.log.error(traceback.format_exc())
        if ve.args[0].endswith('Verify your table is valid and that you '
                               'provided a shallow enough sampling depth.'):
            pipe.log.error('Pipeline Error:')
            pipe.log.error(ve.args[0])
            msg = "From core_metrics: " + ve.args[0]
            pipe.log_to_db(job_id=pipe.job_id,
                           stack=traceback.format_exc(), msg=msg)
        elif ve.args[0].endswith('consist of exactly one value.'):
            pipe.log.error('Pipeline Error:')
            pipe.log.error(ve.args[0])
            msg = "From alpha_group_significance: " + ve.args[0]
            pipe.log_to_db(job_id=pipe.job_id,
                           stack=traceback.format_exc(), msg=msg)
        else:
            pipe.log.error('Pipeline Error')
            pipe.log.error(ve.args[0])
            pipe.log_to_db(job_id=pipe.job_id,
                           stack=traceback.format_exc(),
                           msg=str(ve))
    except BaseException:
        sys.stderr.write(traceback.format_exc())
        pipe.log.error('Pipeline Error:')
        pipe.log.error(traceback.format_exc())
        pipe.log_to_db(job_id=pipe.job_id, stack=traceback.format_exc(),
                       msg=traceback.format_exc())
        exit_status = 1

    finally:
        try:
            sys.stderr = origerror
            sys.stdout = origout
        except BaseException:
            pass
        pipe.log.info(exit_status)
        exit(exit_status)


class DSAnalysis16S(PipeBase):
    """Start pipe class definition.
    """

    def __init__(self, args):
        super().__init__(args, inputs_dir=args.inputs_dir)


if __name__ == '__main__':
    PARSER = argparse.ArgumentParser()
    PARSER.add_argument("--job_id", type=str,
                        help="job_id when running in Nephele. either "
                        "job_id or inputs_dir should be specified.")
    PARSER.add_argument("-i", "--inputs_dir", type=str,
                        help="input directory for running outside Nephele. "
                        "either job_id or inputs_dir should be specified.")
    PARSER.add_argument("-o", "--outputs_dir", type=str,
                        help="output directory for running outside Nephele. optional")
    PARSER.add_argument('--data_type', help='Ignore - placeholder.')
    PARSER.add_argument('-a', '--alpha_group_sig', action='store_true',
                        help="run alpha group significance.")

    REQ_PARAMS = PARSER.add_argument_group('required arguments')

    REQ_PARAMS.add_argument('--biom_fp', type=argparse.FileType('r'),
                            required=True, help='required')
    REQ_PARAMS.add_argument('--map_file', type=argparse.FileType('r'),
                            required=True, help='required')
    REQ_PARAMS.add_argument('--sampling_depth', type=int, required=True,
                            help='required')

    ARGS = PARSER.parse_args()
    main(ARGS)
