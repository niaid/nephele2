#!/usr/bin/env python3

import os
import csv
## need to add qiime2 site-packages to the path
from qiime2 import Artifact
from qiime2.plugins import cutadapt


class Q2Trim:

    @property
    def manifest(self):
        return self._manifest

    @staticmethod
    def gen_manifest(samples, out_dir):
        fname = out_dir+'input_manifest.txt'
        with open(fname, 'w') as mfest:
            m_writer = csv.writer(mfest)
            m_writer.writerow(['sample-id','absolute-filepath','direction'])
            for sample in samples:
                m_writer.writerow([sample.id, sample.fwd_fp, 'forward'])
                if sample.rev_fp:
                    m_writer.writerow([sample.id, sample.rev_fp, 'reverse'])
        return fname

    @staticmethod
    def gen_input_data(manifest_fp, data_type):
        if data_type == 'PE':
            inputdata = Artifact.import_data('SampleData[PairedEndSequencesWithQuality]',
                                             manifest_fp,
                                             view_type='PairedEndFastqManifestPhred33')
        else:
            inputdata = Artifact.import_data('SampleData[SequencesWithQuality]',
                                             manifest_fp,
                                             view_type='SingleEndFastqManifestPhred33')
        return inputdata

