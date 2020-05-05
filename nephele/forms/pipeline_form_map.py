#!/usr/bin/env python
# -*- coding: utf-8 -*-

from nephele2.nephele.forms import qiime_forms, mothur_forms, forms


PIPELINE_NAME_TO_WEB_FORM = {
    "QIIME 16S FASTQ Paired End": qiime_forms.QiimePEForm,
    "mothur 16S FASTQ Paired End": mothur_forms.MothurPEForm,
    "QIIME 16S FASTQ Single End": qiime_forms.QiimeSEForm,
    "DADA2 16S FASTQ Paired End": forms.DadaPEForm,
    "DADA2 16S FASTQ Single End": forms.DadaSEForm,
    "QIIME ITS FASTQ Paired End": qiime_forms.QiimePEITSForm,
    "WGS Paired End FASTQ": forms.WGSbioBakeryForm,
    "WGS Single End FASTQ": forms.WGSbioBakeryForm,
    "Paired End QC": forms.QCPairedEndForm,
    "Single End QC": forms.QCSingleEndForm,
    "Downstream Analysis": forms.DSAnalysisForm}
