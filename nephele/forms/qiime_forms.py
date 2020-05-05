"""The QIIME argument forms"""
from wtforms import SelectField, BooleanField, FloatField, IntegerField
from wtforms.validators import NumberRange, Optional
from nephele2.nephele.forms.validators import greater_than
from nephele2.nephele.forms import qiime_forms_help as qiime_help
from nephele2.nephele.forms.forms import JobDetailsForm


class QiimeBaseOptions(JobDetailsForm):
    """
    Defines the options that are common to (almost) all QIIME pipes.
    """
    phred_quality = IntegerField(
        "Minimum phred quality score",
        default=19,
        description=qiime_help.qiime_base_phred_qual_desc)
    phred_offset = IntegerField(
        "Phred offset",
        default=33,
        description=qiime_help.qiime_base_phred_offset_desc)
    max_n = IntegerField(
        "Maximum ambiguous",
        default=0,
        description=qiime_help.qiime_pe_max_ambiguous_desc)
    max_bad_run = IntegerField(
        "Max bad run length",
        default=3,
        description=qiime_help.qiime_max_bad_run_length_desc)
    ref_db = SelectField("OTU Picking Database",
                         default="sv99",
                         description=qiime_help.qiime_base_database_desc,
                         choices=[("sv99", "SILVA 99 v132"),
                                  ("sv97", "SILVA 97 v132"),
                                  ("gg99", "Greengenes 99 v13.5"),
                                  ("gg97", "Greengenes 97 v13.5"),
                                  ("homd", "HOMD")])
    picrust = BooleanField("Run PICRUSt annotations",
                           default=False,
                           description=qiime_help.qiime_picrust_desc)
    otu_strategy = SelectField(
        "Analysis type",
        default="OPEN_REFERENCE",
        description=qiime_help.qiime_base_analysis_type_desc,
        choices=[("open", "Open Reference"),
                 ("closed", "Closed Reference"),
                 ("de_novo", "De Novo")])
    sampling_depth = IntegerField(
        "Sampling depth for downstream analysis",
        description=qiime_help.sampling_depth_desc,
        validators=[Optional(),
                    greater_than(message="Must be greater than or equal to 0.",
                                 minimum=0)])


class QiimePEForm(QiimeBaseOptions):
    """
    Defines the pipeline options form for QIIME FASTQ Paired-End.
    """
    min_overlap = IntegerField(
        "Minimum overlap length",
        default=10,
        description=qiime_help.qiime_min_overlap_desc)
    perc_max_diff = IntegerField(
        "Percent difference within overlap",
        default=25,
        description=qiime_help.qiime_percent_max_diff_desc,
        validators=[NumberRange(min=1,
                                max=100,
                                message="Must be an integer between 1-100")])


class QiimeSEForm(QiimeBaseOptions):
    """
    Defines the pipeline options for QIIME SE.
    """


class QiimePEITSForm(JobDetailsForm):
    """
    Defines the form for the QIIME ITS PE pipe.
    """
    phred_quality = IntegerField(
        "Minimum phred quality score",
        default=19,
        description=qiime_help.qiime_base_phred_qual_desc)
    phred_offset = IntegerField(
        "Phred offset",
        default=33,
        description=qiime_help.qiime_base_phred_offset_desc)
    max_n = IntegerField(
        "Maximum ambiguous",
        default=0,
        description=qiime_help.qiime_pe_max_ambiguous_desc)
    max_bad_run = IntegerField(
        "Maximum bad run length",
        default=3,
        description=qiime_help.qiime_its_max_bad_run_length_desc)
    min_overlap = IntegerField(
        "Minimum overlap length",
        default=10,
        description=qiime_help.qiime_min_overlap_desc)
    perc_max_diff = IntegerField(
        "Percent difference within overlap",
        default=25,
        description=qiime_help.qiime_percent_max_diff_desc)
    ref_db = SelectField(
        "Reference Database",
        default="ITS_99",
        description=qiime_help.qiime_its_ref_db_desc,
        choices=[("its99", "ITS 99"),
                 ("its97", "ITS 97")])
    sampling_depth = IntegerField(
        "Sampling depth for downstream analysis",
        description=qiime_help.sampling_depth_desc,
        validators=[Optional(),
                    greater_than(
                        message="Must be greater than or equal to 0.",
                        minimum=0)])
    otu_strategy = SelectField(
        "Analysis type",
        default="OPEN_REFERENCE",
        description=qiime_help.qiime_base_analysis_type_desc,
        choices=[("open", "Open Reference"),
                 ("closed", "Closed Reference"),
                 ("de_novo", "De Novo")])
