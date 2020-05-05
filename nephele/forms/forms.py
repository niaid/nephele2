"""
WTForms definitions of the forms used in the site.
"""
import re
from flask_wtf import FlaskForm
from flask_wtf.file import FileField, FileRequired, FileAllowed
from wtforms import StringField, BooleanField, RadioField, IntegerField,\
    SelectField, SelectMultipleField, FloatField, widgets
from wtforms.validators import DataRequired, InputRequired, \
    NumberRange, Optional, Length, Regexp
from nephele2.nephele.forms import validators as neph_validators
from nephele2.nephele.forms import base_job_forms_help as base_help
from nephele2.nephele.forms import dada_forms_help as dada_help
from nephele2.nephele.forms import wgs_forms_help as wgs_help
from nephele2.nephele.forms import qc_forms_help as qc_help


# ############# custom fields ################


class MultiCheckboxField(SelectMultipleField):
    """
    Field that renders a SelectMultipleField as a set
    of checkboxes.
    """
    widget = widgets.ListWidget(prefix_label=False)
    option_widget = widgets.CheckboxInput()


class MultiDashedCheckboxField(MultiCheckboxField):
    """
    A special variant of a MultiCheckboxField that returns
    a string of values concatenated with a dash (-) instead of
    an array.
    """

    def _value(self):
        if self.data:
            return [x.strip() for x in self.data.split('-')]
        return []

    def post_validate(self, form, validation_stopped):
        if self.data:
            self.data = "-".join(self.data)
        else:
            self.data = ""

# ############# forms ################


class JobDetailsForm(FlaskForm):
    """
    Defines the options that need to be displayed in all 16S pipeline options
    forms, but not passed to the pipeline when it's run. These options are
    stored in the database for other reasons.
    """
    job_desc = StringField("Description of the job",
                           default=None,
                           description=base_help.job_name_desc)


class DSAnalysisForm(JobDetailsForm):
    gt_five = "Input is required, must an integer greater than 5."
    desc = "The total frequency that each sample should be rarefied to prior "\
        "to computing diversity metrics. Samples with counts below this "\
        "value will be removed from the analysis."
    sampling_depth = IntegerField(
        "Sampling depth for downstream analysis",
        description=desc,
        validators=[InputRequired(),
                    neph_validators.greater_than(
                        message=gt_five, minimum=4)])

    red_warn = '<span style="color: red;">Requires Treatment Group column in \
        mapping file, after filtering by sampling depth, to contain at least \
        2 groups, and each group to contain at least 2 samples.</span>'
    link = '<a class="" href="da_details/#input-files-and-parameters" \
        target="_blank" rel="noopener noreferrer">More information.</a>'
    alpha_desc = "Run alpha diversity statistical comparisons between groups, \
        and produce alpha diversity plots. {} {}".format(red_warn, link)
    alpha_group_sig = BooleanField("Run Alpha group significance",
                                   default=False,
                                   description=alpha_desc)


class DadaPEForm(JobDetailsForm):
    """
    Defines the form for the DADA2 paired-end pipeline.
    """
    zero_or_g = 'Must be 0 or greater.'
    gt_zero = "Must be greater than 0."
    ion_torrent = BooleanField("Check if Ion Torrent sequence",
                               default=False,
                               description=dada_help.ion_torrent_seq)
    trimleft_fwd = IntegerField(
        "Trim left forward",
        default=20,
        description=dada_help.trim_left_fwd_desc,
        validators=[
            neph_validators.greater_than(
                message=zero_or_g,
                minimum=-1)])
    trimleft_rev = IntegerField(
        "Trim left reverse",
        default=20,
        description=dada_help.trim_left_rev_desc,
        validators=[
            neph_validators.greater_than(
                message=zero_or_g,
                minimum=-1)])
    truncQ = IntegerField(
        "Truncation quality score",
        default=4,
        description=dada_help.tuncq_desc,
        validators=[
            neph_validators.greater_than(
                message=gt_zero,
                minimum=-1)])
    truncLen_fwd = IntegerField(
        "Truncation length forward", default=0,
        description=dada_help.truncLen_fwd_desc,
        validators=[neph_validators.greater_than(message=gt_zero, minimum=-1)])
    truncLen_rev = IntegerField(
        "Truncation length reverse", default=0,
        description=dada_help.truncLen_rev_desc,
        validators=[neph_validators.greater_than(message=gt_zero, minimum=-1)])
    maxEE = IntegerField(
        "Maximum expected errors",
        default=5,
        description=dada_help.maxee_desc,
        validators=[
            neph_validators.greater_than(
                message=gt_zero,
                minimum=-1)])
    just_concatenate = BooleanField("Just concatenate",
                                    default=False,
                                    description=dada_help.just_concatenate)
    maxMismatch = IntegerField(
        "Maximum mismatches",
        default=0,
        description=dada_help.max_mismatches_desc,
        validators=[
            neph_validators.greater_than(
                message=gt_zero,
                minimum=-1)])
    trim_overhang = BooleanField("Trim overhanging sequence",
                                 default=False,
                                 description=dada_help.trim_overhang_desc)
    chimera = BooleanField("Run chimera removal",
                           default=True,
                           description=dada_help.remove_chimera_desc)
    taxmethod = SelectField("Taxonomic assigment",
                            default="rdp",
                            description=dada_help.tax_method,
                            choices=[("rdp", "rdp"), ("idtaxa", "IDTAXA")])
    sampling_depth = IntegerField("Sampling depth for downstream analysis",
                                  description=dada_help.sampling_depth_desc,
                                  validators=[Optional(),
                                              neph_validators.greater_than(
                                                  message=zero_or_g,
                                                  minimum=0)])
    ref_db = SelectField("Reference Database",
                         default="SILVA v132",
                         description=dada_help.ref_db_desc,
                         choices=[("sv99", "SILVA v132"),
                                  ("gg97", "Greengenes v13.8"),
                                  ("homd", "HOMD")])


class DadaSEForm(JobDetailsForm):
    """
    Defines the form for the DADA2 single-end pipeline.
    """
    ion_torrent = BooleanField("Check if Ion Torrent sequence",
                               default=False,
                               description=dada_help.ion_torrent_seq)
    zero_or_g = 'Must be 0 or greater.'
    trimleft_fwd = IntegerField(
        "Trim left",
        default=20,
        description=dada_help.trim_left_fwd_desc,
        validators=[
            neph_validators.greater_than(
                message=zero_or_g,
                minimum=-1)])
    truncQ = IntegerField(
        "Truncation quality score",
        default=4,
        description=dada_help.tuncq_desc,
        validators=[
            neph_validators.greater_than(
                message="Must be greater than 0.",
                minimum=-1)])
    truncLen_fwd = IntegerField(
        "Truncation length",
        default=0,
        description=dada_help.truncLen_fwd_desc,
        validators=[
            neph_validators.greater_than(
                message="Must be greater than 0.",
                minimum=-1)])
    maxEE = IntegerField(
        "Maximum expected errors",
        default=5,
        description=dada_help.maxee_desc,
        validators=[
            neph_validators.greater_than(message="Must be greater than 0.",
                                         minimum=-1)])
    chimera = BooleanField("Run chimera removal",
                           default=True,
                           description=dada_help.remove_chimera_desc)
    taxmethod = SelectField("Taxonomic assigment",
                            default="rdp",
                            description=dada_help.tax_method,
                            choices=[("rdp", "rdp"), ("idtaxa", "IDTAXA")])
    sampling_depth = IntegerField(
        "Sampling depth for downstream analysis",
        description=dada_help.sampling_depth_desc,
        validators=[
            Optional(),
            neph_validators.greater_than(
                message="Must be greater than or equal to 0.",
                minimum=0)])
    ref_db = SelectField(
        "Reference Database", default="SILVA v132",
        description=dada_help.ref_db_desc,
        choices=[("sv99", "SILVA v132"),
                 ("gg97", "Greengenes v13.8"),
                 ("homd", "HOMD")])


class WGSbioBakeryForm(JobDetailsForm):
    """
    Defines the form for the WGS pipeline.
    """
    strainphlan = BooleanField("Run strainphlan",
                               default=False,
                               description=wgs_help.run_strainphlan_desc)
    project_name = StringField("Project name for visualization pipeline.",
                               description=wgs_help.project_name_desc)


class QCPairedEndForm(JobDetailsForm):
    """
    Defines the form for the paired end QC form.

    NOTE: The QC forms UIs are different enough from the other forms in the
    system that they do not inherit the options_base template.
    """
    run_cutadapt = BooleanField("Run adapter trimming",
                                default=False,
                                description=qc_help.run_cutadapt_desc)
    run_qual_trimming = BooleanField(
        "Run quality trimming", default=False,
        description=qc_help.run_qual_trimming_desc)
    run_flash2_merge = BooleanField("Merge read pairs",
                                    default=False,
                                    description=qc_help.run_flash2_merge_desc)
    error_rate = FloatField(
        "Error rate",
        default=0.1,
        description=qc_help.error_rate_desc,
        validators=[
            NumberRange(
                min=0,
                message="Must be greater than or equal to %(min)s."),
            neph_validators.RequiredIf('run_cutadapt')])
    indels = BooleanField("Indels",
                          default=True,
                          description=qc_help.indel_desc)
    overlap = IntegerField(
        "Adapter overlap",
        default=3,
        description=qc_help.overlap_desc,
        validators=[
            NumberRange(
                min=0,
                message="Must be greater than or equal to %(min)s."),
            neph_validators.RequiredIf('run_cutadapt')])
    match_read_wildcards = BooleanField("Match read wildcards",
                                        default=True,
                                        description=qc_help.wildcard_read_desc)
    match_adapter_wildcards = BooleanField(
        "Match adapter wildcards", default=True,
        description=qc_help.wildcard_adapter_desc)
    adapter_f = StringField(
        "Forward 3' adapter",
        description=qc_help.adapter_desc,
        validators=[
            Length(
                max=250,
                message="Length must be less than %(max)d"),
            Regexp(
                regex=r'(?i)^[acgturyswkmbdhvn\.\-\$\^]+$',
                flags=re.IGNORECASE,
                message="Must be a valid IUPAC sequence."),
            Optional()])
    adapter_r = StringField(
        "Reverse 3' adapter",
        description=qc_help.adapter_desc,
        validators=[
            Length(
                max=250,
                message="Length must be less than %(max)d"),
            Regexp(
                regex=r'(?i)^[acgturyswkmbdhvn\.\-\$\^]+$',
                flags=re.IGNORECASE,
                message="Must be a valid IUPAC sequence."),
            Optional()])
    front_f = StringField(
        "Forward front 5' adapter",
        description=qc_help.front_desc,
        validators=[
            Length(
                max=250,
                message="Length must be less than %(max)d"),
            Regexp(
                regex=r'(?i)^[acgturyswkmbdhvn\.\-\$\^]+$',
                flags=re.IGNORECASE,
                message="Must be a valid IUPAC sequence."),
            Optional()])
    front_r = StringField(
        "Reverse front 5' adapter",
        description=qc_help.front_desc,
        validators=[
            Length(
                max=250,
                message="Length must be less than %(max)d"),
            Regexp(
                regex=r'(?i)^[acgturyswkmbdhvn\.\-\$\^]+$',
                flags=re.IGNORECASE,
                message="Must be a valid IUPAC sequence."),
            Optional()])
    anywhere_f = StringField(
        "Forward anywhere adapter",
        description=qc_help.anywhere_desc,
        validators=[
            Length(
                max=250,
                message="Length must be less than %(max)d"),
            Regexp(
                regex=r'(?i)^[acgturyswkmbdhvn\.\-\$\^]+$',
                flags=re.IGNORECASE,
                message="Must be a valid IUPAC sequence."),
            Optional()])
    anywhere_r = StringField(
        "Reverse anywhere adapter",
        description=qc_help.anywhere_desc,
        validators=[
            Length(
                max=250,
                message="Length must be less than %(max)d"),
            Regexp(
                regex=r'(?i)^[acgturyswkmbdhvn\.\-\$\^]+$',
                flags=re.IGNORECASE,
                message="Must be a valid IUPAC sequence."),
            Optional()])
    window_size = IntegerField(
        "Sliding window size",
        default=4,
        description=qc_help.window_size_desc,
        validators=[
            NumberRange(
                min=0,
                message="Must be greater than or equal to %(min)s."),
            neph_validators.RequiredIf('run_qual_trimming')])
    req_qual = IntegerField(
        "Required quality for sliding window",
        default=12,
        description=qc_help.req_qual_desc,
        validators=[
            NumberRange(
                min=0,
                message="Must be greater than or equal to %(min)s."),
            neph_validators.RequiredIf('run_qual_trimming')])
    lead_qual = IntegerField(
        "Leading quality",
        default=3,
        description=qc_help.lead_qual_desc,
        validators=[
            NumberRange(
                min=0,
                message="Must be greater than or equal to %(min)s."),
            neph_validators.RequiredIf('run_qual_trimming')])
    trail_qual = IntegerField(
        "Trailing quality",
        default=3,
        description=qc_help.trail_qual_desc,
        validators=[
            NumberRange(
                min=0,
                message="Must be greater than or equal to %(min)s."),
            neph_validators.RequiredIf('run_qual_trimming')])
    minlen = IntegerField(
        "Minimum length",
        default=30,
        description=qc_help.minlen_desc,
        validators=[
            NumberRange(
                min=0,
                message="Must be greater than or equal to %(min)s."),
            neph_validators.RequiredIf('run_qual_trimming')])
    avg_qual = IntegerField(
        "Average quality",
        default=0,
        description=qc_help.avg_qual_desc,
        validators=[
            NumberRange(
                min=0,
                message="Must be greater than or equal to %(min)s."),
            neph_validators.RequiredIf('run_qual_trimming')])
    f2_min_overlap = IntegerField(
        "Minimum overlap",
        default=10,
        description=qc_help.min_overlap_desc,
        validators=[
            NumberRange(
                min=0,
                message="Must be greater than or equal to %(min)s."),
            neph_validators.RequiredIf('run_flash2_merge')])
    f2_max_overlap = IntegerField(
        "Maximum overlap",
        default=315,
        description=qc_help.max_overlap_desc,
        validators=[
            NumberRange(
                min=0,
                message="Must be greater than or equal to %(min)s."),
            neph_validators.RequiredIf('run_flash2_merge')])
    f2_min_overlap_outie = IntegerField(
        "outie overlap",
        default=35,
        description=qc_help.min_outie_overlap_desc,
        validators=[
            NumberRange(
                min=0,
                message="Must be greater than or equal to %(min)s."),
            neph_validators.RequiredIf('run_flash2_merge')])
    f2_max_mismatch_density = FloatField(
        "Maximum mismatch density", default=0.25,
        description=qc_help.max_mismatch_density_desc,
        validators=[
                    NumberRange(
                        min=0,
                        message="Must be greater than or equal to %(min)s."),
                    neph_validators.RequiredIf('run_flash2_merge')])

    def validate(self):
        if not super(QCPairedEndForm, self).validate():
            return False
        if self.run_cutadapt.data and not self.adapter_f.data and not\
                self.adapter_r.data and not self.front_f.data and not\
                self.front_r.data and not self.anywhere_f.data and not\
                self.anywhere_r.data:
            msg = "A valid IUPAC sequence must be provided for at least one of\
            these fields."
            self.adapter_f.errors.append(msg)
            self.adapter_r.errors.append(msg)
            self.front_f.errors.append(msg)
            self.front_r.errors.append(msg)
            self.anywhere_f.errors.append(msg)
            self.anywhere_r.errors.append(msg)
            return False
        return True


class QCSingleEndForm(JobDetailsForm):
    """
    Defines the form for the single end QC form.

    NOTE: The QC forms UIs are different enough from the other forms in the
    system that they do not inherit the options_base template.
    """
    run_cutadapt = BooleanField("Run adapter trimming",
                                default=False,
                                description=qc_help.run_cutadapt_desc)
    run_qual_trimming = BooleanField(
        "Run quality trimming", default=False,
        description=qc_help.run_qual_trimming_desc)
    error_rate = FloatField(
        "Error rate",
        default=0.1,
        description=qc_help.error_rate_desc,
        validators=[
            NumberRange(
                min=0,
                message="Must be greater than or equal to %(min)s."),
            neph_validators.RequiredIf('run_cutadapt')])
    indels = BooleanField("Allow indels",
                          default=True,
                          description=qc_help.indel_desc)
    overlap = IntegerField(
        "Adapter overlap",
        default=3,
        description=qc_help.overlap_desc,
        validators=[
            NumberRange(
                min=0,
                message="Must be greater than or equal to %(min)s."),
            neph_validators.RequiredIf('run_cutadapt')])

    match_read_wildcards = BooleanField("Match read wildcards",
                                        default=True,
                                        description=qc_help.wildcard_read_desc)
    match_adapter_wildcards = BooleanField(
        "Match adapter wildcards", default=True,
        description=qc_help.wildcard_adapter_desc)
    adapter_f = StringField(
        "3' Adapter",
        description=qc_help.adapter_desc,
        validators=[
            Length(
                max=250,
                message="Length must be less than %(max)d"),
            Regexp(
                regex=r'(?i)^[acgturyswkmbdhvn\.\-\$\^]+$',
                flags=re.IGNORECASE,
                message="Must be a valid IUPAC sequence."),
            Optional()])
    front_f = StringField(
        "Front 5' adapter",
        description=qc_help.front_desc,
        validators=[
            Length(
                max=250,
                message="Length must be less than %(max)d"),
            Regexp(
                regex=r'(?i)^[acgturyswkmbdhvn\.\-\$\^]+$',
                flags=re.IGNORECASE,
                message="Must be a valid IUPAC sequence."),
            Optional()])
    anywhere_f = StringField(
        "Anywhere adapter",
        description=qc_help.anywhere_desc,
        validators=[
            Length(
                max=250,
                message="Length must be less than %(max)d"),
            Regexp(
                regex=r'(?i)^[acgturyswkmbdhvn\.\-\$\^]+$',
                flags=re.IGNORECASE,
                message="Must be a valid IUPAC sequence."),
            Optional()])
    window_size = IntegerField(
        "Sliding window size",
        default=4,
        description=qc_help.window_size_desc,
        validators=[
            NumberRange(
                min=0,
                message="Must be greater than or equal to %(min)s."),
            neph_validators.RequiredIf('run_qual_trimming')])
    req_qual = IntegerField(
        "Required quality for sliding window",
        default=12,
        description=qc_help.req_qual_desc,
        validators=[
            NumberRange(
                min=0,
                message="Must be greater than or equal to %(min)s."),
            neph_validators.RequiredIf('run_qual_trimming')])
    lead_qual = IntegerField(
        "Leading quality",
        default=3,
        description=qc_help.lead_qual_desc,
        validators=[
            NumberRange(
                min=0,
                message="Must be greater than or equal to %(min)s."),
            neph_validators.RequiredIf('run_qual_trimming')])
    trail_qual = IntegerField(
        "Trailing quality",
        default=3,
        description=qc_help.trail_qual_desc,
        validators=[
            NumberRange(
                min=0,
                message="Must be greater than or equal to %(min)s."),
            neph_validators.RequiredIf('run_qual_trimming')])
    minlen = IntegerField(
        "Minimum length",
        default=30,
        description=qc_help.minlen_desc,
        validators=[
            NumberRange(
                min=0,
                message="Must be greater than or equal to %(min)s."),
            neph_validators.RequiredIf('run_qual_trimming')])
    avg_qual = IntegerField(
        "Average quality",
        default=0,
        description=qc_help.avg_qual_desc,
        validators=[
            NumberRange(
                min=0,
                message="Must be greater than or equal to %(min)s."),
            neph_validators.RequiredIf('run_qual_trimming')])

    def validate(self):
        if not super(QCSingleEndForm, self).validate():
            return False
        if self.run_cutadapt.data and not self.adapter_f.data and not\
                self.front_f.data and not self.anywhere_f.data:
            msg = "A valid IUPAC sequence must be provided for at least one of\
            these fields."
            self.adapter_f.errors.append(msg)
            self.front_f.errors.append(msg)
            self.anywhere_f.errors.append(msg)
            return False
        return True


class UploadSeqUrlForm(FlaskForm):
    """
    Defines the form rendered with the upload.html template.
    """
    ALLOWED_PREFIXES = ['ftp']
    # TODO: validate that if the file is an archive it contains files with
    # valid extensions.
    data_file = StringField(
        'Please provide the URL for your file: ',
        [InputRequired("You must provide a valid URL"),
         neph_validators.url_valid,
         neph_validators.url_file_allowed(
             ALLOWED_PREFIXES,
             'Invalid file type. You may only upload files with extensions in '
             + ', '.join(ALLOWED_PREFIXES))])


class UploadBiomForm(FlaskForm):
    """
    Defines the form rendered with the upload.html template.
    """
    ALLOWED_EXTS = ['biom']
    EXT_ERR_MSG = 'Invalid file type. You may only upload %s files.' + \
        ', '.join(ALLOWED_EXTS)
    data_file = FileField(
        'Please select file: ', [
            FileRequired(), FileAllowed(
                ALLOWED_EXTS, EXT_ERR_MSG)])


class UploadMapForm(FlaskForm):
    """
    Defines the form rendered with the upload.html template.
    """
    ALLOWED_EXTENSIONS = ['txt', 'csv', 'xlsx']
    ERR_MSG = 'Invalid file type. You may only upload files '\
        'with extensions: {}'.format(', '.join(ALLOWED_EXTENSIONS))
    data_file = FileField('Please select file: ',
                          [FileRequired(),
                           FileAllowed(ALLOWED_EXTENSIONS, ERR_MSG)])


class UploadQualFileForm(FlaskForm):
    """
    Defines the form rendered with the upload.html template for uploading the
    qual file.
    """
    ALLOWED_EXTENSIONS = ['qual']
    data_file = FileField(
        'Please select file: ',
        [FileRequired(),
         FileAllowed(
             ALLOWED_EXTENSIONS,
             'Invalid file type. You may only upload files with extensions in '
             + ', '.join(ALLOWED_EXTENSIONS))])


class UploadBarcodeFileForm(FlaskForm):
    """
    Defines the form rendered with the upload.html template for uploading a
    barcode file.
    """
    ERR_MSG = 'Invalid file type. You may only upload fastq files'
    data_file = FileField('Please select file: ',
                          [FileRequired(),
                           FileAllowed(['fastq'], ERR_MSG)])


class DataTypeSelectionForm(FlaskForm):
    """
    Defines the form rendered with the datatype_selection.html template for
    16S.
    """
    data_type = RadioField('filetype',
                           choices=[("PE", "Paired End FASTQ"),
                                    ("SE", "Single End FASTQ"),
                                    ],
                           validators=[DataRequired("Please select the type \
                                       of data you will be using.")])


class DataTypeQCSelectionForm(FlaskForm):
    """
    Defines the form rendered with the datatype_selection.html template for
    16S.
    """
    data_type = RadioField('filetype',
                           choices=[("QC_PE", "Paired End FASTQ"),
                                    ("QC_SE", "Single End FASTQ"),
                                    ],
                           validators=[DataRequired("Please select the type \
                                       of data you will be using.")])


class DataTypeWGSSelectionForm(FlaskForm):
    """
    Defines the form rendered with the datatype_selection.html template for
    ITS.
    """
    data_type = RadioField('filetype',
                           choices=[("WGS_PE", "WGS Paired End FASTQ"),
                                    ("WGS_SE", "WGS Single End FASTQ")],
                           validators=[DataRequired("Please select the type \
                                       of data you will be using.")])
