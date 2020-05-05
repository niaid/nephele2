"""The mothur argument forms"""
from wtforms import BooleanField, IntegerField, SelectField
from wtforms.validators import Optional
from nephele2.nephele.forms.validators import greater_than
from nephele2.nephele.forms.forms import JobDetailsForm,\
    MultiDashedCheckboxField
from nephele2.nephele.forms import mothur_forms_help as mothur_help


class MothurPEForm(JobDetailsForm):
    """
    Defines the pipeline options form for Mothur FASTQ Paired-End.
    """
    maxlength = IntegerField(
        "Max length *",
        default=0,
        description=mothur_help.mothur_max_seq_length_desc)
    remove_lineage = MultiDashedCheckboxField(
        "Remove sequences that belong to taxa of certain lineages",
        default=["Chloroplast", "Mitochondria", "unknown", "Eukaryote"],
        description=mothur_help.mothur_remove_lineage_desc,
        choices=[("Archaea", "Archaea"),
                 ("unknown", "unknown"),
                 ("Chloroplast", "Chloroplast"),
                 ("Mitochondria", "Mitochondria"),
                 ("Eukaryote", "Eukaryote")])
    optimize = SelectField(
        "Optimize",
        default="start-end",
        description=mothur_help.mothur_optimize_criteria_desc,
        choices=[("start-end", "start-end"),
                 ("start", "start"),
                 ("end", "end")])
    criteria = IntegerField(
        "Criteria",
        default=90,
        description=mothur_help.mothur_optimize_criteria_desc)
    ref_db = SelectField(
        "OTU Picking Reference Database",
        default="SILVA",
        description=mothur_help.mothur_reference_db_desc,
        choices=[("sv99", "SILVA 99"), ("homd", "HOMD")])
    picrust = BooleanField("Run PICRUSt annotations",
                           default=False,
                           description=mothur_help.mothur_picrust_desc)
    sampling_depth = IntegerField(
        "Sampling depth for downstream analysis",
        description=mothur_help.sampling_depth_desc,
        validators=[Optional(),
                    greater_than(
                        message="Must be greater than or equal to 0.",
                        minimum=0)])
