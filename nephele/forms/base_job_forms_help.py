"""
Defines the help text common to most of our pipeline forms.
"""

job_name_desc = """
    Please provide a brief description of the project for easier reference.
	"""

keep_data_desc = """
    User data is typically stored in our system for <strong>90</strong> days,
    during which you can view and download the data. However, we store data
    from the labs we collaborate with for a longer period of time. If you are
    one of our collaborators, check this box.
	"""

run_hmp_desc = """
    Check this box if you want to compare your analysis to the
    <a href='http://hmpdacc.org/HMR16S/' target='_blank' rel="noopener noreferrer">HMP DACC</a>.<br>
    Note: this will only work with human data.
	"""

dacc_region_desc = """
    This option allows comparing the generated <a href='http://biom-format.org/' target="_blank" rel="noopener noreferrer">
    Biological Observation Matrix (BIOM)</a> file from your analysis to the
    <a href='http://hmpdacc.org/HMR16S/' target='_blank' rel="noopener noreferrer">HMP DACC</a>. Choose
    V1-V3 for regions 1 to 3, V3-V5 for regions 3 to 5, or V6-V9 for regions
    6 to 9.
	"""

hmp_reference_db_desc = """
    <strong>SILVA</strong>: the SILVA database will be used for taxonomic assignments.<br>
    <strong>Greengenes</strong>: the Greengenes database will be used for taxonomic assignments.
    <ul>
        <li>97: sequences were clustered at 97% identity.</li>
        <li>99: sequences were clustered at 99% identity.</li>
    </ul>
	"""

hmp_body_site_desc = """
    This option allows comparing the generated Biological Observation Matrix (BIOM)
    file from your data to the <a href='http://hmpdacc.org/HMR16S/' target='_blank' rel="noopener noreferrer">
    HMP DACC</a>.<br>
    Note: if you use V6-V9, the following body sites are not available:
    <ul>
        <li>Urogenital tract</li>
        <li>Mid Vagina</li>
        <li>Vaginal introitus</li>
        <li>Posterior fornix</li>
    </ul>
	"""

hmp_nearest_n_samples_desc = """
    Please choose number of samples (range from 1 to 10) you would like to compare
    and display. Recommended number is 7.
	"""
