KBYTE = 1024

#: in kilobytes per du man page
MAX_MOTHUR_INPUT = 10000000

HUMAN_READABLE_SIZE = round(MAX_MOTHUR_INPUT/KBYTE**2, 2)

INPUTS_TOO_BIG_MSG = (
    'The inputs for this job are greater than {}GiB. mothur will not be able '
    'to analyse these data due to memory limitations. Please consider another '
    'pipeline.').format(HUMAN_READABLE_SIZE)

#: filename for biom summary
OTU_SUMMARY_TABLE = 'otu_summary_table.txt'

CTG_FNAME = 'combo.trim.contigs.renamed.good.unique.good.filter.unique.precluster'

BIOMFILE = CTG_FNAME + '.pick.opti_mcc.0.03.biom'
#: output files that will be sent to the user
FILES_TO_KEEP = ['logfile.txt',
                 OTU_SUMMARY_TABLE,
                 '{}.denovo.vsearch.pick.pick.count_table'.format(CTG_FNAME),
                 '{}.denovo.vsearch.pick.pick.count.summary'.format(CTG_FNAME),
                 '{}.pick.fasta'.format(CTG_FNAME),
                 '{}.pick.nr_v128.wang.pick.tax.summary'.format(CTG_FNAME),
                 BIOMFILE,
                 '{}.pick.opti_mcc.0.03.cons.tax.summary'.format(CTG_FNAME),
                 '{}.pick.opti_mcc.0.03.cons.taxonomy'.format(CTG_FNAME),
                 '{}.pick.opti_mcc.list'.format(CTG_FNAME),
                 '{}.pick.opti_mcc.shared'.format(CTG_FNAME),
                 '{}.pick.nr_v128.wang.pick.taxonomy'.format(CTG_FNAME),
                 '{}.pick.phylip.dist'.format(CTG_FNAME),
                 '{}.pick.phylip.tre'.format(CTG_FNAME)]

#: dict mapping reference database name to *ref_db* reference fasta and
#: *tax_db* reference taxonomy.  **Add new reference dbs here.**
DBS = {'homd': {'ref_db': 'HOMD_16S_rRNA_RefSeq_V15.1.aligned.fasta',
                'tax_db': 'HOMD_16S_rRNA_RefSeq_V15.1.mothur.taxonomy'},
       'sv99': {'ref_db': 'silva.nr_v128.align',
                'tax_db': 'silva.nr_v128.tax'}}

#: output directories sent to the user
DIRS_TO_KEEP = ['otus_picrust', 'PICRUSt_data', 'graphs']

UNKNOWN_ERR_MSG = 'An unknown error arose during execution of the pipeline'
