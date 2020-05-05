DADA2 Pipeline Documentation
============================

Setup
-----
-  Requires `rpy2 <https://rpy2.bitbucket.io>`__ python library
-  This pipeline runs the dada2nephele R package.  That `README <dada2nephele/README.md>`_ (`sphinx link <pipelines.DADA2.dada2nephele.readme>`) has the information about R deps and databases.
- it also runs the datavis16s R package for viz, so see `datavis16s README <pipelines.datavis16s.readme>` .

Notes
-----

-  See script usage ↓↓↓ for parameter data type.
-  See `user details page <dada2nephele/doc/user_doc.md>`_ (`sphinx link <pipelines.DADA2.dada2nephele.user_doc>`) for the help text, defaults and expected output.
-  The mapping file requirements are the same as those for the other 16S pipelines.


Usage
-----
.. code-block:: bash

        usage: dada2nephele.py [-h] [--job_id str] [-i str] [-o str] --map_file FILE
                       [--data_type {SE,PE}] [--trimleft_fwd int] [--trimleft_rev int]
                       [--maxee int] [--trunclen_fwd int] [--trunclen_rev int]
                       [--truncq int] [--ion_torrent] [--just_concatenate] [--maxmismatch int]
                       [--trim_overhang] [--chimera] [--ref_db {sv99,homd,gg97}]
                       [--taxmethod {rdp,idtaxa}] [--sampling_depth int]


job arguments:
~~~~~~~~~~~~~~
-h, --help                 show this help message and exit
--job_id str               if running inside Nephele
-i str, --inputs_dir str   for running outside of Nephele
-o str, --outputs_dir str  for running outside of Nephele
--map_file FILE            required
--data_type                {SE,PE} single end or paired end. (default: PE)

user options:
~~~~~~~~~~~~~
Filter and Trim
###############
--trimleft_fwd int         dada2::filterAndTrim trimLeft bp of fwd read. (default: 20)
--trimleft_rev int         dada2::filterAndTrim trimLeft bp of rev read. (default: 20)
--maxee int                dada2::filterAndTrim discard reads with EE higher than this value.
                           (default: 5)
--trunclen_fwd int         dada2::filterAndTrim truncate fwd reads at this length. (default: 0)
--trunclen_rev int         dada2::filterAndTrim truncate rev reads at this length. (default: 0)
--truncq int               dada2::filterAndTrim truncate reads at first bp with qual <= this
                           value. (default: 4)

Denoising
#########
--ion_torrent              dada2::dada use recommended parameters for ion torrent data (default: False)

Merge Pairs *paired-end only*
#############################
--just_concatenate         dada2::mergePairs concatenates instead of merging reads (default: False)
--maxmismatch int          dada2::mergePairs max mismatches allowed. (default: 0)
--trim_overhang            dada2::mergePairs trims overhanging seq. (default: False)


Analysis
########
--chimera                  run dada2::removeBimeraDenovo. (default: False)
--ref_db                   {sv99,homd,gg97} reference database (default: sv99)
--taxmethod                {rdp,idtaxa} taxonomic assignment method (default: rdp)
--sampling_depth int       sampling depth for downstream analysis. optional
