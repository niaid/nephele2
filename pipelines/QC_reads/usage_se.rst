Single-end Usage and Help Text
==============================

- **Single-end** runs the same script as `paired-end <./usage_pe.rst>`_ with a reduced set of options.
  
qc_reads.py [-h] [--job_id str] [--inputs_dir str] --map_file filename --data_type {PE,SE}
                   [--run_cutadapt] [--error_rate float] [--indels] [--overlap int]
                   [--match_read_wildcards] [--match_adapter_wildcards] [--adapter_f str]
                   [--front_f str] [--anywhere_f str]
                   [--run_qual_trimming] [--window_size int] [--req_qual int]
                   [--lead_qual int] [--trail_qual int] [--minlen int] [--avg_qual AVGQUAL]
                   
------                

- You can't call this script directly due to qiime2 dependency; instead use qc_reads.sh_ , but the arguments and usage are otherwise the same.
- The help text is taken from the help of the individual tools.
- The text in bold is what the option name should be on the front end, and the rest
  should be in the help text box.

job arguments:
--------------

-h, --help                       show this help message and exit
--job_id str                     job_id when running in Nephele
--inputs_dir str                 input directory for running outside Nephele
--map_file filename              required
--data_type                      {PE,SE} required

Adapter/Primer Trimming:
------------------------

--run_cutadapt                   Master Switch for running QIIME2 cudadapt. If this flag is set,
                                   also need to spec one or more of below params (default: False)
--error_rate float               **Error rate**: Maximum allowed error rate.  (default: 0.1, range: [0,1))
--indels                         **Indels**: Allow insertions or deletions of bases when
                                   matching adapters.  (default: True)
--overlap int                    **Adapter overlap**: Require at least this many bases of overlap
                                   between read and adapter for an adapter to
                                   be found.  (default: 3)
--match_read_wildcards           **Match read wildcards**: Interpret IUPAC wildcards (e.g., N) in
                                   reads.  (default: False)
--match_adapter_wildcards        **Match adapter wildcards**: Interpret IUPAC wildcards (e.g., N) in
                                   adapters.  (default: True)
                                   
Adapter strings
~~~~~~~~~~~~~~~
These strings must validate as valid IUPAC strings with the addition of the ``$`` and ``^`` characters.  **At least
one of the following must be non-empty!**  If left empty, pass *None*.

--adapter_f str                  **3' adapter**: Sequence of an adapter ligated to the 3’ end. The adapter and any
                                   subsequent bases are trimmed. If a ``$`` is appended, the adapter is only found if it
				   is at the end of the read. Maximum 250bp.
--front_f str                    **Front 5' adapter**: Sequence of an adapter or primer ligated to the 5’ end. The
                                   adapter and any preceding bases are trimmed. Partial matches at the 5’ end are allowed.
				   If a ``^`` character is prepended, the adapter is only found if it is at the beginning of the read.
				   Maximum 250bp.
--anywhere_f str                 **Anywhere adapter**: Sequence of an adapter that may be ligated to the 5’ or 3’ end. Both types of matches as described under ``3' adapter`` and ``Front 5' adapter`` are allowed. If the first base of the read is part of the match, the behavior is as with ``Front 5'``, otherwise as with ``3' adapter``. This option is mostly for rescuing failed library preparations - do not use if you know which end your adapter was ligated to. Maximum 250bp.

Quality Trimming:
-----------------
--run_qual_trimming              Master Switch for running Trimmomatic. If this flag is used, also
                                   need to spec one or more of below params (default: False)
--window_size int                **Window size**: Sliding window size
                                   for quality trimming, cutting once the average quality within the window falls
				   below the ``required quality``. Specifies the number of bases to average across. (default: 4) 
--req_qual int                   **Required quality**: specifies the average quality required in the sliding window (default: 12)
--lead_qual int                  **Leading quality**: Remove  low  quality  bases from the beginning (5' end) of the read.  
                                 As long as a base has a value below this threshold
                                 the base is removed and the next base will be investigated. (default: 0)
--trail_qual int                 **Trailing quality**: Remove low quality bases from the end of the read.  As long as a base has a 
                                 value below this threshold the base is  removed and the next base (which as Trimmomatic is starting    
                                 from the 3' end would be base preceding the just removed base) will be investigated. (default: 0)
--minlen int                     **Minimum length**: Remove reads that fall below the specified minimal length. (default: 60)
--avg_qual int                   **Average quality**: Drop the read if the average quality across the entire read is below 
                                 the specified level. (default: 0)

.. _qc_reads.sh: https://github.niaid.nih.gov/bcbb/nephele2/blob/next_release/pipelines/QC_reads/qc_reads.sh
