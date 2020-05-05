QC Pipeline README
==================

-  Usage for the front end: `PE <./usage_pe.rst>`_ | `SE <./usage_se.rst>`_
-  User details: `github <./user_doc.md>`_ | `html <./qc_details.html>`_
.. contents:: :local:
   :depth: 2

Dependencies
------------

-  `qiime2 2019.7 <https://docs.qiime2.org/2019.7/install/>`__ - See `qiime2 backup <../../docs/source/qiime2_backup.rst>`__ (:doc:`sphinx </qiime2_backup>`) for how to install/restore this particular version.
-  `trimmomatic <http://www.usadellab.org/cms/?page=trimmomatic>`__
-  `FastQC <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`__
-  `multiQC <http://multiqc.info/docs/>`__
-  `FLASH <https://ccb.jhu.edu/software/FLASH/>`__
   (`FLASH2 <https://github.com/dstreett/FLASH2>`__ - also available in
   `conda <https://github.com/dstreett/FLASH2/issues/10#issuecomment-364373920>`__)

The above should be available on ami-05ab4c36218c30780

-  .. raw:: html

    <s><a href="https://github.com/linsalrob/fastq-pair">fastq-pair</a></s>

Steps
-----

-  The user can choose which of steps 1-4 to run (though, they will
   always run in this order)
-  In comments, I put actual-ish commands I used for example data in
   Locus:


0. MultiQC directory
~~~~~~~~~~~~~~~~~~~~

-  multiqc can take in log files for some of the tools, so we will make
   a directory for its input

   .. code:: bash

      mkdir ./outputs/multiqc_input

1. FastQC
~~~~~~~~~

-  no options
-  fastqc can take in fastq or fastq.gz without any need to set flag

   -  run fastqc command for each fastq file - either run multiple files
      in parallel (best performance for many files)
   -  or run each one multithreaded with -t (best performance for big
      files)

   .. code:: bash

      ## run for each file or in parallel
      fastqc -t numberofthreads -o outputs/multiqc_input inputs/fastqfile
      ## cat fastqlist | parallel -j 10 fastqc -o ./outputs/multiqc_input
      ## or
      ## cat fastqlist | xargs fastqc -o ./outputs/multiqc_input -t 10

2. Primer/adapter Trimming
~~~~~~~~~~~~~~~~~~~~~~~~~~

-  It might be easier to use the `Artifact
   API <https://docs.qiime2.org/2018.6/interfaces/artifact-api/>`__,
   than to run as shell commands. I started writing a little script
   `trim.py <./trim.py>`__
-  qiime2 import

   -  Should only allow phred+33 data (most common)? Otherwise we can
      have parameter for Phred score type

   -  will need to make manifest from mapping file - `manifest file
      format <https://docs.qiime2.org/2018.6/tutorials/importing/#id5>`__

   -  Paired

   -  .. code:: bash

         ## Older version
         qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path ./pairedmanifest.txt --output-path paired-end-demux.qza --source-format PairedEndFastqManifestPhred33
         ## Newer version of qiime2 - 2018.11 and forward use different syntax
         qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path ./pairedmanifest.txt --output-path paired-end-demux.qza --input-format PairedEndFastqManifestPhred33

   -  Single end
   -  .. code:: bash

         qiime tools import --type 'SampleData[SequencesWithQuality]' --input-path ./se_manifest.txt --output-path single-end-demux.qza --source-format SingleEndFastqManifestPhred33

-  `qiime2
   cutadapt <https://docs.qiime2.org/2018.6/plugins/available/cutadapt/>`__

   -  we can add all of the flags as user options eventually
   -  `Paired-end <https://docs.qiime2.org/2018.6/plugins/available/cutadapt/trim-paired/>`__

      -  The user must submit at least one of the following flags/options ``--p-adapter-f --p-adapter-r --p-front-f --p-front-r --p-anywhere-f --p-anywhere-r``
      -  Each flag takes a string as an argument. The string should validate as a sequence of `IUPAC characters <https://www.bioinformatics.org/sms/iupac.html>`__ (there are no spaces in the IUPAC code!) and the characters ``$`` and ``^``. Let's limit the length to 250bp?
      -  If left empty as a user option, the flag should not be used in the cutadapt command
      -  The other available flags are optional (and have defaults); let's just start with these, so we can test a bit.

      .. code:: bash

         ## examples
         qiime cutadapt trim-paired --verbose --i-demultiplexed-sequences paired-end-demux.qza --p-cores numberofthreads --p-front-f FWDPRIMER --p-front-r REVPRIMER --o-trimmed-sequences ./cutadapt_trimmed_seqs.qza
         # or
         qiime cutadapt trim-paired --verbose --i-demultiplexed-sequences paired-end-demux.qza --p-cores numberofthreads --p-adapter-f ADAPTER --o-trimmed-sequences ./cutadapt_trimmed_seqs.qza

   -  `Single end <https://docs.qiime2.org/2018.6/plugins/available/cutadapt/trim-single/>`__

      -  The user must submit at least one of the following flags/options ``--p-adapter --p-front --p-anywhere``
      -  Other bits are the same as paired-end

   -  .. code:: bash

         qiime cutadapt trim-single --verbose --i-demultiplexed-sequences ./outputs/single-end-demux.qza --p-cores numberofthreads --p-front FWDPRIMER --o-trimmed-sequences ./cutadapt_trimmed_seqs.qza

   -  Options in both pipelines (see `single end api <https://docs.qiime2.org/2018.6/plugins/available/cutadapt/trim-single/#api>`__ for data types)

      -  ``error_rate`` - 0.1 , ``indels`` - TRUE, ``overlap`` - 3, ``match_read_wildcards`` - FALSE, ``match_adapter_wildcards`` - TRUE

   -  it would be good if we could capture **stdout** from cutadapt to a separate file (separate from logfile.txt). Does ``tee`` work with our logging setup??? multiQC can also take in the cutadapt log ("real" cutadapt has command line option to divert output to file, qiime2 does not) and make some graphs

      .. code:: bash

         # qiime cutadapt trim-paired --verbose --i-demultiplexed-sequences paired-end-demux.qza --p-cores 30 --p-front-f TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGAGAGTTTGATCCTGGCTCAG --p-front-r GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGATTACCGCGGCTGCTGG --o-trimmed-sequences cutadapt_trimmed_seqs.qza >outputs/fastqc/cutadapt.log

-  **qiime2 export**

   .. code:: bash

      ## Older version
      qiime tools export --output-dir ./outputs/cutadapt_trimmed_seqs cutadapt_trimmed_seqs.qza
      # qiime2-2018.11 and newer uses the following
      qiime tools export --output-path ./outputs/cutadapt_trimmed_seqs --input-path cutadapt_trimmed_seqs.qza

   -  upon export, qiime names the files:

      sample-id_n_L001_R1_001.fastq.gz and
      sample-id_(n+1)_L001_R2_001.fastq.gz

      The prefixes for R1 and R2 are not identical (!!) - there is a
      and index 0 integer inserted. Luckily, there is also an
      output manifest file ./outputs/cutadapt_trimmed_seqs/MANIFEST
      which you can parse with the sample id and forward/reverse to get
      the corresponding filenames

3. Quality trimming
~~~~~~~~~~~~~~~~~~~

-  Trimmomatic has 5 options which we will let the user set (text copied
   from their website; defaults set by me):

   -  SLIDINGWINDOW ``windowSize:requiredQuality``

      -  windowSize (int): specifies the number of bases to average
         across (default 4)
      -  requiredQuality (int): specifies the average quality required
         (default 12).

   -  LEADING quality (int): Cut bases off the start of a read, if below
      a threshold quality (default 3)
   -  TRAILING quality (int): Cut bases off the end of a read, if below
      a threshold quality (default 0)
   -  MINLEN (int): Drop the read if it is below a specified length
      (default 60)
   -  AVGQUAL (int): Drop the read if the average quality is below the
      specified level (default 0)

.. raw:: html

   -  <s> (there is another parameter MAXINFO:: , but I'm not sure we will
      use it, as I've read it was implemented incorrectly; will check if
      newer version is fixed...) </s> <em>won't use</em>

-  can take in fastq or fastq.gz

-  if the user chooses adapter trimming as well, then this step is run
   on the cutadapt_trimmed_seqs, otherwise on the input fastq files

-  we would like to capture **stderr** from this command to a log file,
   so we can run multiqc at the end

-  we won't output the trimlog for space reasons

-  Paired end - run for each pair of files

   -  ``-baseout SampleID.trimmed.fastq.gz`` is format for output files
      which end up being:

      -  SampleID.trimmed_1P.fastq.gz - for paired forward reads -
         *would be used as input to FLASH2*
      -  SampleID.trimmed_2P.fastq.gz - for paired reverse reads -
         *would be used as input to FLASH2*
      -  SampleID.trimmed_1U.fastq.gz - for unpaired forward reads
      -  SampleID.trimmed_2U.fq.gz - for unpaired reverse reads

   .. code:: bash

      ## I don't think trimmomatic will make the output directory
      mkdir ./outputs/qtrimmed_seqs

      java -jar trimmomatic-0.39.jar PE -threads numberofthreads -phred33 R1.fastq.gz R2.fastq.gz -baseout ./outputs/qtrimmed_seqs/SampleID.trimmed.fastq.gz LEADING:leading TRAILING:trailing SLIDINGWINDOW:windowsize:requiredquality MINLEN:minlen AVGQUAL:avgqual

      ## cat ./outputs/cutadapt_trimmed_seqs/MANIFEST | tail -n +2 | paste -d "," - - | parallel -C "," -j1 trimmomatic PE -threads 30 -phred33 -trimlog ./outputs/qtrimmed_seqs/{1}.log ./outputs/cutadapt_trimmed_seqs/{2} ./outputs/cutadapt_trimmed_seqs/{5} -baseout ./outputs/qtrimmed_seqs/{1}.trimmed.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:12 MINLEN:60 AVGQUAL:2 2>./outputs/multiqc_input/trimmomatic.log

-  Single end - run for each file

   .. code:: bash

      java -jar trimmomatic-0.39.jar SE -threads numberofthreads -phred33 -trimlog logfile.txt R1.fastq.gz SampleID.trimmed.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:60

4. FLASH2 merging pairs
~~~~~~~~~~~~~~~~~~~~~~~

-  **only for paired end data!!**
-  options copied from FLASH2 help

   -  **min-overlap:** (int >= 1) The minimum required overlap length
      between two reads to provide a confident overlap. Default 10bp.
   -  **max-overlap:** (int > 1) Maximum overlap length expected in
      approximately 90% of read pairs. Overlaps longer than the maximum
      overlap parameter are still considered as good overlaps, but the
      mismatch density (explained below) is calculated over the first
      max_overlap bases in the overlapped region rather than the entire
      overlap. Default 300bp
   -  **min-overlap-outie:** (int > 0) The minimum required overlap
      length between two reads to provide a confident overlap in an
      outie scenario. Default 35bp.
   -  **max-mismatch-density:** (0 <= float <= 1) Maximum allowed ratio
      between the number of mismatched base pairs and the overlap
      length. Default 0.25.

   .. code:: bash

      ## produce compressed output, as it's the last step
      ## for each pair of fastq files:
      flash2 -D -z -m min_overlap -t numberofthreads --max-mismatch-density max_mismatch_density -o SampleID_merged -d ./outputs/merged R1.fastq.gz R2.fastq.gz

      ## ls ./outputs/qtrimmed_seqs/*1P.fastq.gz | sed 's/1P\.fastq\.gz//' | parallel -j1 flash2 -D -z -t 30 -o {/} -M 200 -x 0.5 -d ./outputs/merged {}1P.fastq.gz {}2P.fastq.gz

5. multiQC
~~~~~~~~~~

-  This step will run at the end

-  multiqc takes as argument the directory with the fastqc and other log
   output and provides summary report

.. code:: bash

   multiqc -o outputs -c multiqc_config.yaml ./outputs/multiqc_input

Cleanup
~~~~~~~

-  I think we should delete intermediate files; otherwise the download
   will potentially be 3-4x size of the input files
-  If only fastqc is run, then we don't need to throw anything away
-  If any of the other steps is run, we should only keep the final
   fastq.gz files. All others should be discarded (even ones from
   previous steps)
-  Discard all \*.qza files (if not holding in memory in Artifact API)
