Biobakery Pipeline Documentation
================================

Usage
-----

.. code:: bash

     usage: biobakery.py [-h] --job_id JOB_ID --map_file MAP_FILE
			 [--data_type {WGS_SE,WGS_PE}] [--threads THREADS]
			 [--file_ext {fastq.gz,fastq,fq.gz,fq,fasta,fasta.gz}]
			 [--local_jobs LOCAL_JOBS] [--strainphlan]
			 [--project_name PROJECT_NAME]

    optional arguments:
      -h, --help            show this help message and exit
      --job_id JOB_ID       nephele job id (default: None)
      --map_file MAP_FILE   full path to mapping file (default: None)
      --data_type {WGS_SE,WGS_PE}
			    input sequence data type (default: WGS_SE)
      --threads THREADS     number of processors per task (local job).
			    threads*local_jobs <= 8. If none, will use total
			    available processors/local_jobs. (default: None)
      --file_ext {fastq.gz,fastq,fq.gz,fq,fasta,fasta.gz}
			    the input sequence file extension (default: fastq)
      --local_jobs LOCAL_JOBS
			    number of bb tasks to run at a time. (default: 4)
      --strainphlan         user option. run strainphlan. (default: False)
      --project_name PROJECT_NAME
			    user option. project name for visualization pipeline.
			    if none, will use job_id. (default: None)

-  By default the pipeline runs in single end mode and does not run
   `StrainPhlAn <http://segatalab.cibio.unitn.it/tools/strainphlan/>`__
   .
-  User documentation [ `html for Nephele website <https://github.niaid.nih.gov/bcbb/nephele2/tree/master/pipelines/bioBakery/biobakerywgs_pipeline.html>`__ ][
   :doc:`sphinx <pipelines.biobakery.user_doc>` ] - contains language for the strainphlan
   user option

Paired end
~~~~~~~~~~

-  pass ``--data_type WGS_PE``
-  The mapping file is only used if ``--data_type WGS_PE`` is passed.
   Then, the first 3 columns must be:

   1. #SampleID
   2. ForwardFastqFile
   3. ReverseFastqFile

-  The SampleID can contain alphanumeric characters, dashes, and
   underscores. **No periods.**
-  Test data set
   s3://nephele-test-data/poorani_biobakery_datasets/paired_end.tar.gz

Single End
~~~~~~~~~~

-  The mapping file isn’t used, but we require anyway for consistency
   and validation of the upload. The first 2 columns must be:

   1. #SampleID
   2. ForwardFastqFile

-  The SampleID can contain alphanumeric characters, dashes, and
   underscores. **No periods.**
-  Test data set
   s3://nephele-test-data/poorani_biobakery_datasets/single_end.tar.gz

Installation
------------

-  Development was done using the `biobakery workflows docker
   image <https://bitbucket.org/biobakery/biobakery/wiki/biobakery_workflows#rst-header-install-with-docker>`__

   -  There are some notes of how Poorani set things up
      `here <https://github.niaid.nih.gov/bcbb/nephele2/blob/master/pipelines/bioBakery/bbdocker.md>`__

-  The docker image environment can be replicated by installing the
   `biobakery workflows with
   Homebrew <https://github.com/biobakery/homebrew-biobakery>`__ using
   all default settings.

-  Databases will need to be downloaded from inside the docker container
   or instance where biobakery is installed:

   .. code:: bash

      sudo /home/linuxbrew/.linuxbrew/bin/biobakery_workflows_databases --install wmgx

   -  They get downloaded to /opt - that’s where the workflows expect
      them to be. It is possible to specify an alternative database
      location using ``--location`` parameter and
      ``$BIOBAKERY_WORKFLOWS_DATABASES`` env variable, see `Install
      databases <https://bitbucket.org/biobakery/biobakery_workflows/wiki/Home#!installation>`__
   -  The databases are quite large uncompressed ~30G

