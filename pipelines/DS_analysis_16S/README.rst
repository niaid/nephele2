DS Analysis Pipeline README
===========================

.. contents:: :local:

Dependency
----------
`qiime2 2018.6 <https://docs.qiime2.org/2018.6/install/>`__ - See `qiime2 backup </docs/source/qiime2_backup.rst>`__ (:doc:`sphinx <qiime2_backup>`) for how to install/restore this particular version.

Spec
----
Mariam's `Qiime 2.0 Core Diversity spec <Qiime_2.0_Core_Diversity.md>`_ (:doc:`sphinx <Qiime_2.0_Core_Diversity>`)

Dev setup
---------

`downstream_16S.sh <downstream_16S.sh>`_  is used to run `downstream_16S.py <downstream_16S.py>`_ in QIIME 2 anaconda environment.

Script usage::

  downstream_16S.sh [-h] [--job_id JOB_ID] [-i INPUTS_DIR][-o OUTPUTS_DIR]
                    [--data_type DATA_TYPE] [-a] --biom_fp BIOM_FP
                    --map_file MAP_FILE --sampling_depth SAMPLING_DEPTH


required arguments:
  --biom_fp BIOM_FP     required
  --map_file MAP_FILE   required
  --sampling_depth SAMPLING_DEPTH
                        required

optional arguments:
  -h, --help            show this help message and exit
  --job_id JOB_ID       job_id when running in Nephele. either job_id or
                        inputs_dir should be specified.
  -i INPUTS_DIR, --inputs_dir INPUTS_DIR
                        input directory for running outside Nephele. either
                        job_id or inputs_dir should be specified.
  -o OUTPUTS_DIR, --outputs_dir OUTPUTS_DIR
                        output directory for running outside Nephele. optional
  --data_type DATA_TYPE
                        Ignore - placeholder. (default: DS_Analysis)
  -a, --alpha_group_sig
                        run alpha group significance. (default: False)

**Example**:

.. code:: bash

   source ../neph2-envs/dev/environment_vars

   awssume utils/launch_ec2.py -e ../neph2-envs/dev/dev_outputs.yaml -a ami-0baf677a63aee3e2e -t c5.4xlarge -k philip_bcbb

   utils/deploy/deploy_nephele2.py --sync_worker_to_local -d 10.100.1.13

   Creating EC2...
   ec2.Instance(id='i-0a9cfc90b3cfe2f19') has been created.
   To connect type:
   ssh 10.100.1.13

   ssh 10.100.1.13

  /usr/local/src/nephele2/pipelines/DS_analysis_16S/downstream_16S.sh --map_file /mnt/EFS/user_uploads/philip/inputs/Mapping_file_corrected.txt.no_gz --inputs_dir /mnt/EFS/user_uploads/philip/inputs --biom_fp /mnt/EFS/user_uploads/philip/inputs/mcc.0.03.biom --sampling_depth 1000


User doc
--------
- `github <https://github.com/niaid/nephele2/tree/master/pipelines/DS_analysis_16S/user_doc.md>`__,
- `html <https://github.com/niaid/nephele2/tree/master/pipelines/DS_analysis_16S/amplicon_da_pipeline.html>`__
