# WGSA pipeline #
More details can be found here: https://github.niaid.nih.gov/angelovaag/WmGS-Nephele/blob/master/WmGSpipeline-v1_4developers.md

- [Scripts](#scripts--files)
- [WGSA instance and perfomance testing](#wgsa-instance-and-perfomance-testing)

## Scripts & files

The pipeline is implemented using [Snakemake](https://snakemake.readthedocs.io/).

- [Snakefile](Snakefile) - main snakefile with pipeline rules and vars
- [utils.smk](utils.smk) - helper/utility functions for [inclusion](https://snakemake.readthedocs.io/en/stable/tutorial/additional_features.html?highlight=include#modularization) in Snakefile
- [config.yaml](config.yaml) - config paths to scripts/tools
- [WGSA_output_diagram.pdf](WGSA_output_diagram.pdf) - output files diagram for user
- [wgsa.py](wgsa.py) - pipeline script
  - **required:**  
      `--job_id JOB_ID`  job_id & `--map_file MAP_FILE`   Mapping file
  - **optional:**  
      `-i INPUTS_DIR`, `-o OUTPUTS_DIR`  input/output directory for running outside Nephele.  
      `--decontaminate {human,house_mouse}` (default: human)  
      `--include_ted_files` keep trimmed/decontaminated reads fastq files  
      `--data_type DATA_TYPE` Ignore - placeholder.  
      `--stdout_redirect {file,std}` choose where to redirect stdout. Using std to see it
                            from terminal (default: file)  
      `--stderr_redirect {file,std}` choose where to redirect stderr. Using std to see it
                            from terminal (default: file)  
      
  - **trim with fastp** `--trim_filter`  
      `--average_read_quality AVERAGE_READ_QUALITY`
      `--min_read_length MIN_READ_LENGTH`
      `--trim_of_5 TRIM_OF_5`
      `--trim_of_3 TRIM_OF_3`
  
  

# WGSA instance and perfomance testing
1. We are currently using EC2 m5.12xlarge (48 vCPU & 192GB memory) for WGSA pipeline. I was thinking between m5.8x large and m5.12x large. After a couple of tests, I decided to go with m5.12x large. There were 2 set of data, 05M (140MB per sample) and 1M (280MB per sample), the results:
Please note that these tests were running without involving EFS

    - 05M dataset with 1 sample
        - 25 mins (m5 8x)
        - 22 mins (m5 12x)

    - 05M dataset with 3 samples
        - 45 mins (m5 8x)
        - 31 mins (m5 12x)

    - 1M  dataset with 3 samples 
        - 1h10 mins (m5 8x)
        - 40 mins (m5 12x)

2. WGSA pipeline was originally written by shell script and then converted to Snamake for Nephele. Regarding parallelism, Snakemake can help to distribute multiple rules (one or a set of commands) at a time. This helps to reduce the time taken to finish the pipeline when user submit a job with multiple samples. 

3. WGSA uses many tools (fastp, bbmap, metaspades, metaprokka, metabat, verse, checkm, etc) and each tool has different way to utilize the computation resource (threads/cores, memory). By default, we allocate 16 threads/cores if the tool has option for setting the threads/cores. However, since the EC2 instance we chose has 48 vCPU, we need an util function to make use all of those cores in case the job only have 1 or 2 samples. With 1 sample, we use all 48 cores. With 2 samples, we use 24 cores for each job. With more or equal than 3 samples, we use minimum 16 cores for each rule (that has the setting for threading). The number 16 threads being used here after some testings with a set of data. It provides the balance/sweet spot for performance as well as parallelism (there are some tools that might not need up to 16 threads, but for the sake of simplicity, I do not want to have too many configurations). If you are reading this note and trying to optimize the pipeline, please feel free to change the number of threads, EC2 instance type, dataset or another version of the tools. 
Note: You can think of vCPU, cores, threads in this context mean the same thing. 
Some testing examples:
    - metaspades
        - 8  threads: 2m52s
        - 16 threads: 1m58s
        - 32 threads: 1m44s
        - 48 threads: 1m47s
    - bowtie2
        - 8  threads: 60s
        - 16 threads: 30s
        - 32 threads: 30s
    - samtools
        - 4 threads:  15s
        - 8 threads:  10s
        - 16 threads: 10s
    - metaprokka
        - 8 threads: 19 mins
        - 16 threads: 14 mins
        - 32 threads: 12 mins

4. By the time I used bbmap.sh, this tool somehow needs a huge amount of memory. I need to set to run one kind of this job at a time, otherwise, it will fail with heap out of memory issue.

