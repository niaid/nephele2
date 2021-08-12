def picrust2_pipeline_cmd():
    # Generate the command to run picrust2 pipeline
    cmd = (
        "rm -rf {config[output_dir]} && {params.script_path} -s {input.fasta} -i {input.biom} -o {config[output_dir]} -p {threads} "
        "--max_nsti {params.max_nsti} --min_reads {params.min_reads} --min_samples {params.min_samples} --verbose --remove_intermediate"
    )
    if config["stratified"] == "stratified_only":
        cmd += " --stratified"
    elif config["stratified"] == "stratified_per_sequence_contrib":
        cmd += " --stratified --per_sequence_contrib"
    return cmd
