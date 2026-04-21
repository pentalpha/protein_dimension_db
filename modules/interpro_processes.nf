process consult_interpro {
    input:
    path fasta_file

    output:
    path "${fasta_file.baseName}_interpro.tsv", emit: interpro_tsv

    script:
    """
    interpro_request.py ${fasta_file} ${fasta_file.baseName}_interpro.tsv
    """
}

process join_interpro_consults {
    input:
    path tsv_files

    output:
    path "concatenated.tsv", emit: concatenated_tsv

    script:
    """
    cat ${tsv_files} > concatenated.tsv
    """
}

//#TODO: Make run_interproscan_pipeline succeed always, even when it fails
//#TODO: Join run_interproscan_pipeline results and convert to same tsv format as consult_interpro
//#TODO: Pass both results to the train_interpro_autoencoder process

process run_interproscan_pipeline {
    errorStrategy 'retry'
    maxRetries 3
    maxForks 1

    input:
    path input_fasta
    val interpro_data_dir
    val interproscan_tmp_dir

    output:
    path "interpro_out/${input_fasta.baseName}_interpro.tsv", emit: interpro_tsv

    script:
    """
    # 1. Capture the absolute path of the input fasta and current work dir
    # We must do this before changing directories so we don't lose the files!
    FASTA_ABS_PATH=\$(readlink -f ${input_fasta})
    OUTER_WORKDIR=\$(pwd)
    
    # 2. Define and move to a permanent execution directory
    PERSISTENT_DIR="${interproscan_tmp_dir}/interproscan_persistent_run"
    mkdir -p \$PERSISTENT_DIR
    cd \$PERSISTENT_DIR

    
    # Isolate the inner pipeline's work directory to prevent conflicts
    export NXF_WORK=${interproscan_tmp_dir}/interproscan_nxf_work
    export NXF_SINGULARITY_CACHEDIR="${projectDir}/singularity/sif"
    
    # Run the EBI pipeline. It will automatically chunk your Swissprot fasta 
    # based on the --batch-size parameter
    EXIT_CODE=0
    nextflow run ebi-pf-team/interproscan6 -resume \\
        -r 6.0.0 \\
        -profile singularity \\
        --formats tsv \\
        --datadir ${interpro_data_dir} \\
        --input \$FASTA_ABS_PATH \\
        --cpus ${params.interproscan_cpus} \\
        --maxWorkers 5 \\
        --outdir interpro_out \\
        --outprefix ${input_fasta.baseName}_interpro\\
        --batch-size 5000 || EXIT_CODE=\$?

    # If the pipeline failed, generate the empty file
    if [ \$EXIT_CODE -ne 0 ]; then
        echo "WARNING: InterProScan failed for ${input_fasta.baseName}. Creating an empty file to continue."
        mkdir -p interpro_out
        touch interpro_out/${input_fasta.baseName}_interpro.tsv
    fi

    # 5. Return to the outer pipeline's working directory
    cd \$OUTER_WORKDIR
    mkdir -p interpro_out
    
    # 6. Copy the outputs back so the outer pipeline can emit them
    cp \$PERSISTENT_DIR/interpro_out/*.tsv interpro_out/
    """
}

process parse_interpro_raw {
    input:
    path interpro_tsvs

    output:
    path "interpro_parsed.tsv", emit: interpro_parsed_tsv

    script:
    """
    interpro_parse.py interpro_parsed.tsv ${interpro_tsvs}
    """
}

process join_interpro_tsvs {
    storeDir "${params.raw_data_dir}/interpro_all_tsv"

    input:
    path tsv_a
    path tsv_b

    output:
    path "interpro_parsed.tsv", emit: interpro_parsed_tsv

    script:
    """
    cat ${tsv_a} ${tsv_b} > interpro_parsed.tsv
    """
}

process make_interpro_obo {
    storeDir "${params.raw_data_dir}/interpro_obo"

    input:
    path interpro_parsed_tsv

    output:
    path "interpro.obo", emit: interpro_obo

    script:
    """
    wget https://ftp.ebi.ac.uk/pub/databases/interpro/releases/latest/ParentChildTreeFile.txt -O interpro_tree.txt
    interpro_tree_to_obo.py interpro_tree.txt ${interpro_parsed_tsv} interpro.obo
    """
}

process calc_interpro_ia {
    storeDir "${params.raw_data_dir}/interpro_ia"

    input:
    path interpro_obo
    path interpro_parsed_tsv

    output:
    path "IA.txt", emit: interpro_ia_tsv

    script:
    """
    interpro_ann_to_ia_format.py ${interpro_parsed_tsv} ia_format_ann.tsv
    git clone https://github.com/pentalpha/InformationAccretion-Interpro.git
    python InformationAccretion-Interpro/ia.py --annot ia_format_ann.tsv --graph ${interpro_obo} --prop
    """
}

process make_interpro_vocab {
    input:
    path interpro_parsed_tsv
    path interpro_ia_tsv

    output:
    path "interpro_vocab_ia.tsv", emit: interpro_vocab_ia_tests_tsv
    path "interpro_vocab_ia.json", emit: interpro_vocab_ia_json
    path "interpro_vocab_top_k.tsv", emit: interpro_vocab_top_k_tests_tsv
    path "interpro_vocab_top_k.json", emit: interpro_vocab_top_k_json

    script:
    """
    interpro_make_vocab.py interpro interpro_vocab_ia 64 41800 0.95 ia_rich ${interpro_parsed_tsv} ${interpro_ia_tsv}
    interpro_make_vocab.py interpro interpro_vocab_top_k 64 41800 0.95 top_k ${interpro_parsed_tsv} ${interpro_ia_tsv}
    """
}

process train_interpro_autoencoder {
    input:
    path interpro_tsv_consult
    path sorted_vocab_json

    output:
    path "model_512", emit: model_512_dir
    path "model_640", emit: model_640_dir
    path "model_896", emit: model_896_dir
    path "model_1280", emit: model_1280_dir

    script:
    """
    mkdir -p model_512
    mkdir -p model_640
    mkdir -p model_896
    mkdir -p model_1280
    interpro_make_model.py autoencoder 512 model_512 ${interpro_tsv_consult} ${sorted_vocab_json} > model_512/stdout.log 2> model_512/stderr.log
    interpro_make_model.py autoencoder 640 model_640 ${interpro_tsv_consult} ${sorted_vocab_json} > model_640/stdout.log 2> model_640/stderr.log
    interpro_make_model.py autoencoder 896 model_896 ${interpro_tsv_consult} ${sorted_vocab_json} > model_896/stdout.log 2> model_896/stderr.log
    interpro_make_model.py autoencoder 1280 model_1280 ${interpro_tsv_consult} ${sorted_vocab_json} > model_1280/stdout.log 2> model_1280/stderr.log
    """
}

process make_interpro_onehotencoder {
    input:
    path interpro_tsv_consult
    path sorted_vocab_json

    output:
    path "model_800", emit: model_800_dir
    path "model_1600", emit: model_1600_dir
    path "model_3200", emit: model_3200_dir
    path "model_6400", emit: model_6400_dir

    script:
    """
    mkdir -p model_800
    mkdir -p model_1600
    mkdir -p model_3200
    mkdir -p model_6400
    interpro_onehot.py 800 model_800 ${interpro_tsv_consult} ${sorted_vocab_json} > model_800/stdout.log 2> model_800/stderr.log
    interpro_onehot.py 1600 model_1600 ${interpro_tsv_consult} ${sorted_vocab_json} > model_1600/stdout.log 2> model_1600/stderr.log
    interpro_onehot.py 3200 model_3200 ${interpro_tsv_consult} ${sorted_vocab_json} > model_3200/stdout.log 2> model_3200/stderr.log
    interpro_onehot.py 6400 model_6400 ${interpro_tsv_consult} ${sorted_vocab_json} > model_6400/stdout.log 2> model_6400/stderr.log
    """
}
