nextflow.enable.dsl = 2

include {
    create_caches ;
    download_gocheck_do_not_annotate ;
    download_esm ;
    download_uniprot ;
    download_taxallnomy ;
    filter_large_proteins ;
    list_taxids_test ;
    list_taxids_train ;
    process_train_terms ;
    index_go_by_term ;
    taxa_profiles_train ;
    taxa_profiles_test ;
    calc_ankh_v1 ;
    calc_esm2 ;
    process_cafa_annotations ;
    copy_additional_files
} from './modules/cafa_processes.nf'

process download_go {
    input:
    val url

    output:
    path "go-basic.obo", emit: go_basic

    script:
    """
    wget ${url}
    """
}

process download_prot5 {
    input:
    val url

    output:
    path "per-protein.h5", emit: prot5_embs_h5

    script:
    """
    wget ${url}
    """
}

process download_goa {
    input:
    val url

    output:
    path "goa_uniprot_all.gaf.gz", emit: go_annotation_raw

    script:
    """
    wget ${url}
    """
}

process join_old_releases {
    input:
    path protein_ids_path
    path "releases/*"
    path src_dir

    output:
    path "joins_dir/emb.prottrans.parquet"
    path "joins_dir/emb.ankh_large.parquet"
    path "joins_dir/emb.ankh_base.parquet"
    path "joins_dir/emb.esm2_t30.parquet"
    path "joins_dir/emb.esm2_t33.parquet"
    path "joins_dir/emb.esm2_t36.parquet"
    path "joins_dir", emit: joined_dfs_output_dir

    script:
    """
    python -u ${src_dir}/map_releases.py ${protein_ids_path} joins_dir releases/*
    """
}

process split_fasta {
    input:
    path fasta_path
    val max_tokens
    path src_dir

    output:
    path "uniprot.*.fasta", emit: fasta_files

    script:
    """
    python ${src_dir}/split_fasta.py ${fasta_path} ${max_tokens} uniprot
    """
}

process consult_interpro {
    input:
    path fasta_file
    path src_dir

    output:
    path "${fasta_file.baseName}_interpro.tsv", emit: interpro_tsv

    script:
    """
    python ${src_dir}/request_interpros.py ${fasta_file} ${fasta_file.baseName}_interpro.tsv
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
    maxForks 1

    input:
    path input_fasta
    val interpro_data_dir

    output:
    path "interpro_out/*.tsv", emit: interpro_tsvs

    script:
    """
    # 1. Capture the absolute path of the input fasta and current work dir
    # We must do this before changing directories so we don't lose the files!
    FASTA_ABS_PATH=\$(readlink -f ${input_fasta})
    OUTER_WORKDIR=\$(pwd)
    
    # 2. Define and move to a permanent execution directory
    PERSISTENT_DIR="${interpro_data_dir}/interproscan_persistent_run"
    mkdir -p \$PERSISTENT_DIR
    cd \$PERSISTENT_DIR

    
    # Isolate the inner pipeline's work directory to prevent conflicts
    export NXF_WORK=${interpro_data_dir}/interproscan_nxf_work
    export NXF_SINGULARITY_CACHEDIR="${projectDir}/singularity/sif"
    
    # Run the EBI pipeline. It will automatically chunk your Swissprot fasta 
    # based on the --batch-size parameter
    nextflow run ebi-pf-team/interproscan6 -resume \\
        -r 6.0.0 \\
        -profile singularity \\
        --formats tsv \\
        --datadir ${interpro_data_dir} \\
        --input \$FASTA_ABS_PATH \\
        --cpus ${params.interproscan_cpus} \\
        --maxWorkers 5 \\
        --outdir interpro_out \\
        --outprefix swissprot_interpro \\
        --batch-size 5000
    # 5. Return to the outer pipeline's working directory
    cd \$OUTER_WORKDIR
    mkdir -p interpro_out
    
    # 6. Copy the outputs back so the outer pipeline can emit them
    cp \$PERSISTENT_DIR/interpro_out/*.tsv interpro_out/
    """
}



process train_interpro_autoencoder {
    input:
    path interpro_tsvs
    path src_dir

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
    python -u ${src_dir}/interpro_autoencoder.py 512 model_512 ${interpro_tsvs} > model_512/stdout.log 2> model_512/stderr.log
    python -u ${src_dir}/interpro_autoencoder.py 640 model_640 ${interpro_tsvs} > model_640/stdout.log 2> model_640/stderr.log
    python -u ${src_dir}/interpro_autoencoder.py 896 model_896 ${interpro_tsvs} > model_896/stdout.log 2> model_896/stderr.log
    python -u ${src_dir}/interpro_autoencoder.py 1280 model_1280 ${interpro_tsvs} > model_1280/stdout.log 2> model_1280/stderr.log
    """
}

process make_interpro_onehotencoder {
    input:
    path interpro_tsvs
    path src_dir

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
    python -u ${src_dir}/interpro_onehot.py 800 model_800 ${interpro_tsvs} > model_800/stdout.log 2> model_800/stderr.log
    python -u ${src_dir}/interpro_onehot.py 1600 model_1600 ${interpro_tsvs} > model_1600/stdout.log 2> model_1600/stderr.log
    python -u ${src_dir}/interpro_onehot.py 3200 model_3200 ${interpro_tsvs} > model_3200/stdout.log 2> model_3200/stderr.log
    python -u ${src_dir}/interpro_onehot.py 6400 model_6400 ${interpro_tsvs} > model_6400/stdout.log 2> model_6400/stderr.log
    """
}


params.mode = "release"
//params.esm_script_path = "esm/scripts/extract.py"
params.go_basic_url = "https://purl.obolibrary.org/obo/go/go-basic.obo"
params.esm_git_url = "https://github.com/facebookresearch/esm.git"
params.gocheck_url = "https://current.geneontology.org/ontology/subsets/gocheck_do_not_annotate.json"
params.taxallnomy_tsv_url = "https://huggingface.co/datasets/pitagoras-alves/taxallnomy/resolve/main/taxallnomy.tsv.gz"
params.max_protein_len = 1800
params.max_tokens_for_interproscan = 200000
params.evi_not_use_path = projectDir + '/evi_not_to_use.txt'
params.others_dir = projectDir + '/others'
params.old_release_paths_str = ""

params.create_taxon_profiles = false
params.create_plm_embeddings = false
params.create_esm_embeddings = true
params.create_ankh_embeddings = true
params.create_fast_embeddings = true
//params.basic_env_container = "singularity_images/basic_env.sif"
//params.env2_container = "singularity_images/env2.sif"

workflow {
    create_esm_embeddings = params.create_esm_embeddings
    create_ankh_embeddings = params.create_ankh_embeddings
    create_fast_embeddings = params.create_fast_embeddings
    def old_release_files = params.old_release_paths_str
        ? params.old_release_paths_str.split(',').collect { file(it) }
        : []

    if (params.create_plm_embeddings != true) {
        print('Not creating PLM embeddings')
        create_esm_embeddings = false
        create_ankh_embeddings = false
        create_fast_embeddings = false
    }

    goa_all_url = 'https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz'
    goa_test_url = "https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/test_input/goa_uniprot_all.gaf.gz"
    uniprot_all_url = "https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/uniprot_sprot.fasta.gz"
    uniprot_test_url = "https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/test_input/uniprot_sprot.fasta.gz"

    goa_url = Channel.value(goa_all_url)
    //old_release_paths_str = Channel.value(old_release_paths_str)
    uniprot_url = Channel.value(uniprot_all_url)
    if (params.mode == 'test') {
        goa_url = Channel.value(goa_test_url)
        uniprot_url = Channel.value(uniprot_test_url)
    }

    swissprot_path = download_uniprot(uniprot_url)
    goa_raw_path = download_goa(goa_url)
    go_basic_path = download_go(params.go_basic_url)
    gocheck_do_not_annotate = download_gocheck_do_not_annotate(params.gocheck_url)
    taxallnomy_tsv_path = download_taxallnomy(params.taxallnomy_tsv_url)


    src_dir = file(projectDir + '/src')
    split_fasta(swissprot_path, params.max_tokens_for_interproscan, src_dir)
    ch_split_fastas = split_fasta.out.fasta_files.flatten()
    consult_interpro(ch_split_fastas, src_dir)
    ch_all_tsvs = consult_interpro.out.interpro_tsv.collect()
    join_interpro_consults(ch_all_tsvs)

    run_interproscan_pipeline(
        ch_split_fastas,
        params.interpro_data_dir,
    )
    // Filter Train

    filter_large_proteins(swissprot_path, params.max_protein_len, "swissprot", src_dir)

    '''join_old_releases(
        filter_large_proteins.out.ids,
        old_release_files,
        src_dir,
    )'''

    '''if (create_ankh_embeddings || create_esm_embeddings || create_fast_embeddings) {
        train_interpro_autoencoder(
            run_interproscan_pipeline.out.interpro_tsvs,
            src_dir,
        )
        parent_dir = file(params.release_dir).getParent()
        caches_tp = create_caches(parent_dir)
        if (create_ankh_embeddings) {
            calc_ankh_v1(
                filter_large_proteins.out.fasta,
                filter_large_proteins.out.ids,
                create_caches.out.ankh_cache,
                src_dir,
                join_old_releases.out.joined_dfs_output_dir,
            )
        }

        if (create_esm_embeddings) {
            esm_dir = download_esm(params.esm_git_url)

            calc_esm2(
                filter_large_proteins.out.fasta,
                filter_large_proteins.out.ids,
                create_caches.out.fairesm_cache,
                esm_dir,
                params.others_dir,
                src_dir,
                join_old_releases.out.joined_dfs_output_dir,
            )
        }
    }

    train_terms = Channel.fromPath("databases/cafa6/Train/train_terms.tsv")
    go_basic = Channel.fromPath("databases/cafa6/Train/go-basic.obo")
    format_cafa_terms_script = Channel.fromPath("src/format_cafa_terms.py")
    process_cafa_annotations(format_cafa_terms_script, train_terms, go_basic, gocheck_do_not_annotate)
    list_taxids_test(filter_large_proteins_test.out.fasta, filter_large_proteins_test.out.ids)
    list_taxids_train(filter_large_proteins_train.out.fasta, filter_large_proteins_train.out.ids)

    if(params.create_taxon_profiles){
        train_taxonomy = Channel.fromPath("databases/cafa6/Train/train_taxonomy.tsv")
        
        // Train Profiles
        taxa_profiles_train(process_cafa_annotations.out.mf, train_taxonomy, taxallnomy_tsv_path, src_dir)
        
        taxa_profiles_test(
            list_taxids_test.out.taxids, 
            taxallnomy_tsv_path, 
            src_dir,
            taxa_profiles_train.out.top_taxa_256,
            taxa_profiles_train.out.top_taxa_128
        )

    }

    ia_tsv = Channel.fromPath("databases/cafa6/IA.tsv")
    copy_additional_files(go_basic, ia_tsv)'''
}
