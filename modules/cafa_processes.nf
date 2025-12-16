process create_caches {
    input:
    path parent_dir

    output:
    path "${parent_dir}/ankh_caches", emit: ankh_cache
    path "${parent_dir}/fairesm_cache", emit: fairesm_cache

    script:
    """
    mkdir -p ${parent_dir}/ankh_caches
    mkdir -p ${parent_dir}/fairesm_cache
    """
}

process process_cafa_annotations {
    publishDir params.release_dir, mode: 'copy'
    
    input:
        path format_cafa_script
        path train_terms
        path go_basic
        path go_not_use
    
    output:
        path "go.experimental.mf.tsv", emit: mf
        path "go.experimental.bp.tsv", emit: bp
        path "go.experimental.cc.tsv", emit: cc

    script:
    """
    ls ./
    echo $projectDir
    python $format_cafa_script $train_terms go.experimental $go_not_use $go_basic
    """
}

process download_gocheck_do_not_annotate {
    //publishDir "databases", mode: 'copy'
    
    input:
    val url

    output:
    path "gocheck_do_not_annotate.json", emit: gocheck_do_not_annotate

    script:
    """
    wget $url
    """
}

process download_esm{
    //publishDir "libs/", mode: 'copy'
    
    input:
    val esm_git

    output:
    path "esm", emit: esm_dir
    
    script:
    """
    git clone $esm_git
    """
}

process download_uniprot{
    //publishDir "databases", mode: 'copy'

    input:
    val url

    output:
    path "uniprot_sprot.fasta.gz", emit: uniprot_fasta

    script:
    """
    wget $url
    """
}


process download_taxallnomy{
    //publishDir "libs/", mode: 'copy'

    input:
        val url
    
    output:
        path "taxallnomy.tsv.gz", emit: taxallnomy_tsv_path

    script:
    """
    wget $url
    """
}

process sort_train_fasta{
    publishDir params.release_dir, mode: 'copy'
    
    input:
        path original_uniprot
    
    output:
        path "swissprot.txt", emit: ids
        path "swissprot_sorted.fasta.gz", emit: uniprot_sorted

    script:
    """
    python $projectDir/src/sort_uniprot.py $original_uniprot swissprot.txt swissprot_sorted.fasta.gz
    """
}

process filter_large_proteins{
    publishDir params.release_dir, mode: 'copy'
    
    input:
        path input_fastas
        val max_protein_len
        val output_suffix
    
    output:
        path "sequences.${output_suffix}.fasta", emit: fasta
        path "ids.${output_suffix}.txt", emit: ids

    script:
    """
    python $projectDir/src/filter_fasta_by_len.py $input_fastas sequences.${output_suffix}.fasta ids.${output_suffix}.txt $max_protein_len
    """
}

process list_taxids_train{
    publishDir params.release_dir, mode: 'copy'
    
    input:
        path swissprot_fasta
        path sorted_ids
    
    output:
        path "taxid.train.tsv", emit: taxids

    script:
    """
    touch dummy.fasta
    python $projectDir/src/list_uniprot_taxids.py $swissprot_fasta dummy.fasta $sorted_ids taxid.train.tsv
    """
}

process process_train_terms{
    //conda 'conda_envs/env2_wsl.txt'
    publishDir params.release_dir, mode: 'copy'

    input:
        //path sorted_ids
        //path goa_raw_path
        //path gocheck_do_not_annotate
        //path evi_not_to_use_path
        //path go_basic_path

    output:
        path "go.experimental.parquet"
        path "go.experimental_expanded.parquet"
        path "go.by_uniprot.parquet", emit: go_by_uniprot
        //#
        //path "go.experimental.mf.tsv.gz", emit: go_experimental_mf
        //#path "go.experimental.bp.tsv.gz"
        //path "go.experimental.cc.tsv.gz"

    script:
    """
    python $projectDir/src/filter_cafa_terms.py $evi_not_to_use_path $gocheck_do_not_annotate $go_basic_path $goa_raw_path $sorted_ids rerun
    """
}

process index_go_by_term{
    //conda 'conda_envs/env2_wsl.txt'
    publishDir params.release_dir, mode: 'copy'

    input:
        path go_by_uniprot

    output:
        path "go.by_term.parquet", emit: go_by_term

    script:
    """
    python $projectDir/src/ann_by_term.py $go_by_uniprot
    """
}

process taxa_profiles_train{
    publishDir params.release_dir, mode: 'copy'

    input:
        path go_experimental_mf
        path taxids_path
        path taxallnomy_tsv_path
        path src_dir

    output:
        path "onehot.taxa_128.train.parquet"
        path "onehot.taxa_256.train.parquet"
        path "emb.taxa_profile_128.train.parquet"
        path "emb.taxa_profile_256.train.parquet"
        path "top_taxa_128.txt", emit: top_taxa_128
        path "top_taxa_256.txt", emit: top_taxa_256

    script:
    //singularity exec --bind $projectDir/src:/src --bind $taxallnomy_tsv_path:/$taxallnomy_tsv_path --bind $go_experimental_mf:/$go_experimental_mf --bind $taxids_path:/$taxids_path \\
    //$projectDir/$params.basic_env_container \\ 
    """
    python3 src/calc_taxa_profiles.py $taxallnomy_tsv_path $go_experimental_mf $taxids_path train
    """
}

process taxa_profiles_test{
    publishDir params.release_dir, mode: 'copy'

    input:
        path taxids_path
        path taxallnomy_tsv_path
        path src_dir
        path top_taxa_256
        path top_taxa_128

    output:
        path "emb.taxa_profile_128.test.parquet"
        path "emb.taxa_profile_256.test.parquet"
        path "onehot.taxa_128.test.parquet"
        path "onehot.taxa_256.test.parquet"

    script:
    """
    # Create top_taxa_dir and move files there so script finds them
    mkdir top_taxa_dir
    cp $top_taxa_256 top_taxa_dir/top_taxa_256.txt
    cp $top_taxa_128 top_taxa_dir/top_taxa_128.txt
    
    # We pass a dummy file for go_experimental_mf since it won't be used (we provide top taxa)
    touch dummy_go.tsv
    
    python3 src/calc_taxa_profiles.py $taxallnomy_tsv_path dummy_go.tsv $taxids_path test top_taxa_dir
    """
}

process list_taxids_test{
    publishDir params.release_dir, mode: 'copy'
    
    input:
        path input_fasta
        path ids_list
    
    output:
        path "taxid.test.tsv", emit: taxids

    script:
    """
    touch dummy.fasta
    python $projectDir/src/list_uniprot_taxids.py $input_fasta dummy.fasta $ids_list taxid.test.tsv
    """
}

process calc_ankh_embeddings{
    //conda 'conda_envs/ankh_wsl.yml'
    publishDir params.release_dir, mode: 'copy'
    label 'long'
    
    input:
        path sorted_uniprot_not_large
        path all_uniprot_ids
        path ankh_cache_path
        path src_dir
        val output_suffix
    
    output:
        path "emb.ankh_base.${output_suffix}.parquet"

    script:
    """
    ls -la ./
    python $src_dir/ankh_calc.py $sorted_uniprot_not_large $ankh_cache_path $all_uniprot_ids $output_suffix
    """
}

process calc_esm_embeddings{
    //conda 'conda_envs/pytorch2.yml'
    publishDir params.release_dir, mode: 'copy'
    label 'long'
    
    input:
        path sorted_uniprot_not_large
        path all_uniprot_ids
        path esm_cache_path
        path esm_dir
        path others_dir
        path src_dir
        val output_suffix
    
    output:
        path "emb.esm2_t33.${output_suffix}.parquet"

    script:
    """
    python $src_dir/esm_calc.py $sorted_uniprot_not_large $esm_cache_path $all_uniprot_ids $output_suffix
    """
}

process copy_additional_files {
    publishDir params.release_dir, mode: 'copy'

    input:
    path go_basic, stageAs: 'go-basic.source.obo'
    path ia_tsv

    output:
    path "go-basic.obo"
    path "information_acc.tsv"

    script:
    """
    cp $go_basic go-basic.obo
    cp $ia_tsv information_acc.tsv
    """
}
