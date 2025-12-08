
params.mode = "release"
//params.esm_script_path = "esm/scripts/extract.py"
params.go_basic_url = "https://purl.obolibrary.org/obo/go/go-basic.obo"
params.esm_git_url = "https://github.com/facebookresearch/esm.git"
params.gocheck_url = "https://current.geneontology.org/ontology/subsets/gocheck_do_not_annotate.json"
params.taxallnomy_tsv_url = "https://huggingface.co/datasets/pitagoras-alves/taxallnomy/resolve/main/taxallnomy.tsv.gz"
params.max_protein_len = 1800
params.evi_not_use_path = projectDir+'/evi_not_to_use.txt'
params.others_dir = projectDir+'/others'

params.create_taxon_profiles = true
params.create_plm_embeddings = true
params.create_esm_embeddings = true
params.create_ankh_embeddings = true
//params.esm_cache_path = "${params.release_dir}/../fairesm_cache"
//params.ankh_cache_path = "${params.release_dir}/../ankh_caches"

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
        val output_name
    
    output:
        path "${output_name}.fasta", emit: fasta
        path "${output_name}_ids.txt", emit: ids

    script:
    """
    python $projectDir/src/filter_fasta_by_len.py $input_fastas ${output_name}.fasta ${output_name}_ids.txt $max_protein_len
    """
}

process list_taxids{
    publishDir params.release_dir, mode: 'copy'
    
    input:
        path swissprot_fasta
        path trembl_fasta
        path sorted_ids
    
    output:
        path "taxid.tsv", emit: taxids

    script:
    """
    python $projectDir/src/list_uniprot_taxids.py $swissprot_fasta $trembl_fasta $sorted_ids taxid.tsv
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
        path "onehot.taxa_256.parquet"
        path "emb.taxa_profile_256.parquet"
        path "top_taxa_256.txt", emit: top_taxa_256
        path "onehot.taxa_128.parquet"
        path "emb.taxa_profile_128.parquet"
        path "top_taxa_128.txt", emit: top_taxa_128

    script:
    //singularity exec --bind $projectDir/src:/src --bind $taxallnomy_tsv_path:/$taxallnomy_tsv_path --bind $go_experimental_mf:/$go_experimental_mf --bind $taxids_path:/$taxids_path \\
    //$projectDir/$params.basic_env_container \\ 
    """
    python3 src/calc_taxa_profiles.py $taxallnomy_tsv_path $go_experimental_mf $taxids_path
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
        path "onehot.taxa_256_test.parquet"
        path "emb.taxa_profile_256_test.parquet"
        path "onehot.taxa_128_test.parquet"
        path "emb.taxa_profile_128_test.parquet"

    script:
    """
    # Create top_taxa_dir and move files there so script finds them
    mkdir top_taxa_dir
    cp $top_taxa_256 top_taxa_dir/top_taxa_256.txt
    cp $top_taxa_128 top_taxa_dir/top_taxa_128.txt
    
    # We pass a dummy file for go_experimental_mf since it won't be used (we provide top taxa)
    touch dummy_go.tsv
    
    python3 src/calc_taxa_profiles.py $taxallnomy_tsv_path dummy_go.tsv $taxids_path top_taxa_dir "_test"
    """
}

process list_taxids_simple{
    publishDir params.release_dir, mode: 'copy'
    
    input:
        path input_fasta
        path ids_list
    
    output:
        path "taxid_test.tsv", emit: taxids

    script:
    """
    # We pass the same fasta twice because the script expects 2 inputs before ids output etc
    # Usage: python list_uniprot_taxids.py <fasta1> <fasta2> <ids> <output>
    # Actually wait, the script takes fasta1, fasta2. I can pass a dummy empty file for fasta2.
    touch dummy.fasta
    python $projectDir/src/list_uniprot_taxids.py $input_fasta dummy.fasta $ids_list taxid_test.tsv
    """
}

process calc_ankh_embeddings{
    //conda 'conda_envs/ankh_wsl.yml'
    publishDir params.release_dir, mode: 'copy'
    
    input:
        path sorted_uniprot_not_large
        path all_uniprot_ids
        path ankh_cache_path
        val output_suffix
    
    output:
        path "emb.ankh_base${output_suffix}.parquet"

    script:
    """
    python $projectDir/src/ankh_calc.py $sorted_uniprot_not_large $ankh_cache_path $all_uniprot_ids $output_suffix
    """
}

process calc_esm_embeddings{
    //conda 'conda_envs/pytorch2.yml'
    publishDir params.release_dir, mode: 'copy'
    
    input:
        path sorted_uniprot_not_large
        path all_uniprot_ids
        path esm_cache_path
        path esm_dir
        path others_dir
        val output_suffix
    
    output:
        path "emb.esm2_t33${output_suffix}.parquet"

    script:
    """
    python $projectDir/src/esm_calc.py $sorted_uniprot_not_large $esm_cache_path $all_uniprot_ids $output_suffix
    """
}

workflow {
    create_esm_embeddings = params.create_esm_embeddings
    create_ankh_embeddings = params.create_ankh_embeddings
    create_prottrans_embeddings = params.create_prottrans_embeddings

    if(params.create_plm_embeddings != true){
        print('Not creating PLM embeddings')
        create_esm_embeddings = false
        create_ankh_embeddings = false
        create_prottrans_embeddings = false
    }

    print(params)
    
    train_fasta = Channel.fromPath("databases/cafa6/Train/train_sequences.fasta")
    test_fasta = Channel.fromPath("databases/cafa6/Test/testsuperset.fasta")
    
    // Filter Train
    filter_large_proteins_train = filter_large_proteins.alias('filter_large_proteins_train')
    filter_large_proteins_train(train_fasta, params.max_protein_len, "train_filtered")
    
    // Filter Test
    filter_large_proteins_test = filter_large_proteins.alias('filter_large_proteins_test')
    filter_large_proteins_test(test_fasta, params.max_protein_len, "test_filtered")
    
    
    if(create_ankh_embeddings || create_esm_embeddings){
        parent_dir = file(params.release_dir).getParent()
        caches_tp = create_caches(parent_dir)
        
        if(create_ankh_embeddings){
            // Train
            calc_ankh_embeddings_train = calc_ankh_embeddings.alias('calc_ankh_embeddings_train')
            calc_ankh_embeddings_train(
                filter_large_proteins_train.out.fasta, 
                filter_large_proteins_train.out.ids, 
                create_caches.out.ankh_cache,
                "_train"
            )
            // Test
            calc_ankh_embeddings_test = calc_ankh_embeddings.alias('calc_ankh_embeddings_test')
            calc_ankh_embeddings_test(
                filter_large_proteins_test.out.fasta, 
                filter_large_proteins_test.out.ids, 
                create_caches.out.ankh_cache,
                "_test"
            )
        }
        
        if(create_esm_embeddings){
            esm_dir = download_esm(params.esm_git_url)
            
            // Train
            calc_esm_embeddings_train = calc_esm_embeddings.alias('calc_esm_embeddings_train')
            calc_esm_embeddings_train(
                filter_large_proteins_train.out.fasta, 
                filter_large_proteins_train.out.ids, 
                create_caches.out.fairesm_cache, 
                esm_dir, 
                params.others_dir,
                "_train"
            )
            
            // Test
            calc_esm_embeddings_test = calc_esm_embeddings.alias('calc_esm_embeddings_test')
            calc_esm_embeddings_test(
                filter_large_proteins_test.out.fasta, 
                filter_large_proteins_test.out.ids, 
                create_caches.out.fairesm_cache, 
                esm_dir, 
                params.others_dir,
                "_test"
            )
        }
    }

    
    taxallnomy_tsv_path = download_taxallnomy(params.taxallnomy_tsv_url)
    
    src_dir = file(projectDir+'/src')
    if(params.create_taxon_profiles){
        train_terms = Channel.fromPath("databases/cafa6/Train/train_terms.tsv")
        train_taxonomy = Channel.fromPath("databases/cafa6/Train/train_taxonomy.tsv")
        
        // Train Profiles
        taxa_profiles_train(train_terms, train_taxonomy, taxallnomy_tsv_path, src_dir)
        
        // Test Profiles
        // Need to extract taxids for test set first
        list_taxids_simple(filter_large_proteins_test.out.fasta, filter_large_proteins_test.out.ids)
        
        taxa_profiles_test(
            list_taxids_simple.out.taxids, 
            taxallnomy_tsv_path, 
            src_dir,
            taxa_profiles_train.out.top_taxa_256,
            taxa_profiles_train.out.top_taxa_128
        )
    }
    //release_dir_channel = Channel.fromPath(params.release_dir)
}