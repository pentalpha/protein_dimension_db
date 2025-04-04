
params.mode = "release"
//params.esm_script_path = "esm/scripts/extract.py"
params.go_basic_url = "https://purl.obolibrary.org/obo/go/go-basic.obo"
params.esm_git_url = "https://github.com/facebookresearch/esm.git"
params.taxallnomy_git_url = "https://github.com/tetsufmbio/taxallnomy.git"
params.gocheck_url = "https://current.geneontology.org/ontology/subsets/gocheck_do_not_annotate.json"
//ProtT5 swiss_prot
params.prot_t5_embs_url = "https://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/embeddings/uniprot_sprot/per-protein.h5"
params.max_protein_len = 1800
params.evi_not_use_path = projectDir+'/evi_not_to_use.txt'
params.others_dir = projectDir+'/others'
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

process download_go {
    //publishDir "databases", mode: 'copy'
    
    input:
    val url
    
    output:
    path "go-basic.obo", emit: go_basic

    script:
    """
    wget $url
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

process download_prot5{
    //publishDir "databases", mode: 'copy'

    input:
    val url

    output:
    path "per-protein.h5", emit: prot5_embs_h5

    script:
    """
    wget $url
    """
}

process download_goa{
    //publishDir "databases", mode: 'copy'

    input:
        val url
    
    output:
        path "goa_uniprot_all.gaf.gz", emit: go_annotation_raw

    script:
    """
    wget $url
    """
}

process clone_taxallnomy{
    //publishDir "libs/", mode: 'copy'

    input:
        val taxallnomy_git
    
    output:
        path "taxallnomy", emit: taxallnomy_dir

    script:
    """
    git clone $taxallnomy_git
    """
}

process run_taxallnomy{
    input:
        val taxallnomy_dir
    
    output:
        path "taxallnomy_data/taxallnomy_lin.tab", emit: taxallnomy_lin

    script:
    """
    perl $taxallnomy_dir/generate_taxallnomy.pl
    """
}

process compress_and_save_taxallnomy{
    //publishDir "databases", mode: 'copy'

    input:
        val taxallnomy_lin_path
    
    output:
        path "taxallnomy.tsv.gz", emit: taxallnomy_tsv_path

    script:
    """
    gzip -c $taxallnomy_lin_path > taxallnomy.tsv.gz
    """
}

process sort_uniprot{
    //publishDir params.release_dir, mode: 'copy'
    
    input:
        path original_uniprot
    
    output:
        path "ids.txt", emit: ids
        path "uniprot_sorted.fasta.gz", emit: uniprot_sorted

    script:
    """
    python $projectDir/src/sort_uniprot.py $original_uniprot ids.txt uniprot_sorted.fasta.gz
    """
}

process filter_large_proteins{
    publishDir params.release_dir, mode: 'copy'
    
    input:
        path sorted_uniprot
        val max_protein_len
    
    output:
        path "uniprot_sorted.not_large.fasta", emit: uniprot_not_large

    script:
    """
    python $projectDir/src/filter_fasta_by_len.py $sorted_uniprot uniprot_sorted.not_large.fasta $max_protein_len
    """
}

process list_taxids{
    publishDir params.release_dir, mode: 'copy'
    
    input:
        path original_uniprot
        path sorted_ids
    
    output:
        path "taxid.tsv", emit: taxids

    script:
    """
    python $projectDir/src/list_uniprot_taxids.py $original_uniprot $sorted_ids taxid.tsv
    """
}

process process_goa{
    conda 'conda_envs/env2_wsl.txt'
    publishDir params.release_dir, mode: 'copy'

    input:
        path sorted_ids
        path goa_raw_path
        path gocheck_do_not_annotate
        path evi_not_to_use_path
        path go_basic_path

    output:
        path "go.experimental.tsv.gz"
        path "go.expanded.tsv.gz"
        path "go.experimental.mf.tsv.gz", emit: go_experimental_mf
        path "go.experimental.bp.tsv.gz"
        path "go.experimental.cc.tsv.gz"

    script:
    """
    python $projectDir/src/filter_gaf.py $evi_not_to_use_path $gocheck_do_not_annotate $go_basic_path $goa_raw_path $sorted_ids rerun
    """
}

process prottrans_embs{
    conda 'conda_envs/env2_wsl.txt'
    publishDir params.release_dir, mode: 'copy'
    
    input:
        path prot_trans_original
        path sorted_ids
    
    output:
        path "emb.prottrans.parquet", emit: emb_prottrans

    script:
    """
    python $projectDir/src/create_prottrans_embs.py $prot_trans_original $sorted_ids emb.prottrans.parquet
    """
}

process taxa_profiles{
    conda 'conda_envs/basic_env.txt'
    publishDir params.release_dir, mode: 'copy'

    input:
        path go_experimental_mf
        path taxids_path
        path taxallnomy_tsv_path

    output:
        path "onehot.taxa_256.npy.gz"
        path "emb.taxa_profile_256.npy.gz"
        path "top_taxa_256.txt"
        path "onehot.taxa_128.npy.gz"
        path "emb.taxa_profile_128.npy.gz"
        path "top_taxa_128.txt"

    script:
    """
    python $projectDir/src/calc_taxa_profiles.py $taxallnomy_tsv_path $go_experimental_mf $taxids_path
    """
}

process calc_ankh_embeddings{
    conda 'conda_envs/ankh_wsl.yml'
    publishDir params.release_dir, mode: 'copy'
    
    input:
        path sorted_uniprot_not_large
        path all_uniprot_ids
        path ankh_cache_path
    
    output:
        path "emb.ankh_large.parquet"
        path "emb.ankh_base.parquet"

    script:
    """
    python $projectDir/src/ankh_calc.py $sorted_uniprot_not_large $ankh_cache_path $all_uniprot_ids
    """
}

process calc_esm_embeddings{
    conda 'conda_envs/pytorch2.yml'
    publishDir params.release_dir, mode: 'copy'
    
    input:
        path sorted_uniprot_not_large
        path all_uniprot_ids
        path esm_cache_path
        path esm_dir
        path others_dir
    
    output:
        path "emb.esm2_t6.parquet"
        path "emb.esm2_t12.parquet"
        path "emb.esm2_t30.parquet"
        path "emb.esm2_t33.parquet"
        path "emb.esm2_t36.parquet"

    script:
    """
    python $projectDir/src/esm_calc.py $sorted_uniprot_not_large $esm_cache_path $all_uniprot_ids
    """
}

workflow {
    goa_all_url = 'https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz'
    goa_test_url = "https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/test_input/goa_uniprot_all.gaf.gz"
    uniprot_all_url = "https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/uniprot_sprot.fasta.gz"
    uniprot_test_url = "https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/test_input/uniprot_sprot.fasta.gz"
    goa_url = Channel.value(goa_all_url)
    uniprot_url = Channel.value(uniprot_all_url)
    if (params.mode == 'test'){
        goa_url = Channel.value(goa_test_url)
        uniprot_url = Channel.value(uniprot_test_url)
    }
    
    uniprot_path = download_uniprot(uniprot_url)
    goa_raw_path = download_goa(goa_url)
    go_basic_path = download_go(params.go_basic_url)
    gocheck_do_not_annotate = download_gocheck_do_not_annotate(params.gocheck_url)
    esm_dir = download_esm(params.esm_git_url)
    //esm_dir + "/scripts/extract.py"
    taxallnomy_dir = clone_taxallnomy(params.taxallnomy_git_url)
    //taxallnomy_lin = run_taxallnomy(taxallnomy_dir)
    //taxallnomy_tsv_path = compress_and_save_taxallnomy(taxallnomy_lin)

    sort_uniprot(uniprot_path)
    not_large_proteins = filter_large_proteins(sort_uniprot.out.uniprot_sorted, params.max_protein_len)
    taxids = list_taxids(uniprot_path, sort_uniprot.out.ids)
    process_goa(sort_uniprot.out.ids, goa_raw_path, 
        gocheck_do_not_annotate, params.evi_not_use_path, go_basic_path)
    
    prot_trans_path = download_prot5(params.prot_t5_embs_url)
    prottrans_embs(prot_trans_path, sort_uniprot.out.ids)
    //taxa_profiles(process_goa.out.go_experimental_mf, taxids, taxallnomy_tsv_path)

    //release_dir_channel = Channel.fromPath(params.release_dir)
    parent_dir = file(params.release_dir).getParent()
    caches_tp = create_caches(parent_dir)
    //calc_ankh_embeddings(not_large_proteins, sort_uniprot.out.ids, 
    //    create_caches.out.ankh_cache)
    //calc_esm_embeddings(not_large_proteins, sort_uniprot.out.ids, 
    //    create_caches.out.fairesm_cache, esm_dir, params.others_dir)
}