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

process prottrans_embs {
    publishDir params.release_dir, mode: 'copy'

    input:
    path prot_trans_original
    path sorted_ids

    output:
    path "emb.prottrans.parquet", emit: emb_prottrans

    script:
    """
    plm_create_prottrans_embs.py ${prot_trans_original} ${sorted_ids} emb.prottrans.parquet
    """
}

process calc_ankh {
    //conda 'conda_envs/ankh_wsl.yml'
    publishDir params.release_dir, mode: 'copy'
    label 'long'

    input:
    path sorted_uniprot_not_large
    path all_uniprot_ids
    path ankh_cache_path
    path previous_embs_dir

    output:
    path "emb.ankh_base.parquet", emit: emb_ankh_base
    path "emb.ankh_large.parquet", emit: emb_ankh_large

    script:
    """
    ls -la ./
    plm_ankh_calc.py ${sorted_uniprot_not_large} ${ankh_cache_path} ${all_uniprot_ids} . ${previous_embs_dir}
    """
}
