#!/usr/bin/env nextflow

evi_not_to_use_path = projectDir+"/evi_not_to_use.txt"
go_basic_path = projectDir+"/databases/go-basic.obo"
goa_raw_path = projectDir+"/databases/goa_uniprot_all.gaf.gz"
uniprot_path = "${projectDir}/databases/uniprot_sprot.fasta.gz"
prot_trans_path = "${projectDir}/databases/per-protein.h5"
taxallnomy_tsv_path = "${projectDir}/databases/taxallnomy.tsv.gz"
esm_cache_path = "${params.release_dir}/../fairesm_cache"
ankh_cache_path = "${params.release_dir}/../ankh_caches"
max_protein_len = 1800

process sort_uniprot{
    publishDir params.release_dir, mode: 'copy'
    
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
    
    output:
        path "uniprot_sorted.not_large.fasta", emit: uniprot_not_large

    script:
    """
    python $projectDir/src/filter_fasta_by_len.py $sorted_uniprot uniprot_sorted.not_large.fasta $max_protein_len
    """
}

process calc_esm_embeddings{
    conda 'conda_envs/pytorch2.yml'
    publishDir params.release_dir, mode: 'copy'
    
    input:
        path sorted_uniprot_not_large
        path all_uniprot_ids
    
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

process calc_ankh_embeddings{
    conda 'conda_envs/ankh.yml'
    publishDir params.release_dir, mode: 'copy'
    
    input:
        path sorted_uniprot_not_large
        path all_uniprot_ids
    
    output:
        path "emb.ankh_large.parquet"
        path "emb.ankh_base.parquet"

    script:
    """
    python $projectDir/src/ankh_calc.py $sorted_uniprot_not_large $ankh_cache_path $all_uniprot_ids
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

process prottrans_embs{
    conda 'conda_envs/env2.txt'
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

process process_goa{
    conda 'conda_envs/env2.txt'
    publishDir params.release_dir, mode: 'copy'

    input:
        path sorted_ids

    output:
        path "go.experimental.tsv.gz"
        path "go.expanded.tsv.gz"
        path "go.experimental.mf.tsv.gz", emit: go_experimental_mf
        path "go.experimental.bp.tsv.gz"
        path "go.experimental.cc.tsv.gz"

    shell:
    """
    python $projectDir/src/filter_gaf.py $projectDir $goa_raw_path $sorted_ids ./ rerun
    """
}

process taxa_profiles{
    conda 'conda_envs/basic_env.txt'
    publishDir params.release_dir, mode: 'copy'

    input:
        path go_experimental_mf
        path taxids_path

    output:
        path "onehot.taxa_256.npy.gz"
        path "emb.taxa_profile_256.npy.gz"
        path "top_taxa_256.txt"
        path "onehot.taxa_128.npy.gz"
        path "emb.taxa_profile_128.npy.gz"
        path "top_taxa_128.txt"

    shell:
    """
    python $projectDir/src/calc_taxa_profiles.py $taxallnomy_tsv_path $go_experimental_mf $taxids_path
    """
}

workflow {
    sort_uniprot(uniprot_path)
    filter_large_proteins(sort_uniprot.out.uniprot_sorted)
    list_taxids(uniprot_path, sort_uniprot.out.ids)
    process_goa(sort_uniprot.out.ids)
    

    prottrans_embs(prot_trans_path, sort_uniprot.out.ids)
    calc_ankh_embeddings(filter_large_proteins.out.uniprot_not_large, sort_uniprot.out.ids)
    calc_esm_embeddings(filter_large_proteins.out.uniprot_not_large, sort_uniprot.out.ids)
    taxa_profiles(process_goa.out.go_experimental_mf, list_taxids.out.taxids)
}