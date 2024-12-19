#!/usr/bin/env nextflow

evi_not_to_use_path = projectDir+"/evi_not_to_use.txt"
go_basic_path = projectDir+"/databases/go-basic.obo"
goa_raw_path = projectDir+"/databases/goa_uniprot_all.gaf.gz"
uniprot_path = "${projectDir}/databases/uniprot_sprot.fasta.gz"
prot_trans_path = "${projectDir}/databases/per-protein.h5"
taxallnomy_tsv_path = "${projectDir}/databases/taxallnomy.tsv.gz"

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
        path "emb.prottrans.npy.gz", emit: emb_prottrans

    script:
    """
    python $projectDir/src/create_prottrans_embs.py $prot_trans_original $sorted_ids emb.prottrans.npy.gz
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
    list_taxids(uniprot_path, sort_uniprot.out.ids)
    //prottrans_embs(prot_trans_path, sort_uniprot.out.ids)
    process_goa(sort_uniprot.out.ids)
    taxa_profiles(process_goa.out.go_experimental_mf, list_taxids.out.taxids)
}