#!/usr/bin/env nextflow

evi_not_to_use_path = projectDir+"/evi_not_to_use.txt"
go_basic_path = projectDir+"/databases/go-basic.obo"
goa_raw_path = projectDir+"/databases/goa_uniprot_all.gaf.gz"

process download_uniprot{
    //publishDir "${params.release_dir}/databases", mode: 'copy'

    output:
    path "uniprot_sprot.fasta.gz", emit: uniprot_fasta

    script:
    """
    wget https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/uniprot_sprot.fasta.gz
    """
}

process sort_uniprot{
    publishDir "${params.release_dir}/databases", mode: 'copy'
    
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

process process_goa{
    conda 'env.yaml'

    input:
        path evi_not_to_use
        path go_basic
        path uniprot_fasta
        path goa_raw

    output:
        path "goa_parsed_expanded.mf.tsv.gz"
        path "goa_parsed_expanded.bp.tsv.gz"
        path "goa_parsed_expanded.cc.tsv.gz"

    shell:
    """
    python $projectDir/src/download_annotation.py
    """
}

workflow {
    download_uniprot()
    sort_uniprot(download_uniprot.out.uniprot_fasta)
    //process_goa(evi_not_to_use_path, go_basic_path, sort_uniprot.out.uniprot_sorted, goa_raw_path)
}