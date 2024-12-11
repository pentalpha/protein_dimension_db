#!/usr/bin/env nextflow

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

workflow {
    download_uniprot()
    //download_uniprot.out.uniprot_fasta.view()
    sort_uniprot(download_uniprot.out.uniprot_fasta)
    sort_uniprot.out.ids.view()
    sort_uniprot.out.uniprot_sorted.view()
}