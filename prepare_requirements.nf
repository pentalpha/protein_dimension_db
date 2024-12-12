esm_script_path = "esm/scripts/extract.py"
go_basic_url = "https://purl.obolibrary.org/obo/go/go-basic.obo"
esm_git_url = "git@github.com:facebookresearch/esm.git"
gocheck_url = "https://current.geneontology.org/ontology/subsets/gocheck_do_not_annotate.json"
goa_all_url = 'https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz'
uniprot_url = "https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/uniprot_sprot.fasta.gz"
//ProtT5 swiss_prot
prot_t5_embs_url = "https://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/embeddings/uniprot_sprot/per-protein.h5"

process download_go {
    publishDir "databases", mode: 'copy'
    
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
    publishDir "databases", mode: 'copy'
    
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
    publishDir "libs/", mode: 'copy'
    
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
    publishDir "databases", mode: 'copy'

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
    publishDir "databases", mode: 'copy'

    input:
    val url

    output:
    path "per-protein.h5", emit: prot5_embs_h5

    script:
    """
    wget $url
    """
}

process download_uniprot{
    publishDir "databases", mode: 'copy'

    output:
    path "uniprot_sprot.fasta.gz", emit: uniprot_fasta

    script:
    """
    wget https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/uniprot_sprot.fasta.gz
    """
}

process download_goa{
    publishDir "databases", mode: 'copy'

    input:
        val url
    
    output:
        path "goa_uniprot_all.gaf.gz", emit: go_annotation_raw

    script:
    """
    wget $url
    """
}

workflow {
    download_uniprot(uniprot_url)
    download_uniprot()
    //download_goa(goa_all_url)
    download_go(go_basic_url)
    download_esm(esm_git_url)
    download_gocheck_do_not_annotate(gocheck_url)
    download_prot5(prot_t5_embs_url)
}