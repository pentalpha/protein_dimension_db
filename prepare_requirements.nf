esm_script_path = "esm/scripts/extract.py"

process download_go {
    publishDir "databases", mode: 'copy'
    
    output:
    path "go-basic.obo"

    script:
    """
    wget https://purl.obolibrary.org/obo/go/go-basic.obo
    """
}

process download_gocheck_do_not_annotate {
    publishDir "databases", mode: 'copy'
    
    output:
    path "gocheck_do_not_annotate.json"

    script:
    """
    wget https://current.geneontology.org/ontology/subsets/gocheck_do_not_annotate.json
    """
}

process download_esm{
    publishDir "./", mode: 'copy'

    output:
    path "esm"
    
    script:
    """
    git clone git@github.com:facebookresearch/esm.git
    """
}

workflow {
    download_esm()
    download_gocheck_do_not_annotate()
    download_go()
}