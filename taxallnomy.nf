params.taxallnomy_git = "https://github.com/tetsufmbio/taxallnomy.git"

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
    publishDir "./", mode: 'copy'

    input:
        val taxallnomy_lin_path
    
    output:
        path "taxallnomy.tsv.gz", emit: taxallnomy_tsv_path

    script:
    """
    gzip -c $taxallnomy_lin_path > taxallnomy.tsv.gz
    """
}

workflow {
    taxallnomy_dir = clone_taxallnomy(params.taxallnomy_git)
    taxallnomy_lin = run_taxallnomy(taxallnomy_dir)
    compress_and_save_taxallnomy(taxallnomy_lin_path)
}
