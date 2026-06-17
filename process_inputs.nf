nextflow.enable.dsl = 2

include {
    create_caches ;
    download_gocheck_do_not_annotate ;
    download_esm ;
    download_uniprot ;
    download_taxallnomy ;
    filter_large_proteins ;
    index_go_by_term ;
    copy_additional_files
} from './modules/cafa_processes.nf'

include {
    consult_interpro ;
    join_interpro_consults ;
    run_interproscan_pipeline ;
    parse_interpro_raw ;
    join_interpro_tsvs ;
    make_interpro_obo ;
    calc_interpro_ia ;
    make_interpro_vocab
} from './modules/interpro_processes.nf'

include {
    download_prot5
} from './modules/embedders.nf'

process make_taxallnomy_parquet {
    publishDir params.release_dir, mode: 'copy'

    input:
    path taxallnomy_tsv_path

    output:
    path "taxallnomy.parquet", emit: taxallnomy_parquet

    script:
    """
    taxo_taxallnomy_parquet.py ${taxallnomy_tsv_path} taxallnomy.parquet
    """
}

process list_taxids {
    publishDir params.release_dir, mode: 'copy'

    input:
    path swissprot_fasta
    path taxallnomy_parquet
    path protein_ids_path

    output:
    path "taxid.tsv", emit: taxids

    script:
    """
    export HOME=\$PWD
    taxo_list_uniprot_taxids.py ${swissprot_fasta} ${protein_ids_path} ${taxallnomy_parquet} taxid.tsv
    """
}

process make_taxid_obo {
    publishDir params.release_dir, mode: 'copy'

    input:
    path taxid_tsv

    output:
    path "taxid.obo", emit: taxid_obo

    script:
    """
    taxo_tree_to_obo.py ${taxid_tsv} taxid.obo
    """
}

process calc_taxid_ia {
    storeDir "${params.raw_data_dir}/taxid_ia"
    publishDir params.release_dir, mode: 'copy'

    input:
    path taxid_tsv
    path taxid_obo

    output:
    path "taxid_IA.tsv", emit: taxid_ia_tsv

    script:
    """
    interpro_ann_to_ia_format.py ${taxid_tsv} taxid_ann_format.tsv ncbi_taxid
    git clone https://github.com/pentalpha/InformationAccretion-Interpro.git
    python InformationAccretion-Interpro/ia.py --outfile taxid_IA.tsv --annot taxid_ann_format.tsv --graph ${taxid_obo} --prop
    """
}

process make_taxid_vocab {
    storeDir "${params.raw_data_dir}/taxid_vocab"
    publishDir params.release_dir, mode: 'copy'

    input:
    path taxids
    path taxid_ia_tsv

    output:
    path "taxid_vocab_ia.tsv", emit: taxid_vocab_ia_tests_tsv
    path "taxid_vocab_ia.json", emit: taxid_vocab_ia_json

    script:
    """
    interpro_make_vocab.py taxallnomy taxid_vocab_ia 1600 12000 0.99 ia_rich ${taxids} ${taxid_ia_tsv}
    interpro_make_vocab.py taxallnomy taxid_vocab_ia 1600 12000 0.99 top_k ${taxids} ${taxid_ia_tsv}
    """
}

process download_go {
    storeDir "${params.raw_data_dir}/go"
    publishDir params.release_dir, mode: 'copy'

    input:
    val url

    output:
    path "go-basic.obo", emit: go_basic

    script:
    """
    wget ${url}
    """
}

process download_goa {
    storeDir "${params.raw_data_dir}/goa"

    input:
    val url

    output:
    path "goa_uniprot_all.gaf.gz", emit: go_annotation_raw

    script:
    """
    wget ${url}
    """
}



process join_old_releases {
    storeDir "${params.raw_data_dir}/joins"

    input:
    path protein_ids_path
    path "releases/*"

    output:
    path "joins_dir/emb.prottrans.parquet", emit: emb_prottrans
    path "joins_dir/emb.ankh_large.parquet", emit: emb_ankh_large
    path "joins_dir/emb.ankh_base.parquet", emit: emb_ankh_base
    path "joins_dir/emb.esm2_t30.parquet", emit: emb_esm2_t30
    path "joins_dir/emb.esm2_t33.parquet", emit: emb_esm2_t33
    path "joins_dir/emb.esm2_t36.parquet", emit: emb_esm2_t36
    path "joins_dir", emit: joined_dfs_output_dir

    script:
    """
    mkdir -p joins_dir
    plm_map_releases.py ${protein_ids_path} joins_dir releases/*
    """
}

process split_fasta {
    storeDir "${params.raw_data_dir}/splitted_fastas"

    input:
    path fasta_path
    val max_tokens

    output:
    path "uniprot.*.fasta", emit: fasta_files

    script:
    """
    fasta_split.py ${fasta_path} ${max_tokens} uniprot
    """
}




params.mode = "release"
//params.esm_script_path = "esm/scripts/extract.py"
params.go_basic_url = "https://purl.obolibrary.org/obo/go/go-basic.obo"
params.esm_git_url = "https://github.com/facebookresearch/esm.git"
params.gocheck_url = "https://current.geneontology.org/ontology/subsets/gocheck_do_not_annotate.json"
params.taxallnomy_tsv_url = "https://huggingface.co/datasets/pitagoras-alves/taxallnomy/resolve/main/taxallnomy.tsv.gz"
params.prot_t5_embs_url = "https://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/embeddings/uniprot_sprot/per-protein.h5"

params.max_protein_len = 1800
params.max_tokens_for_interproscan = 200000
params.evi_not_use_path = projectDir + '/evi_not_to_use.txt'
params.others_dir = projectDir + '/others'
params.old_release_paths_str = ""

params.create_taxon_profiles = false
params.create_plm_embeddings = false
params.create_esm_embeddings = true
params.create_ankh_embeddings = true
params.create_fast_embeddings = true
params.create_prottrans_embeddings = true
//params.basic_env_container = "singularity_images/basic_env.sif"
//params.env2_container = "singularity_images/env2.sif"

workflow {
    create_esm_embeddings = params.create_esm_embeddings
    create_ankh_embeddings = params.create_ankh_embeddings
    create_fast_embeddings = params.create_fast_embeddings
    create_prottrans_embeddings = params.create_prottrans_embeddings

    joins_dir_path = params.release_dir + "/raw_data/joins/joins_dir"

    def old_release_files = params.old_release_paths_str
        ? params.old_release_paths_str.split(',').collect { file(it) }
        : []

    if (params.create_plm_embeddings != true) {
        print('Not creating PLM embeddings')
        create_esm_embeddings = false
        create_ankh_embeddings = false
        create_fast_embeddings = false
    }

    goa_all_url = 'https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz'
    goa_test_url = "https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/test_input/goa_uniprot_all.gaf.gz"
    uniprot_all_url = "https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/uniprot_sprot.fasta.gz"
    uniprot_test_url = "https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/test_input/uniprot_sprot.fasta.gz"

    goa_url = Channel.value(goa_all_url)
    //old_release_paths_str = Channel.value(old_release_paths_str)
    uniprot_url = Channel.value(uniprot_all_url)
    if (params.mode == 'test') {
        goa_url = Channel.value(goa_test_url)
        uniprot_url = Channel.value(uniprot_test_url)
    }

    swissprot_path = download_uniprot(uniprot_url)
    goa_raw_path = download_goa(goa_url)
    go_basic_path = download_go(params.go_basic_url)
    gocheck_do_not_annotate = download_gocheck_do_not_annotate(params.gocheck_url)
    taxallnomy_tsv_path = download_taxallnomy(params.taxallnomy_tsv_url)

    make_taxallnomy_parquet(taxallnomy_tsv_path)

    src_dir = file(projectDir + '/src')
    split_fasta(swissprot_path, params.max_tokens_for_interproscan)
    ch_split_fastas = split_fasta.out.fasta_files.flatten()
    consult_interpro(ch_split_fastas)
    ch_all_tsvs = consult_interpro.out.interpro_tsv.collect()
    join_interpro_consults(ch_all_tsvs)

    make_interpro_obo(join_interpro_consults.out.concatenated_tsv)
    calc_interpro_ia(
        make_interpro_obo.out.interpro_obo,
        join_interpro_consults.out.concatenated_tsv,
    )

    make_interpro_vocab(
        join_interpro_consults.out.concatenated_tsv,
        calc_interpro_ia.out.interpro_ia_tsv,
    )

    // Filter Train

    filter_large_proteins(swissprot_path, params.max_protein_len, "swissprot")

    list_taxids(
        swissprot_path,
        make_taxallnomy_parquet.out.taxallnomy_parquet,
        filter_large_proteins.out.ids,
    )

    make_taxid_obo(list_taxids.out.taxids)
    calc_taxid_ia(
        list_taxids.out.taxids,
        make_taxid_obo.out.taxid_obo,
    )

    make_taxid_vocab(
        list_taxids.out.taxids,
        calc_taxid_ia.out.taxid_ia_tsv,
    )
}
