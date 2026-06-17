nextflow.enable.dsl = 2

process train_taxid_autoencoder {
    storeDir "${params.raw_data_dir}/train_taxid_autoencoder"
    publishDir params.release_dir, mode: 'copy', pattern: 'taxid_autoencoder'

    input:
    path taxid_obo_path
    path taxid_tsv_path
    path taxall_pq_path
    path taxid_vocab_path

    output:
    path "taxid_autoencoder", emit: model_dir

    script:
    """
    taxo_make_model.py taxid.tsv taxid_vocab_ia.json 75 6e-3
    mv model_final taxid_autoencoder
    """
}

process train_interpro_autoencoder {
    storeDir "${params.raw_data_dir}/train_interpro_autoencoder"
    publishDir params.release_dir, mode: 'copy', pattern: 'interpro_autoencoder'

    input:
    path join_interpro_consults_path
    path interpro_vocab_ia_json_path
    path interpro_obo_path

    output:
    path "interpro_autoencoder", emit: model_dir

    script:
    """
    interpro_make_model.py concatenated.tsv interpro_vocab_ia.json 75 6e-3
    mv model_final interpro_autoencoder
    """
}

params.mode = "release"

workflow {
    interpro_obo_path = params.release_dir + "/interpro.obo"
    join_interpro_consults_path = params.release_dir + "/raw_data/interpro_consults/concatenated.tsv"
    interpro_vocab_ia_json_path = params.release_dir + "/interpro_vocab_ia.json"

    taxid_obo_path = params.release_dir + "/taxid.obo"
    taxid_tsv_path = params.release_dir + "/taxid.tsv"
    taxall_pq_path = params.release_dir + "/taxallnomy.parquet"
    taxid_vocab_path = params.release_dir + "/taxid_vocab_ia.json"

    train_taxid_autoencoder(
        taxid_obo_path,
        taxid_tsv_path,
        taxall_pq_path,
        taxid_vocab_path
    )
}
