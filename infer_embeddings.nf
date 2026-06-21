nextflow.enable.dsl = 2


params.mode = "release"
//params.esm_script_path = "esm/scripts/extract.py"
params.go_basic_url = "https://purl.obolibrary.org/obo/go/go-basic.obo"
params.esm_git_url = "https://github.com/facebookresearch/esm.git"
params.gocheck_url = "https://current.geneontology.org/ontology/subsets/gocheck_do_not_annotate.json"
params.taxallnomy_tsv_url = "https://huggingface.co/datasets/pitagoras-alves/taxallnomy/resolve/main/taxallnomy.tsv.gz"
params.max_protein_len = 1800
params.evi_not_use_path = projectDir + '/evi_not_to_use.txt'
params.others_dir = projectDir + '/others'

params.create_taxon_profiles = true
params.create_plm_embeddings = true
params.create_esm_embeddings = true
params.create_ankh_embeddings = true
params.create_autoencoder_embeddings = true
//params.basic_env_container = "singularity_images/basic_env.sif"
//params.env2_container = "singularity_images/env2.sif"

process infer_taxid_autoencoder {
    publishDir params.release_dir, mode: 'copy'

    input:
    path model_path
    path taxid_tsv_path

    output:
    path "emb.taxid_autoencoded.parquet", emit: emb_pq

    script:
    """
    taxo_encode.py ${model_path} uniprot_id lineage ${taxid_tsv_path}
    """
}

process infer_interpro_autoencoder {
    publishDir params.release_dir, mode: 'copy'

    input:
    path model_path
    path interpro_annot_path

    output:
    path "emb.interpro_autoencoded.parquet", emit: emb_pq

    script:
    """
    interpro_encode.py ${model_path} ids terms ${interpro_annot_path}
    """
}

workflow {
    create_esm_embeddings = params.create_esm_embeddings
    create_ankh_embeddings = params.create_ankh_embeddings
    create_autoencoder_embeddings = params.create_autoencoder_embeddings

    taxid_autoencoder_model_path = params.release_dir + "/taxid_autoencoder"
    interpro_autoencoder_model_path = params.release_dir + "/interpro_autoencoder"
    taxid_tsv_path = params.release_dir + "/taxid.tsv"
    interpro_tsv_path = params.release_dir + "/interpro.tsv"

    if (params.create_plm_embeddings != true) {
        print('Not creating PLM embeddings')
        create_esm_embeddings = false
        create_ankh_embeddings = false
    }

    print(params)

    if (create_autoencoder_embeddings) {
        infer_taxid_autoencoder(
            taxid_autoencoder_model_path,
            taxid_tsv_path,
        )
        infer_interpro_autoencoder(
            interpro_autoencoder_model_path,
            interpro_tsv_path,
        )
    }
}
