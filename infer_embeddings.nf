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

process create_caches {
    input:
    path parent_dir

    output:
    path "${parent_dir}/fastplms_cache", emit: fastplms_cache

    script:
    """
    mkdir -p ${parent_dir}/fastplms_cache
    """
}

process infer_taxid_autoencoder {
    label 'pytorchcpu'
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
    label 'pytorchcpu'
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

process infer_ankh_base {
    label 'pytorchgpu_light'
    publishDir params.release_dir, mode: 'copy'

    input:
    path fasta_seqs
    path cache_dir

    output:
    path "emb.ankh_base_mean.parquet", emit: emb_mean_pq
    path "emb.ankh_base_norm.parquet", emit: emb_norm_pq
    path "emb.ankh_base_max.parquet", emit: emb_max_pq
    path "emb.ankh_base_softmax.parquet", emit: emb_softmax_pq
    path "emb.ankh_base_std.parquet", emit: emb_std_pq
    path "emb.ankh_base_k4p_max.parquet", emit: emb_k4p_max_pq
    path "emb.ankh_base_k8p_max.parquet", emit: emb_k8p_max_pq
    path "emb.ankh_base_k16p_max.parquet", emit: emb_k16p_max_pq
    path "emb.ankh_base_k32p_max.parquet", emit: emb_k32p_max_pq
    path "emb.ankh_base_parti.parquet", emit: emb_parti_pq

    script:
    """
    fasta_encode.py ${fasta_seqs} ${cache_dir} ElnaggarLab/ankh-base emb.ankh_base.parquet mean,max,softmax,std,norm,k4p_max,k8p_max,k16p_max,k32p_max,parti
    """
}

process infer_ankh_large {
    label 'pytorchgpu_light'
    publishDir params.release_dir, mode: 'copy'

    input:
    path fasta_seqs
    path cache_dir

    output:
    path "emb.ankh_large_mean.parquet", emit: emb_mean_pq
    path "emb.ankh_large_max.parquet", emit: emb_max_pq
    path "emb.ankh_large_std.parquet", emit: emb_std_pq
    path "emb.ankh_large_parti.parquet", emit: emb_parti_pq

    script:
    """
    fasta_encode.py ${fasta_seqs} ${cache_dir} ElnaggarLab/ankh-large emb.ankh_large.parquet mean,max,std,parti
    """
}

process infer_ankh2_large {
    label 'pytorch251'
    publishDir params.release_dir, mode: 'copy'

    input:
    path fasta_seqs
    path cache_dir

    output:
    path "emb.ankh2_large_mean.parquet", emit: emb_mean_pq
    path "emb.ankh2_large_max.parquet", emit: emb_max_pq
    path "emb.ankh2_large_std.parquet", emit: emb_std_pq
    path "emb.ankh2_large_parti.parquet", emit: emb_parti_pq

    script:
    """
    fasta_encode.py ${fasta_seqs} ${cache_dir} ElnaggarLab/ankh2-ext2 emb.ankh2_large.parquet
    """
}

process infer_ankh3_large {
    label 'pytorch251'
    publishDir params.release_dir, mode: 'copy'

    input:
    path fasta_seqs
    path cache_dir

    output:
    path "emb.ankh3_large_mean.parquet", emit: emb_mean_pq
    path "emb.ankh3_large_max.parquet", emit: emb_max_pq
    path "emb.ankh3_large_std.parquet", emit: emb_std_pq
    path "emb.ankh3_large_parti.parquet", emit: emb_parti_pq

    script:
    """
    fasta_encode.py ${fasta_seqs} ${cache_dir} ElnaggarLab/ankh3-large emb.ankh3_large.parquet
    """
}

process infer_pfe1_300 {
    label 'pytorchgpu_light'
    publishDir params.release_dir, mode: 'copy'

    input:
    path fasta_seqs
    path cache_dir

    output:
    path "emb.e1_300_mean.parquet", emit: emb_mean_pq
    path "emb.e1_300_max.parquet", emit: emb_max_pq
    path "emb.e1_300_std.parquet", emit: emb_std_pq
    path "emb.e1_300_parti.parquet", emit: emb_parti_pq

    script:
    """
    fasta_encode.py ${fasta_seqs} ${cache_dir} Profluent-Bio/E1-300m emb.e1_300.parquet
    """
}

process infer_pfe1_600 {
    label 'pytorchgpu_light'
    publishDir params.release_dir, mode: 'copy'

    input:
    path fasta_seqs
    path cache_dir

    output:
    path "emb.e1_600_mean.parquet", emit: emb_mean_pq
    path "emb.e1_600_max.parquet", emit: emb_max_pq
    path "emb.e1_600_std.parquet", emit: emb_std_pq
    path "emb.e1_600_parti.parquet", emit: emb_parti_pq

    script:
    """
    fasta_encode.py ${fasta_seqs} ${cache_dir} Profluent-Bio/E1-600m emb.e1_600.parquet
    """
}

//TODO:
//Synthyra/ESM2-650M
//Synthyra/ESM2-3B
//Synthyra/ESMplusplus_small
//Synthyra/ESMplusplus_large
//Synthyra/ESM3_small
//DPLM2-650M
//DPLM2-3B
//Ankh3-XL?

workflow {
    create_esm_embeddings = params.create_esm_embeddings
    create_ankh_embeddings = params.create_ankh_embeddings
    create_autoencoder_embeddings = params.create_autoencoder_embeddings

    taxid_autoencoder_model_path = params.release_dir + "/taxid_autoencoder"
    interpro_autoencoder_model_path = params.release_dir + "/interpro_autoencoder"
    taxid_tsv_path = params.release_dir + "/taxid.tsv"
    interpro_tsv_path = params.release_dir + "/interpro.tsv"
    swissprot_fasta = params.release_dir + "/sequences.swissprot.fasta"

    parent_dir = file(params.release_dir).getParent()
    create_caches(parent_dir)

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

    if (create_ankh_embeddings) {
        infer_ankh_base(
            swissprot_fasta,
            create_caches.out.fastplms_cache,
        )
        infer_ankh_large(
            swissprot_fasta,
            create_caches.out.fastplms_cache,
        )
        infer_ankh2_large(
            swissprot_fasta,
            create_caches.out.fastplms_cache,
        )
        infer_ankh3_large(
            swissprot_fasta,
            create_caches.out.fastplms_cache,
        )
    }

    infer_pfe1_300(
        swissprot_fasta,
        create_caches.out.fastplms_cache,
    )
    infer_pfe1_600(
        swissprot_fasta,
        create_caches.out.fastplms_cache,
    )
}
