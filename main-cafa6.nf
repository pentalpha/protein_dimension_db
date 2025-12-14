nextflow.enable.dsl=2

include { 
    create_caches; 
    download_gocheck_do_not_annotate;
    download_esm;
    download_uniprot;
    download_taxallnomy;
    filter_large_proteins as filter_large_proteins_train;
    filter_large_proteins as filter_large_proteins_test;
    list_taxids_test;
    list_taxids_train;
    process_train_terms;
    index_go_by_term;
    taxa_profiles_train;
    taxa_profiles_test;
    calc_ankh_embeddings as calc_ankh_embeddings_train;
    calc_ankh_embeddings as calc_ankh_embeddings_test;
    calc_esm_embeddings as calc_esm_embeddings_train;
    calc_esm_embeddings as calc_esm_embeddings_test;
    process_cafa_annotations;
} from './modules/cafa_processes.nf'

params.mode = "release"
//params.esm_script_path = "esm/scripts/extract.py"
params.go_basic_url = "https://purl.obolibrary.org/obo/go/go-basic.obo"
params.esm_git_url = "https://github.com/facebookresearch/esm.git"
params.gocheck_url = "https://current.geneontology.org/ontology/subsets/gocheck_do_not_annotate.json"
params.taxallnomy_tsv_url = "https://huggingface.co/datasets/pitagoras-alves/taxallnomy/resolve/main/taxallnomy.tsv.gz"
params.max_protein_len = 1800
params.evi_not_use_path = projectDir+'/evi_not_to_use.txt'
params.others_dir = projectDir+'/others'

params.create_taxon_profiles = true
params.create_plm_embeddings = true
params.create_esm_embeddings = true
params.create_ankh_embeddings = true
//params.basic_env_container = "singularity_images/basic_env.sif"
//params.env2_container = "singularity_images/env2.sif"

workflow {
    create_esm_embeddings = params.create_esm_embeddings
    create_ankh_embeddings = params.create_ankh_embeddings
    
    if(params.create_plm_embeddings != true){
        print('Not creating PLM embeddings')
        create_esm_embeddings = false
        create_ankh_embeddings = false
    }

    print(params)
    
    train_fasta = Channel.fromPath("databases/cafa6/Train/train_sequences.fasta")
    test_fasta = Channel.fromPath("databases/cafa6/Test/testsuperset.fasta")
    
    // Filter Train
    filter_large_proteins_train(train_fasta, params.max_protein_len, "train")
    
    // Filter Test
    filter_large_proteins_test(test_fasta, params.max_protein_len, "test")
    
    
    if(create_ankh_embeddings || create_esm_embeddings){
        parent_dir = file(params.release_dir).getParent()
        caches_tp = create_caches(parent_dir)
        
        src_dir = file(projectDir+'/src')
        if(create_ankh_embeddings){
            // Train
            calc_ankh_embeddings_train(
                filter_large_proteins_train.out.fasta, 
                filter_large_proteins_train.out.ids, 
                create_caches.out.ankh_cache,
                src_dir,
                "train"
            )
            // Test
            calc_ankh_embeddings_test(
                filter_large_proteins_test.out.fasta, 
                filter_large_proteins_test.out.ids, 
                create_caches.out.ankh_cache,
                src_dir,
                "test"
            )
        }
        
        if(create_esm_embeddings){
            esm_dir = download_esm(params.esm_git_url)
            
            // Train
            calc_esm_embeddings_train(
                filter_large_proteins_train.out.fasta, 
                filter_large_proteins_train.out.ids, 
                create_caches.out.fairesm_cache, 
                esm_dir, 
                params.others_dir,
                src_dir,
                "train"
            )
            
            // Test
            calc_esm_embeddings_test(
                filter_large_proteins_test.out.fasta, 
                filter_large_proteins_test.out.ids, 
                create_caches.out.fairesm_cache, 
                esm_dir, 
                params.others_dir,
                src_dir,
                "test"
            )
        }
    }

    taxallnomy_tsv_path = download_taxallnomy(params.taxallnomy_tsv_url)
    train_terms = Channel.fromPath("databases/cafa6/Train/train_terms.tsv")
    process_cafa_annotations(train_terms)
    list_taxids_test(filter_large_proteins_test.out.fasta, filter_large_proteins_test.out.ids)
    list_taxids_train(filter_large_proteins_train.out.fasta, filter_large_proteins_train.out.ids)

    if(params.create_taxon_profiles){
        train_taxonomy = Channel.fromPath("databases/cafa6/Train/train_taxonomy.tsv")
        
        // Train Profiles
        taxa_profiles_train(process_cafa_annotations.out.mf, train_taxonomy, taxallnomy_tsv_path, src_dir)
        
        taxa_profiles_test(
            list_taxids_test.out.taxids, 
            taxallnomy_tsv_path, 
            src_dir,
            taxa_profiles_train.out.top_taxa_256,
            taxa_profiles_train.out.top_taxa_128
        )
    }
}