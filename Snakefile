configfile: "config.yml"

extract_esm_script = "esm/scripts/extract.py"

#quickgo_gaf = "databases/QuickGO-annotations-1697773507369-20231020.gaf"
#quickgo_expanded = "input/quickgo_expanded.tsv.gz"
#goa_parsed = 'databases/goa_parsed.tsv.gz'
#goa_parsed_expanded = 'databases/goa_parsed_expanded.tsv.gz'

go_basic = "databases/go-basic.obo"

release_dir = config['release_dir']

uniprot_fasta = release_dir+"/databases/uniprot_sprot.fasta.gz"
goa_parsed_frequent = release_dir+"/databases/goa_parsed_frequent.tsv.gz"
#proteins_for_learning = "input/proteins.fasta"
annotation_path = release_dir+'/annotation.tsv'
taxon_profile_path = release_dir+'/taxa_profile.tsv.gz'
esm_features_prefix = release_dir+'/esm2_t'
labels_path = release_dir+'/go_labels.tsv'
ids_path = release_dir+'/ids.txt'

conda_run1 = "conda run -n dimension_db --live-stream"

rule download_go:
    output:
        go_basic
    shell:
        "cd databases && wget https://current.geneontology.org/ontology/subsets/gocheck_do_not_annotate.json"
        " && wget https://purl.obolibrary.org/obo/go/go-basic.obo"

rule download_esm:
    output:
        extract_esm_script
    shell:
        "rm -rf esm && git clone git@github.com:facebookresearch/esm.git"

rule create_release_dir:
    output:
        release_dir

rule download_goa:
    input:
        'evi_not_to_use.txt',
        go_basic
    output:
        goa_parsed_frequent
    shell:
        "mkdir -p release_dir && " + conda_run1 + " python src/download_annotation.py"

rule download_uniprot:
    input:
        goa_parsed_frequent
    output:
        uniprot_fasta
    shell:
        "wget -O " + uniprot_fasta + " https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/uniprot_sprot.fasta.gz"

'''rule annotated_protein_list:
    input:
        goa_parsed_frequent,
        uniprot_fasta
    output:
        proteins_for_learning,
        input_annotation_path
    shell:
        "conda run --live-stream -n plm python src/create_train_protein_set.py"

rule create_features:
    input:
        proteins_for_learning,
        extract_esm_script
    output:
        input_features_ids_path
    shell:
        "conda run --live-stream -n plm python src/calc_features.py "+proteins_for_learning+" input/features"
    
rule create_taxon_profiles:
    input:
        input_features_ids_path
    output:
        features_taxon_profile_path
    shell:
        "conda run --live-stream -n plm python src/calc_taxon_dist.py"'''