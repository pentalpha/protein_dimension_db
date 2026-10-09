import os
import sys
import polars as pl
from datasets import Dataset
from huggingface_hub import DatasetCard, DatasetCardData


basecard = r"""---
{{ card_data }}
---
# Protein Dimension DB 2.0 - Functional Annotation Datasets - [ont_name]

Gene Ontology (MF, BP, CC) and DeepLoc 2.1 annotations for all proteins in the UniprotKB/SwissProt database.
This is a parsed version of the Uniprot-GOA dataset, filtered for only SwissProt proteins and with annotations grouped by evidence type.

## Gene Ontology annotations

Columns available in the GO datasets:
- id: Protein ID (Uniprot accession)
- exp: Annotations from experimental evidence
- exp_not: Negative experimental annotations ("NOT has FunctionA")
- phylo: Annotations from phylogenetic evidence
- phylo_not: Negative phylogenetic annotations
- curated: Annotations from curated evidence
- curated_not: Negative curated annotations
- comp: Annotations from computational analysis evidence
- comp_not: Negative computational analysis annotations
- iea: Annotations from inferred from electronic annotation evidence
- derived_not: Derived negative annotations, from the Warwick 2020 study (see references below)

The GO term columns have terms separated by "," as values. For example: "GO:0003674,GO:0008152,GO:0005575".

## DeepLoc 2.1 dataset

Columns in the DeepLoc dataset. These annotation columns are binary (True / False):
- id: Protein ID (Uniprot accession)
- Membrane
- Cytoplasm
- Nucleus
- Extracellular
- Cell membrane
- Mitochondrion
- Plastid
- Endoplasmic reticulum
- Lysosome/Vacuole
- Golgi apparatus
- Peroxisome
- Peripheral
- Transmembrane
- LipidAnchor
- Soluble

## Citation

New manuscript being written. For now, cite the previous release:

@inproceedings{alves2025protein,
  title        = {Protein Dimension DB: A Unified Protein Repository for Representation Learning and Functional Analysis},
  author       = {Alves Sobrinho, Pitágoras de Azevedo and Sakamoto, Tetsu and Figuerola, Wilfredo Blanco},
  booktitle    = {Bioinformatics and Computational Biology: 21st Brazilian Conference, X-Meeting 2025, João Pessoa, Brazil, June 3--6, 2025, Proceedings},
  year         = {2025},
  organization = {Springer}
}

## References

- Ødum, M.T.; Teufel, F.; Thumuluri, V.; Almagro Armenteros, J.J.; Johansen, A.R.; Winther, O.; Nielsen, H. DeepLoc 2.1: multi-label membrane protein type prediction using protein language models. Nucleic Acids Research 2024, 52, W215–W220;
- Huntley, R.P.; Sawford, T.; Mutowo-Meullenet, P.; Shypitsyna, A.; Bonilla, C.; Martin, M.J.; O’Donovan, C. The GOA database: gene ontology annotation updates for 2015. Nucleic acids research 2015, 43, D1057–D1063;
- Warwick Vesztrocy, A.; Dessimoz, C. Benchmarking gene ontology function predictions using negative annotations. Bioinformatics 2020, 36, i210–i218;

"""

datasets = {
    "biological_process-uniprot_swissprot-pddb2": {
        "file": "go.bp.parquet",
        "ont_name": "Biological Process"
    },
    "molecular_function-uniprot_swissprot-pddb2": {
        "file": "go.mf.parquet",
        "ont_name": "Molecular Function"
    },
    "cellular_component-uniprot_swissprot-pddb2": {
        "file": "go.cc.parquet",
        "ont_name": "Cellular Component"
    },
    "deeploc_2.1-uniprot_swissprot-pddb2": {
        "file": "deeploc.parquet",
        "ont_name": "DeepLoc 2.1"
    },
}

for rawline in open(".env", "r"):
    line = rawline.strip()
    if line.startswith("#") or line == "":
        continue
    parts = line.split("=", 1)
    if len(parts) == 2:
        os.environ[parts[0]] = parts[1]

USER_NAME = "pitagoras-alves"
HF_TOKEN = os.getenv("HF_TOKEN")

def upload_to_huggingface(hf_dataset, repo_name):
    print("Iniciando upload para o Hugging Face Hub...")

    # CONFIGURAÇÕES DO UPLOAD
    REPO_ID = USER_NAME + "/" + repo_name
    IS_PRIVATE = False
    url = f"https://huggingface.co/datasets/{REPO_ID}"
    print(f"\tURL: {url}")

    # O método push_to_hub do dataset cuida de criar o repo se não existir
    hf_dataset.push_to_hub(
        REPO_ID,
        token=HF_TOKEN,
        private=IS_PRIVATE
    )

    print(
        f"Upload concluído com sucesso! Veja em: https://huggingface.co/datasets/{REPO_ID}"
    )

if len(sys.argv) >= 2:
    pddb_release_dir = sys.argv[1]
else:
    pddb_release_dir = "."

for reponame, details in datasets.items():
    print("Starting to upload", reponame)
    readme = basecard.replace("[ont_name]", details["ont_name"])
    ont_tag = details["ont_name"].lower().replace(" ", "-").replace("-2.1", "")
    go_tags = ["biology", "uniprot", "swissprot", "protein", "go", "gene-ontology"]
    deeploc_tags = ["protein", "deeploc", "localization", "cellular-localization"]

    parquet_path = details['file']

    df = Dataset.from_parquet(pddb_release_dir+'/'+parquet_path)
    #upload_to_huggingface(df, reponame)

    if "deeploc" in ont_tag:
        tags = deeploc_tags
    else:
        tags = go_tags
    tags += [ont_tag]
    card_data = DatasetCardData(
        language="en",
        license="mit",
        task_categories=[
            "tabular-classification",
        ],
        tags=tags,
    )
    template_path = "dataset_card.md"
    with open(template_path, "w") as f:
        f.write(readme)
    card = DatasetCard.from_template(
        card_data=card_data, template_path=template_path
    )
    card.push_to_hub(USER_NAME + "/" + reponame, token=HF_TOKEN)

