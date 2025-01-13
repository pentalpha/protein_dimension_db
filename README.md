# 🧬🖥 Protein Dimension DB 🖥🧬

Datasets with PLM embeddings, GO annotations and taxonomy representations for all proteins in Uniprot/Swiss-Prot

## Current Release

Proteins are sorted by length. All files contain the same sequence of proteins, so the "ids.txt" file can be used as the row names.

### Uniprot/Swiss-Prot 🔬

|           Name          |                  Content                  | Download Links 🔗 |
|:-----------------------:|:-----------------------------------------:|:-------------:|
|         ids.txt         |           Uniprot Accession IDs           |       [UFRN](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_1/ids.txt)       |
| uniprot_sorted.fasta.gz | Aminoacid sequences of SwissProt proteins |       [UFRN](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_1/uniprot_sorted.fasta.gz)       |
| taxid.tsv               | NCBI taxon ID of each protein             |       [UFRN](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_1/taxid.tsv)       |

### Protein Annotations 📚

All Gene Ontology annotations of Swiss-Prot proteins, excluding computational, non-traceable and no-data annotations. The full list of ignored evidence codes is available at [evi_not_to_use.txt](evi_not_to_use.txt). Annotations have been "expanded upwards": parent terms of existing annotations have been included in these files.

|            Name           |        Content       | Download Links 🔗 |
|:-------------------------:|:--------------------:|:-------------:|
| go.expanded.tsv.gz        | MF, BP and CC annotations in simplified GAF format |  [UFRN](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_1/go.experimental.cc.tsv.gz)       |
| go.experimental.mf.tsv.gz |  Molecular Functions |       [UFRN](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_1/go.experimental.cc.tsv.gz)       |
| go.experimental.bp.tsv.gz | Biological Processes |       [UFRN](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_1/go.experimental.bp.tsv.gz)       |
| go.experimental.cc.tsv.gz | Cellular Components  |       [UFRN](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_1/go.experimental.mf.tsv.gz)       |

### Protein Language Model Embeddings 🔢

Several models are used to create computational descriptions of the Swiss-Prot proteins:

|             Name            |                                          Model 🤖                                         | Vector Length 📏 | Download Links 🔗 |
|:---------------------------:|:--------------------------------------------------------------------------------------:|:-------------:|:-------------:|
|     emb.prottrans.npy.gz    | prottrans_t5_xl_u50 (calculated by [Uniprot](https://www.uniprot.org/help/embeddings)) |      1024     |   [UFRN](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_1/emb.prottrans.npy.gz)  |
|     emb.ankh_large.npy    | ankh-large |      1536     |   [UFRN](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_1/emb.ankh_large.npy.gz)  |
|     emb.ankh_base.npy    | ankh-base |      768     |   [UFRN](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_1/emb.ankh_base.npy.gz)  |
|     emb.esm2_t36.npy.gz     |                                   esm2_t36_3B_UR50D                                  |      2560     |    Upcoming   |
|     emb.esm2_t33.npy.gz     |                                   esm2_t33_650M_UR50D                                  |      1280     |    [UFRN](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_1/emb.esm2_33.npy.gz)   |
|     emb.esm2_t30.npy.gz     |                                   esm2_t30_150M_UR50D                                  |      640      |    [UFRN](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_1/emb.esm2_30.npy.gz)   |
|     emb.esm2_t12.npy.gz     |                                   esm2_t12_35M_UR50D                                   |      480      |    [UFRN](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_1/emb.esm2_t12.npy.gz)   |
|      emb.esm2_t6.npy.gz     |                                    esm2_t6_8M_UR50D                                    |      320      |    [UFRN](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_1/emb.esm2_t6.npy.gz)   |

### Taxonomy 🔢

Numerical representations of the NCBI taxon IDs of each protein. Instead of the original NCBI taxonomy tree, we use the custom taxonomy created by [taxallnomy](https://github.com/tetsufmbio/taxallnomy) project, because it attributes the same number of parent taxa (genus, family, order...) to each species ID. 

|             Name            |                                          Description                                         | Vector Length 📏 | Download Links 🔗 |
|:---------------------------:|:--------------------------------------------------------------------------------------:|:-------------:|:-------------:|
| emb.taxa_profile_256.npy.gz |                                     Taxa Proximity [0.0, 1.0] to each one of the 256 most annotated taxa                                    |      256      |    [UFRN](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_1/emb.taxa_profile_256.npy.gz)   |
| emb.taxa_profile_128.npy.gz |                                     Taxa Proximity [0.0, 1.0] to each one of the 128 most annotated taxa                                    |      128      |    [UFRN](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_1/emb.taxa_profile_128.npy.gz)   |
|    onehot.taxa_256.npy.gz   |                                  Taxa One-Hot Encoding                                 |      256      |    [UFRN](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_1/onehot.taxa_256.npy.gz)   |
|    onehot.taxa_128.npy.gz   |                                  Taxa One-Hot Encoding                                 |      128      |    [UFRN](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_1/onehot.taxa_128.npy.gz)   |

## File Formats 🗃️

| Files                    | Format Descriptions                                                                                                |
|--------------------------|--------------------------------------------------------------------------------------------------------------------|
| ids.txt                  | One UniprotID per line                                                                                             |
| taxid.tsv                | Tab-separated table with columns: UniprotID, NCBI Taxon ID                                                         |
| go.expanded.tsv.gz       | Tab-separated table with columns: UniprotID, GO ID, ECO ID, NCBI Taxon ID, GO Ontology Code                        |
| go.experimental.*.tsv.gz | Tab-separated table with columns: UniprotID, GO IDs separated by ','                                               |
| emb.*.npy.gz             | Numpy matrix compressed with gzip. For rows where an embedding could not be defined, a vector of np.NaN is placed. |

## Create Release

Requirements to generate the datasets from scratch:
- Nextflow >= 24
- Mamba package manager
- Fast and stable internet connection to download original datasets
- At least 16GB of RAM

```
$ nextflow run prepare_requirements.nf
$ nextflow run new_release.nf --release_dir <path to generate database at>
```