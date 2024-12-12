# Protein Dimension DB

Datasets with embeddings and other representations for all proteins in Uniprot/Swiss-Prot

## Current Release

Proteins are sorted by length. All files contain the same sequence of proteins, so the "ids.txt" file can be used as the row names.

### Uniprot/Swiss-Prot

|           Name          |                  Content                  | Download Links |
|:-----------------------:|:-----------------------------------------:|:-------------:|
|         ids.txt         |           Uniprot Accession IDs           |       [UFRN](http://bioinformatics-brazil.org/~pitagoras/protein_dimension_db/release_1/ids.txt), [Zenodo](#)       |
| uniprot_sorted.fasta.gz | Aminoacid sequences of SwissProt proteins |       [UFRN](http://bioinformatics-brazil.org/~pitagoras/protein_dimension_db/release_1/uniprot_sorted.fasta.gz), [Zenodo](#)       |
| taxid.tsv               | NCBI taxon ID of each protein             |       [UFRN](http://bioinformatics-brazil.org/~pitagoras/protein_dimension_db/release_1/taxid.tsv), [Zenodo](#)       |

### Gene Ontology

All Gene Ontology annotations of Swiss-Prot proteins, excluding computational, non-traceable and no-data annotations. The full list of ignored evidence codes is available at [evi_not_to_use.txt](evi_not_to_use.txt). Annotations have been "expanded upwards": parent terms of existing annotations have been included in these files.

|            Name           |        Content       | Download Links |
|:-------------------------:|:--------------------:|:-------------:|
| go.expanded.tsv.gz        | MF, BP and CC annotations in simplified GAF format |  [UFRN](http://bioinformatics-brazil.org/~pitagoras/protein_dimension_db/release_1/go.experimental.cc.tsv.gz), [Zenodo](go.expanded.tsv.gz)       |
| go.experimental.mf.tsv.gz |  Molecular Functions |       [UFRN](http://bioinformatics-brazil.org/~pitagoras/protein_dimension_db/release_1/go.experimental.cc.tsv.gz), [Zenodo](#)       |
| go.experimental.bp.tsv.gz | Biological Processes |       [UFRN](http://bioinformatics-brazil.org/~pitagoras/protein_dimension_db/release_1/go.experimental.bp.tsv.gz), [Zenodo](#)       |
| go.experimental.cc.tsv.gz | Cellular Components  |       [UFRN](http://bioinformatics-brazil.org/~pitagoras/protein_dimension_db/release_1/go.experimental.mf.tsv.gz), [Zenodo](#)       |

### Embeddings

Several models are used to create computational descriptions of the Swiss-Prot proteins:

|             Name            |                                          Model                                         | Vector Length | Download Links |
|:---------------------------:|:--------------------------------------------------------------------------------------:|:-------------:|:-------------:|
|     emb.prottrans.npy.gz    | prottrans_t5_xl_u50 (calculated by [Uniprot](https://www.uniprot.org/help/embeddings)) |      1024     |   [UFRN](http://bioinformatics-brazil.org/~pitagoras/protein_dimension_db/release_1/emb.prottrans.npy.gz), [Zenodo](#)  |
|     emb.esm2_t33.npy.gz     |                                   esm2_t33_650M_UR50D                                  |      1280     |    Upcoming   |
|     emb.esm2_t30.npy.gz     |                                   esm2_t30_150M_UR50D                                  |      640      |    Upcoming   |
|     emb.esm2_t12.npy.gz     |                                   esm2_t12_35M_UR50D                                   |      480      |    Upcoming   |
|      emb.esm2_t6.npy.gz     |                                    esm2_t6_8M_UR50D                                    |      320      |    Upcoming   |
| emb.taxa_profile_256.npy.gz |                                     Taxa Profile v1                                    |      256      |    Upcoming   |
| emb.taxa_profile_128.npy.gz |                                     Taxa Profile v1                                    |      128      |    Upcoming   |
|    onehot.taxa_128.npy.gz   |                                  Taxa One-Hot Encoding                                 |      128      |    Upcoming   |

## File Formats

| Files                    | Format Descriptions                                                                                                |
|--------------------------|--------------------------------------------------------------------------------------------------------------------|
| ids.txt                  | One UniprotID per line                                                                                             |
| taxid.tsv                | Tab-separated table with columns: UniprotID, NCBI Taxon ID                                                         |
| go.expanded.tsv.gz       | Tab-separated table with columns: UniprotID, GO ID, ECO ID, NCBI Taxon ID, GO Ontology Code                        |
| go.experimental.*.tsv.gz | Tab-separated table with columns: UniprotID, GO IDs separated by ','                                               |
| emb.*.npy.gz             | Numpy matrix compressed with gzip. For rows where an embedding could not be defined, a vector of np.NaN is placed. |

## Create Release

Requirements:
- Nextflow >= 24
- Mamba package manager
- Internet connection to download original datasets

```
$ nextflow run prepare_requirements.nf
$ nextflow run new_release.nf --release_dir <path to generate database at>
```