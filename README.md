# 🧬🖥 Protein Dimension DB 🖥🧬

Scientific data lake with PLM embeddings, GO annotations and taxonomy representations for all proteins in Uniprot/Swiss-Prot

## Current Release (2)

Proteins are sorted by length. All files contain the same sequence of proteins, to make joins and merge operations easier.

### Protein Language Model Embeddings 🔢

Several models are used to create computational descriptions (embeddings) of the Swiss-Prot proteins:

| Model 🤖 | Vector Length 📏 | Download Links (By Pooling Method) 🔗 |
| :--- | :---: | :--- |
| ElnaggarLab/ankh-base | 768 | [Mean (1.5G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.ankh_base_mean.parquet), [Parti (3.0G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.ankh_base_parti.parquet), [Max (1.4G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.ankh_base_max.parquet), [STD (1.4G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.ankh_base_std.parquet) |
| ElnaggarLab/ankh-large | 1536 | [Mean (2.9G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.ankh_large_mean.parquet), [Parti (6.0G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.ankh_large_parti.parquet), [Max (2.7G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.ankh_large_max.parquet), [STD (2.7G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.ankh_large_std.parquet) |
| ElnaggarLab/ankh2-ext2 | 1536 | [Mean (2.9G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.ankh2_large_mean.parquet), [Parti (6.0G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.ankh2_large_parti.parquet), [Max (2.7G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.ankh2_large_max.parquet), [STD (2.8G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.ankh2_large_std.parquet) |
| ElnaggarLab/ankh3-large | 1536 | [Mean (2.9G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.ankh3_large_mean.parquet), [Parti (6.0G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.ankh3_large_parti.parquet), [Max (2.7G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.ankh3_large_max.parquet), [STD (2.7G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.ankh3_large_std.parquet) |
| facebook/esm2_t30_150M_UR50D | 640 | [Mean (1.3G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.esm2_150_mean.parquet), [Parti (2.5G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.esm2_150_parti.parquet), [Max (1.2G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.esm2_150_max.parquet), [STD (1.2G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.esm2_150_std.parquet) |
| facebook/esm2_t33_650M_UR50D | 1280 | [Mean (2.4G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.esm2_650_mean.parquet), [Parti (5.0G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.esm2_650_parti.parquet), [Max (2.3G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.esm2_650_max.parquet), [STD (2.3G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.esm2_650_std.parquet) |
| biohub/ESMC-300M | 960 | [Mean (1.8G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.esmc_300_mean.parquet), [Parti (3.7G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.esmc_300_parti.parquet), [Max (1.7G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.esmc_300_max.parquet), [STD (1.7G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.esmc_300_std.parquet) |
| biohub/ESMC-600M | 1152 | [Mean (2.2G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.esmc_600_mean.parquet), [Parti (4.5G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.esmc_600_parti.parquet), [Max (2.1G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.esmc_600_max.parquet), [STD (2.1G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.esmc_600_std.parquet) |
| Profluent-Bio/E1-150m | 768 | [Mean (1.5G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.e1_150_mean.parquet), [Parti (3.0G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.e1_150_parti.parquet), [Max (587M)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.e1_150_max.parquet), [STD (1.4G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.e1_150_std.parquet) |
| Profluent-Bio/E1-300m | 1024 | [Mean (2.0G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.e1_300_mean.parquet), [Parti (4.0G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.e1_300_parti.parquet), [Max (757M)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.e1_300_max.parquet), [STD (1.8G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.e1_300_std.parquet) |
| Profluent-Bio/E1-600m | 1280 | [Mean (2.4G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.e1_600_mean.parquet), [Parti (5.0G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.e1_600_parti.parquet), [Max (929M)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.e1_600_max.parquet), [STD (2.3G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.e1_600_std.parquet) |
| flair-bio/amplify-120m | 640 | [Mean (1.3G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.amplify_120_mean.parquet), [Parti (2.5G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.amplify_120_parti.parquet), [Max (1.2G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.amplify_120_max.parquet), [STD (1.2G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.amplify_120_std.parquet) |
| flair-bio/amplify-350m | 960 | [Mean (1.9G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.amplify_350_mean.parquet), [Parti (3.7G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.amplify_350_parti.parquet), [Max (1.8G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.amplify_350_max.parquet), [STD (1.8G)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.amplify_350_std.parquet) |

### Autoencoder Embeddings and One-hot Encodings 🔢

Numerical representations of the NCBI taxon IDs and InterproScan categories of each protein. Instead of the original NCBI taxonomy tree, we use the custom taxonomy created by [taxallnomy](https://github.com/tetsufmbio/taxallnomy) project, because it attributes the same number of parent taxa (genus, family, order...) to each species ID. 

| Name | Description | Vector Length 📏 | Download Links 🔗 |
| :--- | :--- | :---: | :--- |
| Interpro Autoencoder | Encoding of the top 17000 most common InterproScan categories in SwissProt proteins | 32 | [Embeddings (36M)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.interpro_autoencoded.parquet) |
| TaxID Autoencoder | Encoding of the top 6006 most common Taxonomic IDs in SwissProt proteins | 32 | [Embeddings (5.1M)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.taxid_autoencoded.parquet) |
| emb.taxa_profile_256.parquet | Taxa One-Hot Encoding | 256 | [Encodings (62M)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.taxa_profile_256.parquet) |
| emb.taxa_profile_128.parquet | Taxa One-Hot Encoding | 128 | [Encodings (33M)](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_2/emb.taxa_profile_128.parquet) |

### Protein Annotations 📚

Gene Ontology annotations of Swiss-Prot proteins, separated by evidence type. We also have made available a parsed version of the DeepLoc dataset.

Gene ontology evidence code groups included:
- exp: EXP, IMP, IGI, IPI, IDA, IEP, HTP, HDA, HMP, HGI, HEP
- phylo: IBA, IBD, IKR, IRD
- curated: IC, TAS
- comp: ISS, ISO, ISA, ISM, IGP, RCA
- iea: IEA

Gene ontology annotation columns:
- id: Uniprot ID;
- exp, phylo, curated, comp, iea: GO ID list of positive annotations;
- exp_not, phylo_not, curated_not, comp_not, iea_not: Negative annotations (NOTs in GOA);
- derived_not: NOTs derived in the *Warwick and Dessimoz (2020)* article;

DeepLoc annotation columns (subcellular locations and membrane protein types) from the *Ødum et al. (2024)* article:
- id: Uniprot ID;
- Membrane, Cytoplasm, Nucleus, Extracellular, Cell membrane, Mitochondrion, Plastid, Endoplasmic reticulum, Lysosome/Vacuole, Golgi apparatus, Peroxisome, Peripheral, Transmembrane, LipidAnchor, Soluble: True / False values;


|            Name           |        Content       | Download Links 🔗 |
|:-:|:-:|:-:|
| go.expanded.tsv.gz        | Simplified version of GOA. Columns: Uniprot ID, GO ID, Evidence Code, Taxon ID and Ontology |  [25M]() |
| go.mf.parquet | Molecular Functions | [14M]() |
| go.bp.parquet | Biological Processes | [21M]() |
| go.cc.parquet | Cellular Components  | [9.9M]() |
| deeploc.parquet | Subcellular Locations (DeepLoc dataset) | [173KB]() |
| interpro.tsv | InterproScan categories of SwissProt proteins | [27M]() |
| taxid.tsv | Taxonomic ID of each protein in SwissProt. Columns: uniprot_id, taxid, lineage | [46M]() |

### Uniprot/Swiss-Prot 🔬

|           Name          |                  Content                  | Download Links 🔗 |
|:-----------------------:|:-----------------------------------------:|:-------------:|
|         ids.txt         |           Uniprot Accession IDs           |       [10K]()       |
| sequences.swissprot.fasta | Aminoacid sequences of SwissProt proteins |       [201M]()       |


### Others

|            Name           |        Content       | Download Links 🔗 |
|:-------------------------:|:--------------------:|:-------------:|
| taxallnomy.parquet | Parent TaxonIDs in Taxallnomy for each NCBI taxon ID |       [306M]()       |
| taxid.obo | NCBI taxonomy graph in OBO format |       [2.0M]()       |
| interpro.obo | Interpro categories graph in OBO format | [4.1M]() |


## Citation

Please cite the following work:

Bibtext:
```bibtex
@inproceedings{AlvesSobrinho2025ProteinDimensionDB,
  author       = {Pit{\'{a}}goras de Azevedo Alves Sobrinho and Tetsu Sakamoto and Wilfredo Blanco Figuerola},
  title        = {Protein Dimension DB: A Unified Protein Repository for Representation Learning and Functional Analysis},
  booktitle    = {BioInformatics: 21st Brazilian Congress, X-Meeting 2025, João Pessoa, Brazil, June 3–6, 2025, Proceedings},
  series       = {Lecture Notes in Computer Science},
  volume       = {16037},
  year         = {2025},
  editor       = {Marcio Dorn and Fabricio Martins Lopes},
  publisher    = {Springer Cham},
  isbn         = {978-3-032-09335-6},
  eisbn        = {978-3-032-09336-3},
  address      = {Cham, Switzerland}
}
```

APA reference:
> Alves Sobrinho, P. de A., Sakamoto, T., & Blanco Figuerola, W. (2025). Protein Dimension DB: A unified protein repository for representation learning and functional analysis. BioInformatics: 21st Brazilian Congress, X-Meeting 2025, João Pessoa, Brazil, June 3–6, 2025, Proceedings (Lecture Notes in Computer Science, Vol. 16037). Springer Cham.

## References

[1] Alex Warwick Vesztrocy and Christophe Dessimoz. "Benchmarking gene ontology function predictions using
negative annotations", Bioinformatics, 36, 2020, i210–i218,
[doi: 10.1093/bioinformatics/btaa466](https://academic.oup.com/bioinformatics/article/36/Supplement_1/i210/5870480);

[2] Marius Thrane Ødum, Felix Teufel, Vineet Thumuluri, et al. "DeepLoc 2.1: multi-label membrane protein type prediction using protein language models", Nucleic Acids Research, Volume 52, Issue W1, 5 July 2024, Pages W215–W220, [doi: 10.1093/nar/gkae237](https://academic.oup.com/nar/article/52/W1/W215/7642068);

[3] Tetsu Sakamoto and Miguel Ortega. "Taxallnomy Database", Laboratório de Biodados, UFMG. [URL](http://bioinfo.icb.ufmg.br/taxallnomy/);