# 🧬🖥 Protein Dimension DB 🖥🧬

Scientific data lake with PLM embeddings, GO annotations and taxonomy representations for all proteins in Uniprot/Swiss-Prot

## Release 2

Proteins are sorted by length. All files contain the same sequence of proteins, to make joins and merge operations easier.

### Protein Language Model Embeddings 🔢

Several models are used to create computational descriptions (embeddings) of the Swiss-Prot proteins:

|                                          Model 🤖                                         | Vector Length 📏 | File Size | Download Links (By Pooling Method) 🔗 |
:--------------------------------------------------------------------------------------:|:------------:|:------------:|:------------:|
ElnaggarLab/ankh-base |      768     | 1.4G to 3G | [Mean](), [Parti](), [Max](), [STD]()  |
ElnaggarLab/ankh-large |      1536     | 2.7G to 6G | [Mean](), [Parti](), [Max](), [STD]()  |
ElnaggarLab/ankh2-ext2 |      1536     | 3.4G | [Mean]( ), [Parti]( ), [Max](), [STD]()  |
ElnaggarLab/ankh3-large |      1536     | 3.4G | [Mean]( ), [Parti]( ), [Max](), [STD]()  |
facebook/esm2_t30_150M_UR50D | X | XG | [Mean]( ), [Parti]( ), [Max](), [STD]()  |
| facebook/esm2_t33_650M_UR50D |      1280     | 2.8G | [Mean]( ), [Parti]( ), [Max](), [STD]()  |
|     facebook/esm2_t36_3B_UR50D | X | XG | [Mean]( ), [Parti]( ), [Max](), [STD]()  |
biohub/ESMC-300M | X | XG | [Mean]( ), [Parti]( ), [Max](), [STD]()  |
biohub/ESMC-600M | X | XG | [Mean]( ), [Parti]( ), [Max](), [STD]()  |
biohub/esm3-sm-open-v1 | X | XG | [Mean]( ), [Parti]( ), [Max](), [STD]()  |
Profluent-Bio/E1-300m | 768 | 587M to 3G | [Mean]( ), [Parti]( ), [Max](), [STD]()  |
Profluent-Bio/E1-600m | 768 | 587M to 3G | [Mean]( ), [Parti]( ), [Max](), [STD]()  |
flair-bio/amplify-350m | X | XG | [Mean]( ), [Parti]( ), [Max](), [STD]()  |
hugohrban/progen2-medium | X | XG | [Mean]( ), [Parti]( ), [Max](), [STD]()  |
oriel9p/protsent-esm2-150M | X | XG | [Mean]( ), [Parti]( ), [Max](), [STD]()  |
AI4PD/ProtGPT3-1.3B | X | XG | [Mean]( ), [Parti]( ), [Max](), [STD]()  |

### Autoencoder Embeddings and One-hot Encodings 🔢

Numerical representations of the NCBI taxon IDs and InterproScan categories of each protein. Instead of the original NCBI taxonomy tree, we use the custom taxonomy created by [taxallnomy](https://github.com/tetsufmbio/taxallnomy) project, because it attributes the same number of parent taxa (genus, family, order...) to each species ID. 

|             Name            |                                          Description                                         | Vector Length 📏 | Download Links 🔗 |
|:---------------------------:|:--------------------------------------------------------------------------------------:|:-------------:|:-------------:|
| Interpro Autoencoder | Encoding of the top 17000 most common InterproScan categories in SwissProt proteins |      32      |    [Embeddings (36M)](), [Model (637M)]()   |
| TaxID Autoencoder | Encoding of the top 6006 most common Taxonomic IDs in SwissProt proteins |      32      |    [Embeddings (5M)](), [Model (181M)]()   |
|    onehot.taxa_256.parquet   |                                  Taxa One-Hot Encoding                                 |      256      |    [Encodings (3.8MB)](), [Taxon IDs (2.2K)]()   |
|    onehot.taxa_128.parquet   |                                  Taxa One-Hot Encoding                                 |      128      |    [Encodings (3.2MB)](), [Taxon IDs (1.1K)]()   |


### Uniprot/Swiss-Prot 🔬

|           Name          |                  Content                  | Download Links 🔗 |
|:-----------------------:|:-----------------------------------------:|:-------------:|
|         ids.txt         |           Uniprot Accession IDs           |       [UFRN](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_1/ids.txt)       |
| sequences.swissprot.fasta | Aminoacid sequences of SwissProt proteins |       [201M]()       |
| taxid.tsv               | NCBI taxon ID of each protein             |       [UFRN](https://ucrania.imd.ufrn.br/~pitagoras/protein_dimension_db/release_1/taxid.tsv)       |

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
|:-------------------------:|:--------------------:|:-------------:|
| go.expanded.tsv.gz        | Simplified version of GOA. Columns: Uniprot ID, GO ID, Evidence Code, Taxon ID and Ontology |  [25M]()       |
| go.mf.parquet | Molecular Functions |       [14M]()       |
| go.bp.parquet | Biological Processes |       [21M]()       |
| go.cc.parquet | Cellular Components  |       [9.9M]()       |
| deeploc.parquet | Subcellular Locations (DeepLoc dataset) | [173KB]()
| interpro.tsv | InterproScan categories of SwissProt proteins | [27M]()
### Others

|            Name           |        Content       | Download Links 🔗 |
|:-------------------------:|:--------------------:|:-------------:|
| taxallnomy.parquet | Parent TaxonIDs in Taxallnomy for each NCBI taxon ID |       [306M]()       |
| taxid.obo | NCBI taxonomy graph in OBO format |       [2.0M]()       |
| interpro.obo | Interpro categories graph in OBO format | [4.1M]() |

## References

[1] Alex Warwick Vesztrocy and Christophe Dessimoz. "Benchmarking gene ontology function predictions using
negative annotations", Bioinformatics, 36, 2020, i210–i218,
[doi: 10.1093/bioinformatics/btaa466](https://academic.oup.com/bioinformatics/article/36/Supplement_1/i210/5870480);

[2] Marius Thrane Ødum, Felix Teufel, Vineet Thumuluri, et al. "DeepLoc 2.1: multi-label membrane protein type prediction using protein language models", Nucleic Acids Research, Volume 52, Issue W1, 5 July 2024, Pages W215–W220, [doi: 10.1093/nar/gkae237](https://academic.oup.com/nar/article/52/W1/W215/7642068);

[3] Tetsu Sakamoto and Miguel Ortega. "Taxallnomy Database", Laboratório de Biodados, UFMG. [URL](http://bioinfo.icb.ufmg.br/taxallnomy/);