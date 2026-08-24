import requests

tree_text_url = "https://ftp.ebi.ac.uk/pub/databases/interpro/releases/latest/ParentChildTreeFile.txt"
"""
Tree raw text format:
IPR000040::Acute myeloid leukemia 1 protein (AML1)/Runt::
--IPR016554::Runt-related transcription factor RUNX::
IPR000053::Thymidine/pyrimidine-nucleoside phosphorylase::
--IPR013466::Thymidine phosphorylase/AMP phosphorylase::
----IPR017713::AMP phosphorylase::
----IPR028579::Putative thymidine phosphorylase::
--IPR018090::Pyrimidine-nucleoside phosphorylase, bacterial/eukaryotic::
----IPR013465::Thymidine phosphorylase::
IPR000056::Ribulose-phosphate 3-epimerase-like::
--IPR026019::Ribulose-phosphate 3-epimerase::
--IPR043677::D-allulose-6-phosphate 3-epimerase::

obo format example:
[Term]
id: GO:0000008
name: obsolete thioredoxin
namespace: molecular_function
alt_id: GO:0000013
def: "OBSOLETE. A small disulfide-containing redox protein that serves as a general protein disulfide oxidoreductase. Interacts with a broad range of proteins by a redox mechanism, based on the reversible oxidation of 2 cysteine thiol groups to a disulfide, accompanied by the transfer of 2 electrons and 2 protons. The net result is the covalent interconversion of a disulfide and a dithiol." [GOC:kd]
comment: This term was made obsolete because it represents gene products.
synonym: "thioredoxin" EXACT []
is_obsolete: true
consider: GO:0003756
consider: GO:0015036

[Term]
id: GO:0000009
name: alpha-1,6-mannosyltransferase activity
namespace: molecular_function
def: "Catalysis of the transfer of a mannose residue to an oligosaccharide, forming an alpha-(1->6) linkage." [GOC:mcc, PMID:2644248]
synonym: "1,6-alpha-mannosyltransferase activity" EXACT []
xref: Reactome:R-HSA-449718 "Addition of a third mannose to the N-glycan precursor by ALG2"
is_a: GO:0000030 ! mannosyltransferase activity

[Term]
id: GO:0000010
name: heptaprenyl diphosphate synthase activity
namespace: molecular_function
alt_id: GO:0036422
def: "Catalysis of the reaction: (2E,6E)-farnesyl diphosphate + 4 isopentenyl diphosphate = 4 diphosphate + all-trans-heptaprenyl diphosphate." [PMID:9708911, RHEA:27794]
synonym: "all-trans-heptaprenyl-diphosphate synthase activity" EXACT [EC:2.5.1.30]
synonym: "HepPP synthase activity" RELATED [EC:2.5.1.30]
synonym: "heptaprenyl pyrophosphate synthase activity" RELATED [EC:2.5.1.30]
synonym: "heptaprenyl pyrophosphate synthetase activity" RELATED [EC:2.5.1.30]
synonym: "trans-hexaprenyltranstransferase activity" EXACT []
xref: EC:2.5.1.30
xref: MetaCyc:TRANS-HEXAPRENYLTRANSTRANSFERASE-RXN
xref: RHEA:27794
is_a: GO:0120531 ! prenyl diphosphate synthase activity
property_value: skos:exactMatch EC:2.5.1.30
property_value: skos:exactMatch RHEA:27794
"""
