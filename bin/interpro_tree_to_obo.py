#!/usr/bin/env python3

import sys

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


if __name__ == "__main__":
    rawtree_path = sys.argv[1]
    annotation_tsv_path = sys.argv[2]
    obo_path = sys.argv[3]

    entities_tree_dict = {}
    last_seen_at_level = {}

    previous_id = "Root"
    previous_level = -1

    # Add fictional root
    entities_tree_dict["Root"] = {
        "name": "InterPro Logical Root",
        "parent": None,
        "children": [],
    }

    for rawline in open(rawtree_path, "r"):
        parts = rawline.strip().split("::")
        if len(parts) < 2:
            continue
        node_id = parts[0]
        node_name = parts[1]
        actual_id = node_id.lstrip("-")

        current_level = (len(node_id) - len(actual_id)) // 2

        if current_level == 0:
            parent_id = "Root"
        else:
            parent_id = last_seen_at_level.get(current_level - 1)

        entities_tree_dict[actual_id] = {
            "name": node_name,
            "parent": parent_id,
            "children": [],
        }

        # 4. Link this new node to its parent's 'children' list
        if parent_id and parent_id in entities_tree_dict:
            entities_tree_dict[parent_id]["children"].append(actual_id)

        # 5. Update the tracker so future children know this node exists at this depth
        last_seen_at_level[current_level] = actual_id
    print(f"{len(entities_tree_dict)} nodes in tree")
    # Look for unmentioned terms
    with open(annotation_tsv_path, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            protein_id = parts[0]
            interpro_ids = parts[1].split(";")
            for interpro_id in interpro_ids:
                if interpro_id not in entities_tree_dict:
                    entities_tree_dict[interpro_id] = {
                        "name": protein_id,
                        "parent": "Root",
                        "children": [],
                    }
    print(f"{len(entities_tree_dict)} nodes in tree after adding unmentioned terms")

    """for node_id, data in entities_tree_dict.items():
        print(f"ID: {node_id}")
        print(f"  Name:   {data['name']}")
        print(f"  Parent: {data['parent']}")
        print(f"  Children: {data['children']}\n")"""

    # Create a new OBO graph
    obo_output_stream = open(obo_path, "w")
    obo_output_stream.write(
        "format-version: 1.2\ndata-version: releases/2025-03-16\nontology: interpro\n"
        + "idspace: dc http://purl.org/dc/elements/1.1/ \n"
        + "idspace: oboInOwl http://www.geneontology.org/formats/oboInOwl# \n"
        + "idspace: terms http://purl.org/dc/terms/ \n"
        + "property_value: owl:versionInfo 2025-03-16 xsd:string\n"
        + "property_value: has_ontology_root_term Root\n"
        + "property_value: terms:license http://creativecommons.org/licenses/by/4.0/\n\n"
    )

    for node_id, data in entities_tree_dict.items():
        obo_output_stream.write("[Term]\n")
        obo_output_stream.write(f"id: {node_id}\n")
        obo_output_stream.write(f"name: {data['name']}\n")
        obo_output_stream.write(f"namespace: interpro\n")
        if data["parent"] is not None:
            obo_output_stream.write(
                f"is_a: {data['parent']} ! {entities_tree_dict[data['parent']]['name']}\n"
            )
        obo_output_stream.write("\n")
    obo_output_stream.close()
