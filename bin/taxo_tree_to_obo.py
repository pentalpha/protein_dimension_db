#!/usr/bin/env python3

import sys
import pandas as pd

"""
Tree raw text format:
uniprot_id	taxid	lineage
P0DPR3	9606	2759;33208;6072;33511;7711;89593;7742;8287;40674;32525;9347;1437010;314146;9443;376913;314293;9526;314295;9604;207598;9605;9606
P83570	6610	2759;33208;6072;2697495;6447;6605;6606;215450;551287;551351;6608;6609;6610
P62968	9823	2759;33208;6072;33511;7711;89593;7742;8287;40674;32525;9347;1437010;314145;91561;35497;9821;9822;9823
P62969	9940	2759;33208;6072;33511;7711;89593;7742;8287;40674;32525;9347;1437010;314145;91561;9845;35500;9895;9963;9935;9940
P62970	8346	2759;33208;6072;33511;7711;89593;7742;8287;8292;41666;8342;30312;8344;8346

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
    obo_path = sys.argv[2]

    entities_tree_dict = {}

    # Add fictional root
    entities_tree_dict["1"] = {
        "name": "NCBI Taxa Logical Root (LUCA)",
        "parents": set(),
    }

    rawtree_df = pd.read_csv(rawtree_path, sep="\t")

    for _, row in rawtree_df.iterrows():
        node_id = row["uniprot_id"]
        parts = row["lineage"].split(";")

        for current_level in range(len(parts)):
            if current_level > 0:
                parent_ids = [parts[current_level - 1]]
            else:
                parent_ids = []

            node_name = parts[current_level]

            if node_name not in entities_tree_dict:
                entities_tree_dict[node_name] = {
                    "name": node_name,
                    "parents": set(parent_ids),
                }
            else:
                entities_tree_dict[node_name]["parents"].update(parent_ids)

    print(f"{len(entities_tree_dict)} nodes in tree")

    for node_id, data in entities_tree_dict.items():
        if len(data["parents"]) == 0 and node_id != "1":
            data["parents"] = ["1"]
        """print(f"  Name:   {data['name']}")
        print(f"  Parent: {data['parent']}\n")"""

    # Create a new OBO graph
    obo_output_stream = open(obo_path, "w")
    obo_output_stream.write(
        "format-version: 1.2\ndata-version: releases/2025-03-16\nontology: ncbi_taxid\n"
        + "idspace: dc http://purl.org/dc/elements/1.1/ \n"
        + "idspace: oboInOwl http://www.geneontology.org/formats/oboInOwl# \n"
        + "idspace: terms http://purl.org/dc/terms/ \n"
        + "property_value: owl:versionInfo 2025-03-16 xsd:string\n"
        + "property_value: has_ontology_root_term LUCA\n"
        + "property_value: terms:license http://creativecommons.org/licenses/by/4.0/\n\n"
    )

    for node_id, data in entities_tree_dict.items():
        obo_output_stream.write("[Term]\n")
        obo_output_stream.write(f"id: {node_id}\n")
        obo_output_stream.write(f"name: {data['name']}\n")
        obo_output_stream.write(f"namespace: ncbi_taxid\n")
        for parent_id in data["parents"]:
            obo_output_stream.write(
                f"is_a: {parent_id} ! {entities_tree_dict[parent_id]['name']}\n"
            )
        obo_output_stream.write("\n")
    obo_output_stream.close()
