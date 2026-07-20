#!/usr/bin/env python
import gzip
from os import path
from collections import defaultdict
import sys
import pandas as pd
from tqdm import tqdm
import polars as pl
import json

from bioinfo_utils.gene_ontology import (
    gos_not_to_use,
    load_go_graph_strict,
    expand_go_list,
    list_ontology_members,
)
from bioinfo_utils.util_base import count_lines_large, open_file, write_file

# manual urls:
# quickgo:
# https://www.ebi.ac.uk/QuickGO/annotations?aspect=molecular_function&evidenceCode=ECO:0000352,ECO:0000269,ECO:0000314,ECO:0000315,ECO:0000316,ECO:0000353,ECO:0000270,ECO:0007005,ECO:0007001,ECO:0007003,ECO:0007007,ECO:0006056,ECO:0000318,ECO:0000320,ECO:0000321,ECO:0000304,ECO:0000305&evidenceCodeUsage=descendants&withFrom=UniProtKB
# geneontology.org
# http://current.geneontology.org/annotations/filtered_goa_uniprot_all_noiea.gaf.gz

"""
go name,code
Inferred from Experiment,EXP
Inferred from Direct Assay,IDA
Inferred from Physical Interaction,IPI
Inferred from Mutant Phenotype,IMP
Inferred from Genetic Interaction,IGI
Inferred from Expression Pattern,IEP
Inferred from High Throughput Experiment,HTP
Inferred from High Throughput Direct Assay,HDA
Inferred from High Throughput Mutant Phenotype,HMP
Inferred from High Throughput Genetic Interaction,HGI
Inferred from High Throughput Expression Pattern,HEP
Inferred from Biological aspect of Ancestor,IBA
Inferred from Biological aspect of Descendant,IBD
Inferred from Key Residues,IKR
Inferred from Rapid Divergence,IRD
Traceable Author Statement,TAS
Inferred by Curator,IC
"""

EVI_GROUPS = {
    "exp": [
        "EXP",
        "IMP",
        "IGI",
        "IPI",
        "IDA",
        "IEP",
        "HTP",
        "HDA",
        "HMP",
        "HGI",
        "HEP",
    ],
    "phylo": ["IBA", "IBD", "IKR", "IRD"],
    "curated": ["IC", "TAS"],
    "comp": ["ISS", "ISO", "ISA", "ISM", "IGP", "RCA"],
    "iea": ["IEA"],
}

EVI_TO_GROUP = {}
for group_name, evis in EVI_GROUPS.items():
    for evi in evis:
        EVI_TO_GROUP[evi] = group_name


def gafparsed_to_id2go(go_expanded_ann: pl.DataFrame, swissprot_ids, output_file):
    indexed_ann = {}
    onto_to_prefix = {"F": "mf", "P": "bp", "C": "cc"}

    print("Indexing by protein id")
    for row in tqdm(go_expanded_ann.rows(named=True), total=go_expanded_ann.height):
        protid = row["prot_id"]
        goid = row["goid"]
        evi = row["evi"]
        taxid = row["taxonid"]
        onto = row["aspect"]
        rel_str = row["relation"]
        if protid not in indexed_ann:
            indexed_ann[protid] = {
                "id": protid,
            }
            for group_name in EVI_GROUPS.keys():
                indexed_ann[protid]["mf_" + group_name] = set()
                indexed_ann[protid]["mf_" + group_name + "_not"] = set()
                indexed_ann[protid]["bp_" + group_name] = set()
                indexed_ann[protid]["bp_" + group_name + "_not"] = set()
                indexed_ann[protid]["cc_" + group_name] = set()
                indexed_ann[protid]["cc_" + group_name + "_not"] = set()

        prot_ann = indexed_ann[protid]
        if onto in ["F", "P", "C"]:
            prefix = onto_to_prefix[onto] + "_"
            if evi in EVI_TO_GROUP:
                group_name = EVI_TO_GROUP[evi]
                prefix = prefix + group_name
                if "NOT" in rel_str:
                    prefix = prefix + "_not"
                prot_ann[prefix].add(goid)

    lines = [v for k, v in indexed_ann.items()]

    print("Post-processing actions")
    for l in lines:
        if l["id"] in swissprot_ids:
            l["ProteinSet"] = "SwissProt"
        else:
            l["ProteinSet"] = "TrEMBL"
        for col in l.keys():
            if col != "id" and col != "ProteinSet":
                l[col] = list(l[col])
    print("Creating final dataframe")
    df = pl.DataFrame(lines)
    df.write_parquet(output_file)
    return df


if __name__ == "__main__":
    go_not_use_path = sys.argv[1]
    go_basic_path = sys.argv[2]
    goa_parsed = sys.argv[3]
    ids_path = sys.argv[4]
    derivated_negatives_path = sys.argv[5]
    output_dir = "./"

    goa_by_uniprot_path = output_dir + "/go.by_uniprot.parquet"

    print("Reading swissprot ids")
    uniprots_list = open(ids_path).read().split("\n")

    uniprots_set = set(uniprots_list)

    goa_parsed = pl.read_parquet(goa_parsed)

    goa_by_uniprot = gafparsed_to_id2go(goa_parsed, uniprots_set, goa_by_uniprot_path)

    goes_to_not_use = gos_not_to_use(go_not_use_path)
    go_graph = load_go_graph_strict(go_basic_path)
    gos_by_ontology, go_alt_ids = list_ontology_members(go_basic_path)

    print(f"Loading derived negatives")

    # an	GO_ID	infer_from_an	type	fam	uniprot
    derivated_negatives_df = pd.read_csv(derivated_negatives_path, sep="\t")
    derivated_negatives_df = derivated_negatives_df[
        derivated_negatives_df["type"] == -1
    ]

    # fill up GO_ID (only digits) so that 140096 -> GO:0140096 and 977 -> GO:0000977
    def add_zeros(goid, min_len: int = 7):
        return "GO:" + str(goid).zfill(min_len)

    derivated_negatives_df["GO_ID"] = derivated_negatives_df["GO_ID"].apply(add_zeros)

    print(f"Filtering by uniprot list")
    before = len(derivated_negatives_df)
    derivated_negatives_df = derivated_negatives_df[
        derivated_negatives_df["uniprot"].isin(uniprots_set)
    ]
    after = len(derivated_negatives_df)
    print(f"Removed {before - after} rows from {before}")
    # derivated_negatives_df = derivated_negatives_df[
    #    derivated_negatives_df["GO_ID"].isin(all_ont_goids)
    # ]
    print(derivated_negatives_df.head())
    derivated_neg_annots = defaultdict(set)
    derivated_negatives_df = derivated_negatives_df[["uniprot", "GO_ID"]]
    # Convert to tuples:
    uniprot_and_goid = derivated_negatives_df.itertuples(index=False, name=None)

    new_parsed = set()

    new_rows = {"mf": [], "bp": [], "cc": []}

    all_evi_groups = []
    for group_name in EVI_GROUPS.keys():
        all_evi_groups.append(group_name)
        all_evi_groups.append(group_name + "_not")
    # all_evi_groups.append("derived_not")

    annot_counts = {
        "mf": {n: 0 for n in all_evi_groups},
        "bp": {n: 0 for n in all_evi_groups},
        "cc": {n: 0 for n in all_evi_groups},
    }
    annot_counts["mf"]["derived_not"] = 0
    annot_counts["bp"]["derived_not"] = 0
    annot_counts["cc"]["derived_not"] = 0

    redundancy_counts = {
        "mf": {n: 0 for n in all_evi_groups},
        "bp": {n: 0 for n in all_evi_groups},
        "cc": {n: 0 for n in all_evi_groups},
    }
    redundancy_counts["mf"]["derived_not"] = 0
    redundancy_counts["bp"]["derived_not"] = 0
    redundancy_counts["cc"]["derived_not"] = 0

    derivated_negatives = {
        "mf": defaultdict(set),
        "bp": defaultdict(set),
        "cc": defaultdict(set),
    }

    for uniprot, goid in tqdm(uniprot_and_goid, desc="Listing derivated negatives"):
        if goid in go_alt_ids:
            goid = go_alt_ids[goid]
        ont = None
        if goid in gos_by_ontology["MF"]:
            ont = "mf"
        elif goid in gos_by_ontology["BP"]:
            ont = "bp"
        elif goid in gos_by_ontology["CC"]:
            ont = "cc"
        else:
            raise ValueError(f"GO: {goid} not found in any ontology")
        derivated_negatives[ont][uniprot].add(goid)

    print("Expanding GOs")
    uniprot_to_index = {"mf": {}, "bp": {}, "cc": {}}
    for row in tqdm(goa_by_uniprot.rows(named=True), total=goa_by_uniprot.height):
        for ont in ["mf", "bp", "cc"]:
            new_row = {"id": row["id"], "ProteinSet": row["ProteinSet"]}
            for evi_group in all_evi_groups:
                if "_not" in evi_group:
                    new_row[evi_group] = set(
                        expand_go_list(
                            row[ont + "_" + evi_group],
                            go_graph,
                            goes_to_not_use,
                            is_negative=True,
                        )
                    )
                else:
                    new_row[evi_group] = set(
                        expand_go_list(
                            row[ont + "_" + evi_group], go_graph, goes_to_not_use
                        )
                    )

            if row["id"] in derivated_negatives[ont]:
                new_row["derived_not"] = derivated_negatives[ont][row["id"]]
                del derivated_negatives[ont][row["id"]]
            else:
                new_row["derived_not"] = set()

            """to_remove = [
                {
                    "go_col": "negative",
                    "from": ["exp", "auto", "phylo", "derived_negative"],
                },
                {"go_col": "exp", "from": ["auto", "phylo", "derived_negative"]},
                {"go_col": "phylo", "from": ["auto", "derived_negative"]},
                {"go_col": "derived_negative", "from": ["auto"]},
            ]

            for removal in to_remove:
                source_set_name = removal["go_col"]
                source_set = new_row[source_set_name]
                for other_set_name in removal["from"]:
                    inter = source_set.intersection(new_row[other_set_name])
                    if len(inter) > 0:
                        new_row[other_set_name] -= inter
                        if (
                            not other_set_name
                            in redundancy_counts[ont][source_set_name]
                        ):
                            redundancy_counts[ont][source_set_name][other_set_name] = 0
                        redundancy_counts[ont][source_set_name][other_set_name] += len(
                            inter
                        )"""

            """counts_by_col2 = {
                c: len(new_row[c]) for c in ["auto", "exp", "phylo", "negative"]
            }
            total2 = sum(counts_by_col2.values())

            if total2 < total1:
                print(f"{row['id']}")
                print(f"  Before: {counts_by_col}")
                print(f"  After: {counts_by_col2}")"""

            for col in all_evi_groups + ["derived_not"]:
                new_row[col] = sorted(list(new_row[col]))
                annot_counts[ont][col] += len(new_row[col])

            new_rows[ont].append(new_row)
            index = len(new_rows[ont]) - 1
            uniprot_to_index[ont][row["id"]] = index

    print("Adding derivated negatives of unmentioned proteins to row lists")
    for ont, by_uniprot in derivated_negatives.items():
        for uniprot, set_of_goids in by_uniprot.items():
            new_row = {
                "id": uniprot,
                "ProteinSet": "SwissProt",
                "derived_not": sorted(list(set_of_goids)),
            }
            for col in all_evi_groups:
                new_row[col] = []
            annot_counts[ont]["derived_not"] += len(new_row["derived_not"])
            new_rows[ont].append(new_row)
            index = len(new_rows[ont]) - 1
            uniprot_to_index[ont][uniprot] = index

    print("Annotation counts:")
    print(json.dumps(annot_counts, indent=4))

    for ont in new_rows.keys():
        print("Creating dataframe for", ont)
        df = pl.DataFrame(new_rows[ont])
        df.write_parquet(output_dir + "/go." + ont + ".parquet")

    json.dump(annot_counts, open(output_dir + "/go.annot_counts.json", "w"), indent=4)
    json.dump(
        redundancy_counts,
        open(output_dir + "/go.redundancy_counts.json", "w"),
        indent=4,
    )
