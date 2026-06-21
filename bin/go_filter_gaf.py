#!/usr/bin/env python
import gzip
from os import path
import sys
import pandas as pd
from tqdm import tqdm
import polars as pl
import json

from bioinfo_utils.gene_ontology import (
    gos_not_to_use,
    load_go_graph_strict,
    expand_go_list,
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

cafa_evis = [
    "EXP",
    "IDA",
    "IPI",
    "IMP",
    "IGI",
    "IEP",
    "HTP",
    "HDA",
    "HMP",
    "HGI",
    "HEP",
    "TAS",
    "IC",
]
phylo_evis = ["IBA", "IBD", "IKR", "IRD"]
auto_evis = ["ISS", "ISO", "ISA", "ISM", "IGP", "RCA", "IEA"]


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
                "mf_exp": set(),
                "bp_exp": set(),
                "cc_exp": set(),
                "mf_phylo": set(),
                "bp_phylo": set(),
                "cc_phylo": set(),
                "mf_auto": set(),
                "bp_auto": set(),
                "cc_auto": set(),
                "mf_negative": set(),
                "bp_negative": set(),
                "cc_negative": set(),
            }
        prot_ann = indexed_ann[protid]
        if onto in ["F", "P", "C"]:
            prefix = onto_to_prefix[onto] + "_"
            if evi in cafa_evis or evi in phylo_evis:
                if "NOT" in rel_str:
                    prot_ann[prefix + "negative"].add(goid)
                else:
                    if evi in cafa_evis:
                        prot_ann[prefix + "exp"].add(goid)
                    else:
                        prot_ann[prefix + "phylo"].add(goid)
            elif evi in auto_evis:
                prot_ann[prefix + "auto"].add(goid)

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
    output_dir = "./"

    goa_by_uniprot_path = output_dir + "/go.by_uniprot.parquet"

    print("Reading swissprot ids")
    uniprots_list = open(ids_path).read().split("\n")

    uniprots_set = set(uniprots_list)

    goa_parsed = pl.read_parquet(goa_parsed)

    goa_by_uniprot = gafparsed_to_id2go(goa_parsed, uniprots_set, goa_by_uniprot_path)

    goes_to_not_use = gos_not_to_use(go_not_use_path)
    go_graph = load_go_graph_strict(go_basic_path)
    new_parsed = set()

    new_rows = {"mf": [], "bp": [], "cc": []}

    annot_counts = {
        "mf": {"auto": 0, "exp": 0, "phylo": 0, "negative": 0},
        "bp": {"auto": 0, "exp": 0, "phylo": 0, "negative": 0},
        "cc": {"auto": 0, "exp": 0, "phylo": 0, "negative": 0},
    }

    print("Expanding GOs")
    for row in tqdm(goa_by_uniprot.rows(named=True), total=goa_by_uniprot.height):
        for ont in ["mf", "bp", "cc"]:
            new_row = {
                "id": row["id"],
                "ProteinSet": row["ProteinSet"],
                "exp": set(
                    expand_go_list(row[ont + "_exp"], go_graph, goes_to_not_use)
                ),
                "phylo": set(
                    expand_go_list(row[ont + "_phylo"], go_graph, goes_to_not_use)
                ),
                "auto": set(
                    expand_go_list(row[ont + "_auto"], go_graph, goes_to_not_use)
                ),
                "negative": set(
                    expand_go_list(
                        row[ont + "_negative"],
                        go_graph,
                        goes_to_not_use,
                        is_negative=True,
                    )
                ),
            }

            """counts_by_col = {
                c: len(new_row[c]) for c in ["auto", "exp", "phylo", "negative"]
            }
            total1 = sum(counts_by_col.values())"""

            to_remove = [
                {"go_col": "negative", "from": ["exp", "auto", "phylo"]},
                {"go_col": "exp", "from": ["auto", "phylo"]},
                {"go_col": "phylo", "from": ["auto"]},
            ]

            for removal in to_remove:
                source_set = new_row[removal["go_col"]]
                for other_set_name in removal["from"]:
                    new_row[other_set_name] -= source_set

            """counts_by_col2 = {
                c: len(new_row[c]) for c in ["auto", "exp", "phylo", "negative"]
            }
            total2 = sum(counts_by_col2.values())

            if total2 < total1:
                print(f"{row['id']}")
                print(f"  Before: {counts_by_col}")
                print(f"  After: {counts_by_col2}")"""

            for col in ["auto", "exp", "phylo", "negative"]:
                new_row[col] = sorted(list(new_row[col]))
                annot_counts[ont][col] += len(new_row[col])

            new_rows[ont].append(new_row)
    print("Annotation counts:")
    print(json.dumps(annot_counts, indent=4))

    for ont in new_rows.keys():
        print("Creating dataframe for", ont)
        df = pl.DataFrame(new_rows[ont])
        df.write_parquet(output_dir + "/go." + ont + ".parquet")

    json.dump(annot_counts, open(output_dir + "/go.annot_counts.json", "w"), indent=4)
