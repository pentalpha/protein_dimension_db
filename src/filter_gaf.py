import gzip
from os import path
import sys
import pandas as pd
from tqdm import tqdm
import polars as pl

from bioinfo_utils.gene_ontology import expand_go_set, gos_not_to_use, load_go_graph
from bioinfo_utils.util_base import count_lines_large, open_file, write_file

# manual urls:
# quickgo:
# https://www.ebi.ac.uk/QuickGO/annotations?aspect=molecular_function&evidenceCode=ECO:0000352,ECO:0000269,ECO:0000314,ECO:0000315,ECO:0000316,ECO:0000353,ECO:0000270,ECO:0007005,ECO:0007001,ECO:0007003,ECO:0007007,ECO:0006056,ECO:0000318,ECO:0000320,ECO:0000321,ECO:0000304,ECO:0000305&evidenceCodeUsage=descendants&withFrom=UniProtKB
# geneontology.org
# http://current.geneontology.org/annotations/filtered_goa_uniprot_all_noiea.gaf.gz


def gafparsed_to_id2go(go_expanded_ann: pl.DataFrame, swissprot_ids, output_file):
    indexed_ann = {}

    print("Indexing by protein id")
    for row in tqdm(go_expanded_ann.rows(named=True), total=go_expanded_ann.height):
        protid = row["prot_id"]
        goid = row["goid"]
        evi = row["evi"]
        taxid = row["taxonid"]
        onto = row["aspect"]
        if protid not in indexed_ann:
            indexed_ann[protid] = {
                "id": protid,
                "mf": [],
                "bp": [],
                "cc": [],
                "mf_evi": [],
                "bp_evi": [],
                "cc_evi": [],
            }
        if onto in ["F", "P", "C"]:
            if onto == "F":
                ann_col = "mf"
                evi_col = "mf_evi"
            elif onto == "P":
                ann_col = "bp"
                evi_col = "bp_evi"
            else:
                ann_col = "cc"
                evi_col = "cc_evi"
            if goid not in indexed_ann[protid][ann_col]:
                indexed_ann[protid][ann_col].append(goid)
                indexed_ann[protid][evi_col].append(evi)

    lines = [v for k, v in indexed_ann.items()]
    for l in lines:
        if l["id"] in swissprot_ids:
            l["ProteinSet"] = "SwissProt"
        else:
            l["ProteinSet"] = "TrEMBL"
    print("Creating final dataframe")
    df = pl.DataFrame(lines)
    df.write_parquet(output_file)


def parse_gaf(goa_gaf, evi_not_use, bar, uniprots_set):
    parsed = ["prot_id\tgoid\tevi\ttaxonid\taspect"]
    droped = 0
    other_ontos = 0
    quickgolines = 0
    other_dbs = 0
    incorrect_line_number = 0
    other_dbs_names = set()
    not_swissprot = 1
    # cols used: 1, protein_id, go_id, evi_type, taxid, 8

    protids = []
    goids = []
    evis = []
    taxids = []
    aspects = []

    try:
        for line in goa_gaf:
            quickgolines += 1
            if line.startswith("UniProtKB"):
                cells = line.rstrip("\n").split("\t")
                if len(cells) >= 13:
                    if cells[1] not in uniprots_set:
                        not_swissprot += 1
                    protid = cells[1]
                    goid = cells[4]
                    evi = cells[6]
                    taxid = cells[12]
                    if len(evi) < 2 or len(evi) > 3:
                        print("strange evidence:", evi)
                    else:
                        if not evi in evi_not_use:
                            for tx in taxid.split("|"):
                                protids.append(protid)
                                goids.append(goid)
                                evis.append(evi)
                                taxids.append(tx)
                                aspects.append(cells[8])
                            # parsed.append('\t'.join([protid,goid,evi,taxid,cells[8]]))
                        else:
                            droped += 1
                else:
                    incorrect_line_number += 1
            else:
                other_dbs += 1
            bar.update(1)
    except EOFError as err:
        print("GAF file download incomplete")
        print(err)
    bar.close()
    goa_parsing_report = f"""
        {quickgolines} lines in goa_uniprot_all original
        {other_dbs} from other dbs: {other_dbs_names}
        {other_ontos} from other ontologies
        {droped} with evi codes we cant use
        {not_swissprot} not swissprot, but included
        {incorrect_line_number} incorrect_line_number
        {quickgolines - other_dbs - other_ontos - droped - incorrect_line_number} final lines
        {len(protids)} final lines in protids
    """
    print(goa_parsing_report)

    df = pl.DataFrame(
        {
            "prot_id": protids,
            "goid": goids,
            "evi": evis,
            "taxonid": taxids,
            "aspect": aspects,
        }
    )

    return df, goa_parsing_report


if __name__ == "__main__":
    evi_not_use_path = sys.argv[1]
    go_not_use_path = sys.argv[2]
    go_basic_path = sys.argv[3]
    go_annotation_raw = sys.argv[4]
    ids_path = sys.argv[5]
    output_dir = "./"
    calc_go_expanded = sys.argv[6] == "rerun"

    # evi_not_use_path = proj_dir+'/evi_not_to_use.txt'
    # go_not_use_path = proj_dir+"/databases/gocheck_do_not_annotate.json"
    # go_basic_path = proj_dir+"/databases/go-basic.obo"

    goa_parsed = output_dir + "/go.experimental.parquet"
    goa_expanded = output_dir + "/go.experimental_expanded.parquet"
    goa_parsed_expanded_final = output_dir + "/go.by_uniprot.parquet"

    print("Reading swissprot ids")
    uniprots_list = open(ids_path).read().split("\n")

    uniprots_set = set(uniprots_list)

    if not path.exists(goa_parsed):
        print("Counting length of ", go_annotation_raw)
        count_ann_total_lines = count_lines_large(go_annotation_raw)
        print(count_ann_total_lines, "lines of annotation")

        print("opening", go_annotation_raw)
        bar = tqdm(total=count_ann_total_lines)
        goa_gaf = open_file(go_annotation_raw)
        evi_df = pd.read_csv(evi_not_use_path, sep=",")
        evi_not_use = set(evi_df["code"].tolist())
        goa_experimental_df, goa_parsing_report = parse_gaf(
            goa_gaf, evi_not_use, bar, uniprots_set
        )
        report_path = goa_parsed.replace(".parquet", ".log")
        open(report_path, "w").write(goa_parsing_report)
        goa_experimental_df.write_parquet(goa_parsed)

    if calc_go_expanded or not path.exists(goa_expanded):
        parsed = pl.read_parquet(goa_parsed)

        print("Loading GO")
        goes_to_not_use = gos_not_to_use(go_not_use_path)
        go_graph = load_go_graph(go_basic_path)
        new_parsed = set()

        print("Expanding GO")
        expanded_lines = {
            "prot_id": [],
            "goid": [],
            "evi": [],
            "taxonid": [],
            "aspect": [],
        }
        for row in tqdm(parsed.rows(named=True), total=parsed.height):
            protid = row["prot_id"]
            goid = row["goid"]
            evi = row["evi"]
            taxids = row["taxonid"]
            onto = row["aspect"]
            expanded_set = expand_go_set(goid, go_graph, goes_to_not_use)
            for goid2 in expanded_set:
                expanded_lines["prot_id"].append(protid)
                expanded_lines["goid"].append(goid2)
                expanded_lines["evi"].append(evi)
                expanded_lines["taxonid"].append(taxids)
                expanded_lines["aspect"].append(onto)

        print(len(parsed), "annotations from goa_uniprot_all")
        print(len(new_parsed), "annotations with expansion")
        print("Writing expanded parquet", goa_expanded)
        goa_expanded_df = pl.DataFrame(expanded_lines)
        goa_expanded_df = goa_expanded_df.unique()
        goa_expanded_df.write_parquet(goa_expanded)
        # parsed_all = ['\t'.join(x) for x in new_parsed]
        # write_file(goa_expanded).write('\n'.join(sorted(parsed_all)))

    print("Reading", goa_expanded)
    go_expanded_ann = pl.read_parquet(goa_expanded)
    gafparsed_to_id2go(go_expanded_ann, uniprots_set, goa_parsed_expanded_final)

    """mf_lines = [x for x in go_expanded_ann if x[-1] == 'F']
    bp_lines = [x for x in go_expanded_ann if x[-1] == 'P']
    cc_lines = [x for x in go_expanded_ann if x[-1] == 'C']

    print('Writing id2go files')
    gafparsed_to_id2go(mf_lines, uniprots_list, goa_parsed_expanded_mf)
    gafparsed_to_id2go(bp_lines, uniprots_list, goa_parsed_expanded_bp)
    gafparsed_to_id2go(cc_lines, uniprots_list, goa_parsed_expanded_cc)"""

    """mf_parsed = ['\t'.join(x).rstrip('\tF') for x in new_parsed if x[-1] == 'F']
    bp_parsed = ['\t'.join(x).rstrip('\tP') for x in new_parsed if x[-1] == 'P']
    cc_parsed = ['\t'.join(x).rstrip('\tC') for x in new_parsed if x[-1] == 'C']
    run_command(['mkdir', path.dirname(goa_parsed_expanded_mf)])
    write_file(goa_parsed_expanded_mf).write('\n'.join(sorted(mf_parsed)))
    write_file(goa_parsed_expanded_bp).write('\n'.join(sorted(bp_parsed)))
    write_file(goa_parsed_expanded_cc).write('\n'.join(sorted(cc_parsed)))"""
