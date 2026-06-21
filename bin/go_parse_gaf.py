#!/usr/bin/env python3
import sys
import polars as pl
from tqdm import tqdm

from bioinfo_utils.util_base import count_lines_large, open_file


def parse_gaf(goa_gaf, bar, uniprots_set):
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
    rels = []
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
                    if cells[1] in uniprots_set:
                        protid = cells[1]
                        relation = cells[3]
                        goid = cells[4]
                        evi = cells[6]
                        taxid = cells[12]
                        if len(evi) < 2 or len(evi) > 3:
                            print("strange evidence:", evi)
                        else:
                            for tx in taxid.split("|"):
                                protids.append(protid)
                                rels.append(relation)
                                goids.append(goid)
                                evis.append(evi)
                                taxids.append(tx)
                                aspects.append(cells[8])
                            # parsed.append('\t'.join([protid,goid,evi,taxid,cells[8]]))
                    else:
                        not_swissprot += 1

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
            "relation": rels,
            "goid": goids,
            "evi": evis,
            "taxonid": taxids,
            "aspect": aspects,
        }
    )

    return df, goa_parsing_report


go_annotation_raw = sys.argv[1]
ids_path = sys.argv[2]
output_dir = "./"
goa_parsed = output_dir + "/goa.uniprot.parquet"

uniprots_list = open(ids_path).read().split("\n")
uniprots_set = set(uniprots_list)

print("Counting length of ", go_annotation_raw)
count_ann_total_lines = count_lines_large(go_annotation_raw)
print(count_ann_total_lines, "lines of annotation")

print("opening", go_annotation_raw)
bar = tqdm(total=count_ann_total_lines)
goa_gaf = open_file(go_annotation_raw)
goa_experimental_df, goa_parsing_report = parse_gaf(goa_gaf, bar, uniprots_set)
report_path = goa_parsed.replace(".parquet", ".log")
open(report_path, "w").write(goa_parsing_report)
goa_experimental_df.write_parquet(goa_parsed)
