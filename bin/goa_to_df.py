import sys
import polars as pl
import os
import numpy as np
from tqdm import tqdm
import gzip
import obonet
from mf_swarm_lib.core.node import go_descendants_dict, go_ancestors_dict


def make_goa_matrix(swarm_dir, ids, goa_path, goa_scores_path):
    print("Loading GOA annotations for swarm:", swarm_dir)
    labels_path = f"{swarm_dir}/go_ids.txt"
    labels = open(labels_path).read().splitlines()

    # Create zeros dataframe with len(ids) rows and len(labels) columns with numpy zeros
    # Create ids index: id  -> index and labels index
    # Run through every line of goa scores path, saving the score in the matrix
    # Run through every line of goa path, finding negative statements and assigning them a score of -2
    # Make numpy df. Columns: ids, scores (goa matrix)
    # Save as parquet

    """
    goa_submission
    A0A009IHW8\tGO:0016787\t1.0
    A0A009IHW8\tGO:0050135\t1.0
    A0A009IHW8\tGO:0003824\t1.0
    A0A009IHW8\tGO:0016798\t1.0
    A0A009IHW8\tGO:0003674\t1.0
    A0A009IHW8\tGO:0016799\t1.0
    A0A009IHW8\tGO:0008150\t1.0
    A0A009IHW8\tGO:0005575\t1.0
    A0A009IHW8\tGO:0061809\t0.999991493
    A0A009IHW8\tGO:0019677\t0.99962947
    A0A009IHW8\tGO:0008152\t0.99962947

    goa_raw:
    What is a GAF?
    GAF (Gene Association File) is a tab-separated standard used to distribute GO annotations. 
        Each row links a database object (e.g., a UniProt accession) to a GO term, alongside evidence and metadata. 
        In GAF v2.2, the key fields are:
            DB_Object_ID (column 2): the stable identifier for the entity (e.g., UniProt ID like P12345).
            Qualifier (column 4): may include relationship labels (e.g., enables, involved_in) and negations using NOT.
            GO_ID (column 5): the Gene Ontology term (e.g., GO:0005515).
        Comment lines start with ! and should be ignored.
        Why care about NOT? A NOT qualifier explicitly states that the gene/product is not associated 
            with that term/relationship. For leaderboard alignment, I removed these negated associations 
            before building the submission.
    Bery large TSV, must scan line by line
    """
    id_index = {id: i for i, id in enumerate(ids)}
    # Create labels index: label  -> index
    label_index = {label: i for i, label in enumerate(labels)}
    goa_matrix = np.zeros((len(ids), len(labels)))
    labels_used = set()
    ids_used = set()
    bar = tqdm(total=len(ids))
    for line in open(goa_scores_path):
        idx, label, score = line.strip().split("\t")
        if idx in id_index and label in label_index:
            goa_matrix[id_index[idx]][label_index[label]] = float(score)
            if not idx in ids_used:
                bar.update(1)
            ids_used.add(idx)
            labels_used.add(label)
    bar.close()

    print(f"Ids used: {len(ids_used)}")
    print(f"Labels used: {len(labels_used)}")

    # Scan GOA
    labels_used = set()
    ids_used = set()
    bar = tqdm(total=len(ids))
    input_stream = gzip.open(goa_path, "rt") if ".gz" in goa_path else open(goa_path)
    for rawline in input_stream:
        if rawline.startswith("!"):
            continue
        cells = rawline.strip().split("\t")
        uniprotid = cells[1].split("|")[-1]
        qualifier = cells[3]
        go_id = cells[4]

        if uniprotid in id_index and "NOT|" in qualifier:
            # print(uniprotid, qualifier, go_id)
            if go_id in label_index:
                prot_idx = id_index[uniprotid]
                label_idx = label_index[go_id]
                goa_matrix[prot_idx][label_idx] = -2
                if not uniprotid in ids_used:
                    bar.update(1)
                ids_used.add(uniprotid)
                labels_used.add(go_id)
    bar.close()

    print(f"Ids used: {len(ids_used)}")
    print(f"Labels used: {len(labels_used)}")

    return goa_matrix, labels


def propagate_goa(goa_matrix, obo_path, go_id_sequence):
    # Propagates positive scores up and negative scores down the GO hierarchy
    print("Loading GO network")
    go_network = obonet.read_obo(obo_path)
    print("Calculating GO descendants and ancestors")
    _, go_descendants_indexes = go_descendants_dict(
        go_id_sequence, go_network, only_in_sequence=True
    )
    _, go_ancestors_indexes = go_ancestors_dict(
        go_id_sequence, go_network, only_in_sequence=True
    )

    print("Propagating GOA scores")
    for prot_idx in tqdm(range(len(goa_matrix))):
        scores_row = goa_matrix[prot_idx]
        positive_cols = [n for n, x in enumerate(scores_row) if x > 0]
        negative_cols = [n for n, x in enumerate(scores_row) if x < 0]

        for col_idx in positive_cols:
            col_name = go_id_sequence[col_idx]
            ancestors_list = go_ancestors_indexes[col_name]
            for ancestor_idx in ancestors_list:
                if scores_row[ancestor_idx] < scores_row[col_idx]:
                    scores_row[ancestor_idx] = scores_row[col_idx]

        neg_indexes = set()
        for neg_idx in negative_cols:
            col_name = go_id_sequence[neg_idx]
            descendants_list = go_descendants_indexes[col_name]
            neg_indexes.update(descendants_list)
        for neg_idx in neg_indexes:
            scores_row[neg_idx] = -2

    return goa_matrix


if __name__ == "__main__":
    mf_swarm_dir = sys.argv[1]
    bp_swarm_dir = sys.argv[2]
    cc_swarm_dir = sys.argv[3]
    cafa6_dir = sys.argv[4]
    output_path = sys.argv[5]

    goa_path = f"{cafa6_dir}/additional_features/goa_uniprot_all-not.gaf"
    goa_scores_path = f"{cafa6_dir}/additional_features/goa_submission.tsv"
    protein_ids_path = f"{cafa6_dir}/ids.test.txt"
    go_basic_path = f"{cafa6_dir}/go-basic.obo"
    ids = open(protein_ids_path).read().splitlines()

    mf_matrix, mf_labels = make_goa_matrix(mf_swarm_dir, ids, goa_path, goa_scores_path)
    mf_matrix = propagate_goa(mf_matrix, go_basic_path, mf_labels)
    bp_matrix, bp_labels = make_goa_matrix(bp_swarm_dir, ids, goa_path, goa_scores_path)
    bp_matrix = propagate_goa(bp_matrix, go_basic_path, bp_labels)
    cc_matrix, cc_labels = make_goa_matrix(cc_swarm_dir, ids, goa_path, goa_scores_path)
    cc_matrix = propagate_goa(cc_matrix, go_basic_path, cc_labels)

    df = pl.DataFrame({"id": ids, "mf": mf_matrix, "bp": bp_matrix, "cc": cc_matrix})
    df.write_parquet(output_path)
    open(output_path + ".labels.mf.txt", "w").write("\n".join(mf_labels))
    open(output_path + ".labels.bp.txt", "w").write("\n".join(bp_labels))
    open(output_path + ".labels.cc.txt", "w").write("\n".join(cc_labels))
