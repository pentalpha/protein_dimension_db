import sys
from glob import glob
import os

from tqdm import tqdm
import polars as pl
import numpy as np

allowed_models = [
    "prottrans",
    "ankh_large",
    "ankh_base",
    "esm2_t12",
    "esm2_t30",
    "esm2_t33",
    "esm2_t36",
]


def make_joined_df(emb_paths, id_sets, original_order_ids, output_path):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    print("Joining dataframes:", emb_paths, "to", output_path)

    print("Mapping from which emb file to load proteins")
    place_to_load_from = []
    for uniprot_id in original_order_ids:
        emb_indexes = [n for n, id_set in enumerate(id_sets) if uniprot_id in id_set]
        if len(emb_indexes) == 0:
            place_to_load_from.append(None)
        else:
            place_to_load_from.append(emb_indexes[0])

    id_sets2 = []
    for n, id_set in enumerate(id_sets):
        subset = [
            original_order_ids[id_n]
            for id_n, emb_path_index in enumerate(place_to_load_from)
            if emb_path_index == n
        ]
        id_sets2.append(subset)

    print("Selectively loading embeddings")
    loaded_dfs = {}
    for n, id_set in tqdm(enumerate(id_sets2)):
        if len(id_set) == 0:
            continue
        emb_path = emb_paths[n]
        print(emb_path)
        loaded_dfs[emb_path] = (
            pl.scan_parquet(emb_path).filter(pl.col("id").is_in(id_set)).collect()
        )
        print(loaded_dfs[emb_path])

    print("Joining embeddings")
    all_embs = []
    for id_n, uniprot_id in tqdm(enumerate(original_order_ids)):
        emb_path_index = place_to_load_from[id_n]
        if emb_path_index is None:
            all_embs.append(None)
        else:
            emb_path = emb_paths[emb_path_index]
            correct_df = loaded_dfs[emb_path]
            embs = correct_df.filter(pl.col("id") == uniprot_id)["emb"].to_list()
            # emb = np.array(embs[0])
            all_embs.append(embs[0])

    print("finishing")
    emb_size = next(len(e) for e in all_embs if e is not None)
    schema = {
        "id": original_order_ids,
        "emb": pl.Series(all_embs).cast(pl.Array(pl.Float64, emb_size)),
    }
    df = pl.DataFrame(schema)
    print(df)
    df.write_parquet(output_path)


if __name__ == "__main__":
    protein_ids_path = sys.argv[1]
    joined_dfs_dir = sys.argv[2]
    release_dir_paths = sys.argv[3:]

    ids = open(protein_ids_path, "r").read().strip().split("\n")
    emb_paths = set()
    for release_dir_path in release_dir_paths:
        new_paths = glob(f"{release_dir_path}/emb.*.parquet")
        new_paths = [p for p in new_paths if not "taxa_profile" in p]
        emb_paths.update(new_paths)

    emb_paths = list(emb_paths)
    emb_paths.sort()

    emb_techs = []
    for n, p in enumerate(emb_paths):
        print(p)
        emb_techs.append(os.path.basename(p).split(".")[1])

    id_sets = []
    for n, p in enumerate(emb_paths):
        df = pl.scan_parquet(p).filter(pl.col("id").is_in(ids)).collect()
        id_set = set()
        for row in df.iter_rows(named=True):
            emb = row["emb"]
            if emb == emb and emb is not None:
                id_set.add(row["id"])
        id_set = id_set.intersection(set(ids))
        id_sets.append(id_set)
        print(p, ":", len(id_set), "proteins")
        del df

    unique_emb_types = set(emb_techs)
    for emb_type in unique_emb_types:
        if emb_type not in allowed_models:
            print("Skipping", emb_type, "not configured for it.")
            continue
        emb_paths_of_type = [p for p, t in zip(emb_paths, emb_techs) if t == emb_type]
        id_sets_of_type = [id_sets[n] for n, t in enumerate(emb_techs) if t == emb_type]
        make_joined_df(
            emb_paths_of_type,
            id_sets_of_type,
            ids,
            f"{joined_dfs_dir}/emb.{emb_type}.parquet",
        )
