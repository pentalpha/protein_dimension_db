import json
from multiprocessing import Pool
from os import path, remove
import sys
from time import sleep

import requests


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i : i + n]


# Locked method to access file for writing, so only one process writes at a time
def open_locked_file(file_path: str):
    lock_path = file_path + ".lock"
    while path.exists(lock_path):
        sleep(0.25)
    output = open(file_path, "a")
    return output, lock_path


def close_locked_file(output, lock_path: str):
    output.close()
    if path.exists(lock_path):
        remove(lock_path)


def get_mf_annots(uniprot_id: str):
    """
    GET /api/entry/interpro/protein/uniprot/p99999
    HTTP 200 OK
    Allow: GET, HEAD
    Cached: true
    Content-Type: application/json
    InterPro-Version: 107.0
    InterPro-Version-Minor: 0
    Server-Timing:
    Vary: Accept

    {
        "count": 3,
        "next": null,
        "previous": null,
        "results": [
            {
                "metadata": {
                    "accession": "IPR002327",
                    "name": "Cytochrome c, class IA/ IB",
                    "source_database": "interpro",
                    "type": "family",
                    "integrated": null,
                    "member_databases": {
                        "panther": {
                            "PTHR11961": "CYTOCHROME C"
                        },
                        "prints": {
                            "PR00604": "CYTCHRMECIAB"
                        }
                    },
                    "go_terms": [
                        {
                            "identifier": "GO:0009055",
                            "name": "electron transfer activity",
                            "category": {
                                "code": "F",
                                "name": "molecular_function"
                            }
                        },
                        {
                            "identifier": "GO:0020037",
                            "name": "heme binding",
                            "category": {
                                "code": "F",
                                "name": "molecular_function"
                            }
                        }
                    ]
                },
                "proteins": [
                    {
                        "accession": "p99999",
                        "protein_length": 105,
                        "source_database": "reviewed",
                        "organism": "9606",
                        "in_alphafold": true,
                        "in_bfvd": false,
                        "entry_protein_locations": [
                            {
                                "fragments": [
                                    {
                                        "start": 2,
                                        "end": 103,
                                        "dc-status": "CONTINUOUS"
                                    }
                                ],
                                "representative": false,
                                "model": null,
                                "score": null
                            }
                        ]
                    }
                ]
            },
        [...]
    """

    url = f"https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/{uniprot_id}"
    headers = {"Accept": "application/json"}
    response = requests.get(url, headers=headers)
    annots = set()
    if response.status_code == 200:
        data = response.json()
        try:
            for result in data["results"]:
                if "metadata" in result:
                    metadata = result["metadata"]
                    if metadata["source_database"] == "interpro":
                        go_terms = metadata["go_terms"]
                        if type(go_terms) == list:
                            for go_term in go_terms:
                                if go_term.get("category", {}).get("code") == "F":
                                    annots.add(go_term.get("identifier"))
        except Exception as err:
            print(json.dumps(data, indent=4))
            return None

        return list(annots)
    return None


def get_interpro_families(uniprot_id: str):
    url = f"https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/{uniprot_id}"
    headers = {"Accept": "application/json"}
    response = requests.get(url, headers=headers)
    annots = set()
    if response.status_code == 200:
        data = response.json()
        try:
            for result in data["results"]:
                if "metadata" in result:
                    metadata = result["metadata"]
                    if metadata["source_database"] == "interpro":
                        fam = metadata["accession"]
                        if type(fam) == str:
                            annots.add(fam)
        except Exception as err:
            print(json.dumps(data, indent=4))
            return None

        return list(annots)
    return None


def consult_uniprot_ids(
    ids: list, cache_path: str, n_per_chunk: int = 150, info_type="mf"
):
    ids_consulted = set()
    if path.exists(cache_path):
        with open(cache_path, "r") as cache_file:
            for line in cache_file:
                cells = line.strip().split("\t")
                if len(cells) > 1:
                    ids_consulted.add(cells[0])

    ids_to_consult = [
        uniprot_id for uniprot_id in ids if uniprot_id not in ids_consulted
    ]
    print(f"Consulting {len(ids_to_consult)} uniprot ids...")
    print(f"{len(ids_consulted)} ids already in cache.")

    chunk_list = list(chunks(ids_to_consult, n_per_chunk))
    print(f"Processing {len(chunk_list)} chunks of up to {n_per_chunk} ids each.")

    if info_type == "mf":
        consult_func = get_mf_annots
    elif info_type == "family":
        consult_func = get_interpro_families
    else:
        raise ValueError(f"Invalid info_type: {info_type}")

    for i, chunk in enumerate(chunk_list):
        perc = (i + 1) / len(chunk_list) * n_per_chunk
        print(f"Processing at {perc:.2f}% ({i+1}/{len(chunk_list)})...")
        with Pool(processes=10) as pool:
            results_list = pool.map(consult_func, chunk)
        results = {
            uniprot_id: annots for uniprot_id, annots in zip(chunk, results_list)
        }
        sleep(0.5)

        output, lock_path = open_locked_file(cache_path)
        for uniprot_id, annots in results.items():
            if annots != None:
                annot_str = ";".join(annots)
                output.write(f"{uniprot_id}\t{annot_str}\n")
            else:
                output.write(f"{uniprot_id}\t\n")
        close_locked_file(output, lock_path)
