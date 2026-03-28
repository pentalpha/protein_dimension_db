def parse_interpro_raw(input_file: str, output_file: str):
    """Raw interpro tsv format:
    1. Protein accession (e.g. P51587)
    2. Sequence MD5 digest (e.g. 14086411a2cdf1c4cba63020e1622579)
    3. Sequence length (e.g. 3418)
    4. Analysis (e.g. Pfam / PRINTS / Gene3D)
    5. Signature accession (e.g. PF09103 / G3DSA:2.40.50.140)
    6. Signature description (e.g. BRCA2 repeat profile)
    7. Start location
    8. Stop location
    9. Score - is the e-value (or score) of the match reported by member database method (e.g. 3.1E-52)
    10. Status - is the status of the match (T: true)
    11. Date - is the date of the run
    12. InterPro annotations - accession (e.g. IPR002093)
    13. InterPro annotations - description (e.g. BRCA2 repeat)
    14. GO annotations with their source(s), e.g. GO:0005515(InterPro)|GO:0006302(PANTHER)|GO:0007195(InterPro,PANTHER). This is an optional column; only displayed if the --goterms option is switched on
    15. Pathways annotations, e.g. REACT_71. This is an optional column; only displayed if the --pathways option is switched on
    """

    col_names = [
        "protein_accession",
        "sequence_md5",
        "sequence_length",
        "analysis",
        "signature_accession",
        "signature_description",
        "start_location",
        "stop_location",
        "score",
        "status",
        "date",
        "accession",
        "desc",
        "go_annotations",
        "pathways_annotations",
    ]

    interpro_annotations = {}
    for rawline in open(input_file, "r"):
        line = rawline.strip().split("\t")
        row = {col_names[i]: line[i] for i in range(len(col_names))}
        if not row["protein_accession"] in interpro_annotations:
            interpro_annotations[row["protein_accession"]] = set()
        if len(row["accession"]) > 4:
            interpro_annotations[row["protein_accession"]].add(row["accession"])

    with open(output_file, "w") as out_file:
        for protein_accession, annotations in interpro_annotations.items():
            annots_list = ";".join(sorted(annotations))
            out_file.write(f"{protein_accession}\t{annots_list}\n")
