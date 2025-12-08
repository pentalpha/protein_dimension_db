import polars as pl
import sys

from tqdm import tqdm

if __name__ == '__main__':
    #Convert annotation file from by protein id to by term
    go_by_uniprot_path = sys.argv[1]
    go_by_term_path = './go.by_term.parquet'

    terms = []
    term_to_index = {}
    index = 0

    input_df = pl.read_parquet(go_by_uniprot_path)
    aspects = ['mf', 'bp', 'cc']
    print('Reading go terms used')
    for row in tqdm(input_df.rows(named=True), total=input_df.height):
        for c in aspects:
            goids = row[c]
            if goids is None:
                continue
            for goid in goids:
                if goid not in term_to_index:
                    term_to_index[goid] = index
                    terms.append({
                        'id': goid,
                        'aspect': c,
                        'proteins': [],
                        'evidences': []
                    })
                    index += 1

    print('Indexing by term')
    for row in tqdm(input_df.rows(named=True), total=input_df.height):
        protid = row['id']
        for c in aspects:
            goids = row[c]
            evis = row[c+'_evi']
            if goids is None:
                continue
            for i, goid in enumerate(goids):
                term_index = term_to_index[goid]
                terms[term_index]['proteins'].append(protid)
                terms[term_index]['evidences'].append(evis[i])
    
    print('Creating final dataframe')
    df = pl.DataFrame(terms)
    df.write_parquet(go_by_term_path)