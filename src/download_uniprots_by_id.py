import requests
import time
import sys

input_file_path = sys.argv[1]     # seu arquivo com IDs (um por linha)
output_file_path = sys.argv[2]     # nome do arquivo de saída
batch_size = 99                   # tamanho de cada lote de requisição
sleep_time = 1                     # intervalo entre requisições (segundos)

BASE_URL = "https://rest.uniprot.org/uniprotkb/search"

def read_ids(filename):
    with open(filename) as f:
        ids = [line.strip() for line in f if line.strip()]
    return ids

def fetch_fasta_batch(id_batch):
    query = " OR ".join(id_batch)
    params = {
        "query": query,
        "format": "fasta",
        "size": len(id_batch)
    }
    response = requests.get(BASE_URL, params=params)
    if response.status_code != 200:
        print(f"⚠️ Erro {response.status_code}: {response.text[:200]}")
        return ""
    return response.text

def main():
    all_ids = read_ids(input_file_path)
    print(f"Total de IDs: {len(all_ids)}")
    with open(output_file_path, "w") as out:
        for i in range(0, len(all_ids), batch_size):
            batch = all_ids[i:i+batch_size]
            print(f"Baixando {i+1}-{i+len(batch)} ...")
            fasta_data = fetch_fasta_batch(batch)
            out.write(fasta_data)
            time.sleep(sleep_time)
    print("✅ Download concluído!")

if __name__ == "__main__":
    main()
