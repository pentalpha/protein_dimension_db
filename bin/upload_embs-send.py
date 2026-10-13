import os
import sys
import glob
from huggingface_hub import HfApi, DatasetCard, DatasetCardData

basecard = r"""---
{{ card_data }}
---
# Protein Dimension DB 2.0 - [PLM_NAME] - Protein Embeddings

This dataset provides pre-computed sequence embeddings for all proteins in the UniprotKB/SwissProt database, generated using the **[PLM_NAME]** Protein Language Model (PLM). 

These embeddings represent a comprehensive structural and functional vector space for proteins, meant to be used for downstream machine learning tasks, representation learning, and functional analysis. 
This dataset is part of the Protein Dimension DB 2.0 release.

Poolings available:
- **PoolParti**: Partitioned pooling
- **Mean**: Average pooling over the sequence length
- **Max**: Max pooling over the sequence length
- **STD**: Standard deviation pooling over the sequence length

## How to use this dataset

Because the dataset is heavily partitioned, you can stream or download it using the Hugging Face `datasets` library, or load it directly into Polars for high-performance, out-of-core data manipulation.

```python
# --- 1. Using the datasets library ---
from datasets import load_dataset

# Load the entire dataset (will download the partitioned parquet files)
dataset = load_dataset("pitagoras-alves/[REPO_NAME]", split="train")
print(dataset[0])

# --- 2. Using Polars to read directly from the Hub (Recommended) ---
import polars as pl

# Lazy load the dataset directly from Hugging Face without downloading all files
lf = pl.scan_parquet("hf://datasets/pitagoras-alves/[REPO_NAME]/*.parquet")

# Filter for specific proteins and poolings before collecting into memory
df = (lf.filter(pl.col("id").is_in(["P12345", "O60260"]))
        .select(["id", "mean", "max"])
        .collect())
print(df)
```

## Citation

New manuscript being written. For now, cite the previous release:

@inproceedings{alves2025protein,
  title        = {Protein Dimension DB: A Unified Protein Repository for Representation Learning and Functional Analysis},
  author       = {Alves Sobrinho, Pitágoras de Azevedo and Sakamoto, Tetsu and Figuerola, Wilfredo Blanco},
  booktitle    = {Bioinformatics and Computational Biology: 21st Brazilian Conference, X-Meeting 2025, João Pessoa, Brazil, June 3--6, 2025, Proceedings},
  year         = {2025},
  organization = {Springer}
}

## References

- Fournier, Q.; Vernon, R.M.; Van Der Sloot, A.; Schulz, B.; Chandar, S.; Langmead, C.J. Protein language models: is scaling necessary? BioRxiv 2024, pp. 2024–09;
- Elnaggar, A.; Essam, H.; Salah-Eldin, W.; Moustafa, W.; Elkerdawy, M.; Rochereau, C.; Rost, B. Ankh: Optimized protein language model unlocks general-purpose modelling. BioRxiv 2023, pp. 2023–01;
- Heinzinger, M.; Weissenow, K.; Sanchez, J.G.; Henkel, A.; Mirdita, M.; Steinegger, M.; Rost, B. Bilingual language model for protein sequence and structure. NAR Genomics and Bioinformatics 2024, 6, lqae150;
- Lin, Z.; Akin, H.; Rao, R.; Hie, B.; Zhu, Z.; Lu, W.; Smetanin, N.; Verkuil, R.; Kabeli, O.; Shmueli, Y.; et al. Evolutionary-scale prediction of atomic-level protein structure with a language model. Science 2023, 379, 1123–1130;
- Alsamkary, H.; Elshaffei, M.; Elkerdawy, M.; Elnaggar, A. Ankh3: Multi-task pretraining with sequence denoising and completion enhances protein representations. arXiv preprint arXiv:2505.20052 2025;
- Jain, S.; Beazer, J.; Ruffolo, J.A.; Bhatnagar, A.; Madani, A. E1: Retrieval-augmented protein encoder models. bioRxiv 2025, pp. 2025–11;

"""

if os.path.exists(".env"):
    for rawline in open(".env", "r"):
        line = rawline.strip()
        if line.startswith("#") or line == "":
            continue
        parts = line.split("=", 1)
        if len(parts) == 2:
            os.environ[parts[0]] = parts[1]

USER_NAME = "pitagoras-alves"
HF_TOKEN = os.getenv("HF_TOKEN")

if not HF_TOKEN:
    print("Error: HF_TOKEN not found in environment variables or .env file.")
    sys.exit(1)

if __name__ == "__main__":
    api = HfApi(token=HF_TOKEN)

    if len(sys.argv) >= 2:
        DATA_DIR = sys.argv[1]
    else:
        DATA_DIR = "."
    
    # Discover all directories that contain chunked parquet files
    # (e.g., ankh3_large, amplify_120, etc.)
    all_items = os.listdir(DATA_DIR)
    model_dirs = []

    for item in all_items:
        full_path = os.path.join(DATA_DIR, item)
        if os.path.isdir(full_path):
            parquets = glob.glob(os.path.join(full_path, "train-*.parquet"))
            if len(parquets) > 0:
                model_dirs.append(item)
                
    model_dirs.sort()

    if not model_dirs:
        print(f"No partitioned model directories found in {DATA_DIR}.")
        sys.exit(0)
        
    print(f"Found {len(model_dirs)} models to upload: {model_dirs}\n")

    for model_name in model_dirs:
        print(f"=====================================================")
        print(f"Starting upload process for model: {model_name}")
        
        repo_name = f"{model_name}-embeddings-uniprot_swissprot-pddb2"
        repo_id = f"{USER_NAME}/{repo_name}"
        model_path = os.path.join(DATA_DIR, model_name)
        
        # 1. Create Repository (if it doesn't exist)
        print(f"  [1/3] Ensuring repository exists: {repo_id}")
        api.create_repo(
            repo_id=repo_id,
            repo_type="dataset",
            private=False,
            exist_ok=True
        )
        
        # 2. Upload the Pre-chunked Parquet files
        # We use upload_folder because it pushes the parquets as-is, which is massively 
        # faster than letting the `datasets` library reload and re-chunk them.
        print(f"  [2/3] Uploading partitioned parquet files from {model_path}/...")
        api.upload_folder(
            folder_path=model_path,
            repo_id=repo_id,
            repo_type="dataset",
            allow_patterns="*.parquet"
        )
        
        # 3. Create and push the Dataset Card
        print(f"  [3/3] Generating and pushing Dataset Card...")
        readme = basecard.replace("[PLM_NAME]", model_name).replace("[REPO_NAME]", repo_name)
        
        # Extract base PLM family for tags (e.g., "ankh3_large" -> "ankh3")
        plm_family = model_name.split("_")[0]
        tags = [
            "protein", 
            "embeddings", 
            "language-model", 
            "uniprot", 
            "swissprot",
            "pddb",
            plm_family
        ]
        
        card_data = DatasetCardData(
            language="en",
            license="mit",
            task_categories=["feature-extraction", "tabular-classification"],
            tags=tags,
        )
        
        template_path = "dataset_card_temp.md"
        with open(template_path, "w", encoding="utf-8") as f:
            f.write(readme)
            
        card = DatasetCard.from_template(
            card_data=card_data, 
            template_path=template_path
        )
        
        card.push_to_hub(repo_id, token=HF_TOKEN)
        
        print(f"Upload complete for {model_name}! View at: https://huggingface.co/datasets/{repo_id}\n")

    # Clean up the temporary template file
    if os.path.exists("dataset_card_temp.md"):
        os.remove("dataset_card_temp.md")
        
    print("All models uploaded successfully!")