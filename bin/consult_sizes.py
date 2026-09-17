import json
from huggingface_hub import HfApi
from huggingface_hub.utils import RepositoryNotFoundError

def get_model_sizes(model_names):
    api = HfApi()
    model_sizes = {}

    print(f"{'Model Name':<35} | {'Parameters (Millions)':<20}")
    print("-" * 60)

    for model_ids in model_names:
        model_size_millions = None
        for model_id in model_ids:
            if model_size_millions is None:
                try:
                    # Fetch only the model metadata from the API (very fast, no weights downloaded)
                    info = api.model_info(model_id)
                    
                    # Check if safetensors metadata is available
                    if info.safetensors and info.safetensors.parameters:
                        # info.safetensors.parameters is a dict like {'F16': 7000000000, 'I8': ...}
                        total_params = sum(info.safetensors.parameters.values())
                        model_size_millions = total_params / 1_000_000
                        print(f"{model_id:<35} | {model_size_millions:,.2f} M")
                    else:
                        print(f"{model_id:<35} | N/A (No safetensors metadata)")
                        
                except RepositoryNotFoundError:
                    print(f"{model_id:<35} | Error: Model not found")
                except Exception as e:
                    print(f"{model_id:<35} | Error: {str(e)}")
        model_sizes[model_ids[0]] = model_size_millions
    return model_sizes

models = [
    ["ElnaggarLab/ankh-base", "Synthyra/ANKH_base"],
    ["ElnaggarLab/ankh-large", "Synthyra/ANKH_large"],
    ["ElnaggarLab/ankh2-ext2", "Synthyra/ANKH2_large"],
    ["ElnaggarLab/ankh3-large", "Synthyra/ANKH3_Large"],
    ["Profluent-Bio/E1-150m"],
    ["Profluent-Bio/E1-300m"],
    ["Profluent-Bio/E1-600m"],
    ["facebook/esm2_t30_150M_UR50D"],
    ["facebook/esm2_t33_650M_UR50D",],
    ["facebook/esm2_t36_3B_UR50D", "Synthyra/ESM2-3B"],
    ["biohub/ESMC-300M-hf", "biohub/ESMC-600M"],
    ["biohub/ESMC-600M-hf", "biohub/ESMC-600M"],
    ["flair-bio/amplify-120m"],
    ["flair-bio/amplify-350m"],
]

sizes = get_model_sizes(models)

with open("input_data/model_sizes.json", "w") as f:
    json.dump(sizes, f, indent=4, ensure_ascii=False)
    