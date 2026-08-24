import os
os.environ["USE_FLASH_ATTN"] = "0"
import numpy as np

import torch
from E1.batch_preparer import E1BatchPreparer
from E1.modeling import E1ForMaskedLM

def unpad(seq_original_lens, batch, embeddings_np, pooled_attention_np):
    single_special_tokens = set([0])
    tuple_special_tokens = set([(1,6), (7,2)])

    input_id_list = batch["input_ids"].cpu().tolist()
    unpadded_embeddings = []
    unpadded_attentions = []

    for i in range(len(seq_original_lens)):
        seq_ids = input_id_list[i]
        residue_idx = []
        
        idx = 0
        while idx < len(seq_ids):
            # 1. Check for multi-token special boundaries (e.g., (1, 6) or (7, 2))
            matched_tuple = False
            for special_tuple in tuple_special_tokens:
                t_len = len(special_tuple)
                # Look ahead to see if the next tokens match the tuple exactly
                if idx + t_len <= len(seq_ids) and tuple(seq_ids[idx:idx+t_len]) == special_tuple:
                    idx += t_len  # Skip over the entire boundary sequence
                    matched_tuple = True
                    break
                    
            if matched_tuple:
                continue
                
            # 2. Check for single-token special characters (e.g., 0 for PAD)
            if seq_ids[idx] in single_special_tokens:
                idx += 1
                continue
                
            # 3. If it's not a special token, it's a valid biological residue!
            residue_idx.append(idx)
            idx += 1
            
        # Sanity check to ensure our mask perfectly matches the input sequence length
        # Note: If your original sequence contains commas (which aren't tokenized), 
        # seq_original_lens[i] might need to be len(seqs[i].replace(",", ""))
        assert len(residue_idx) == seq_original_lens[i], f"Mask length mismatch at seq {i}: expected {seq_original_lens[i]}, got {len(residue_idx)}"
        
        # Slicing the 2D Embedding matrix (L, D)
        valid_emb = embeddings_np[i, residue_idx, :]
        
        # Slicing the 2D Attention matrix (L, L) along both axes
        valid_attn = pooled_attention_np[i][np.ix_(residue_idx, residue_idx)]
        
        unpadded_embeddings.append(valid_emb)
        unpadded_attentions.append(valid_attn)
        
    # Return both for your downstream pipeline
    return unpadded_embeddings, unpadded_attentions

model = E1ForMaskedLM.from_pretrained("Profluent-Bio/E1-150m", attn_implementation="eager").to("cuda:0")
model.eval()

sequences = ["AAA?C", "CASDF,MFCC?S,F", "DFSMF"]
batch_preparer = E1BatchPreparer()
batch = batch_preparer.get_batch_kwargs(sequences, device="cuda:0")
seq_original_lens = [len(seq.replace(',', '')) for seq in sequences]
print("Token input IDs:")
for i, ids in enumerate(batch["input_ids"]):    
    seq = sequences[i]
    print(f"Sequence {i}: {list(seq)}")
    print(f"IDs: {[int(x) for x in list(ids)]}")

    

with torch.autocast("cuda", dtype=torch.bfloat16, enabled=True):
    outputs = model(
        input_ids=batch["input_ids"],
        within_seq_position_ids=batch["within_seq_position_ids"],
        global_position_ids=batch["global_position_ids"],
        sequence_ids=batch["sequence_ids"],
        past_key_values=None,
        use_cache=False,
        output_attentions=True,
        output_hidden_states=False,
    )

    logits: torch.Tensor = outputs.logits  # (B, L, V)
    embeddings: torch.Tensor = outputs.embeddings  # (B, L, E)
    embeddings = embeddings.detach().cpu().float().numpy()
    pooled_attention = None
            
    for layer_attn in outputs.attentions:
        # layer_attn shape is (Batch, Heads, SeqLen, SeqLen)
        # 1. Cast to FP32 for precision
        # 2. Max pool across heads (dim=1) for this specific layer
        layer_max = layer_attn.to(torch.float32).amax(dim=1) 
        
        # 3. Update the global running max across layers
        if pooled_attention is None:
            pooled_attention = layer_max
        else:
            pooled_attention = torch.max(pooled_attention, layer_max)
            
    # Move ONLY the final aggregated 3D tensor to the CPU
    pooled_attention_np = pooled_attention.detach().cpu().numpy()
    print("Embeddings shape:", embeddings.shape)
    print(embeddings)
    print("Pooled attention shape:", pooled_attention_np.shape)
    print(pooled_attention_np)

    unpadded_embeddings, unpadded_attentions = unpad(seq_original_lens, batch, embeddings, pooled_attention_np)

    for seq_i in range(len(sequences)):
        print("Sequence", seq_i, ":", sequences[seq_i], ", Unpadded len:", seq_original_lens[seq_i])
        print("Unpadded embedding:", unpadded_embeddings[seq_i].shape)
        print("\t", unpadded_embeddings[seq_i])
        print("Unpadded attention:", unpadded_attentions[seq_i].shape)
        print("\t", unpadded_attentions[seq_i])
    
