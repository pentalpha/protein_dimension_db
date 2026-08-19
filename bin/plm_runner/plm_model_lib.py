from typing import List, Tuple

import numpy as np
import torch
from transformers import (
    T5EncoderModel,
    #T5ForConditionalGeneration,
    AutoTokenizer,
    #TFT5EncoderModel,
    #TFT5ForConditionalGeneration,
    T5Tokenizer,
)

from plm_runner.plm_model import PLMModel, define_plm_class

BF16_SUPPORT = torch.cuda.is_available() and torch.cuda.is_bf16_supported()

class ANKHModel(PLMModel):
    def __init__(self, model_name: str, cache_path):
        super().__init__(model_name, cache_path)
        self.model, self.tokenizer = self.load_model_and_tokenizer(model_name)

    def load_model_and_tokenizer(self, model_name: str) -> Tuple[T5EncoderModel, AutoTokenizer]:
        """Downloads and returns the base model and its tokenizer
        Returns:
            Tuple[Union[T5EncoderModel, T5ForConditionalGeneration],
            AutoTokenizer]: Returns T5 Model and its tokenizer.
        """
        if "ankh3" in model_name:
            tokenizer = T5Tokenizer.from_pretrained(model_name, token=self.token)
        else:
            tokenizer = AutoTokenizer.from_pretrained(model_name, token=self.token)
        model = T5EncoderModel.from_pretrained(
            model_name, output_attentions=False, token=self.token,
            attn_implementation="eager", 
            #use_safetensors=True
        )
        model.to(device=self.device)
        model.eval()
        self.special_token_ids = set(tokenizer.all_special_ids)
        #model.compile()
        return model, tokenizer

    def extract_vanilla(self, seqs: List[str]):
        seqs = [list(seq) for seq in seqs]
        seq_original_lens = [len(seq) for seq in seqs]
        dtype = torch.bfloat16 if BF16_SUPPORT else torch.float16
        with torch.no_grad(), torch.autocast("cuda", dtype=dtype, enabled=True):
            tokenized = self.tokenizer.batch_encode_plus(
                seqs,
                add_special_tokens=True,
                padding=True,
                is_split_into_words=True,
                return_tensors="pt",
            )
            '''for inputs_vec in tokenized["input_ids"]:
                print(inputs_vec.shape)
                print(inputs_vec)'''
            input_ids = tokenized["input_ids"].to(self.device)
            attention_mask= tokenized['attention_mask'].to(self.device)
            embeddings = self.model(input_ids=input_ids, attention_mask=attention_mask)
            embeddings = embeddings.last_hidden_state.cpu().to(torch.float32).numpy()
            # Dynamically slice out the exact amino acids
            unpadded_embeddings = []
            for i in range(len(seqs)):
                orig_len = seq_original_lens[i]
                # Start at 0 (no prepended start token) 
                # End at orig_len (drops the EOS token and all padding)
                valid_emb = embeddings[i][0 : orig_len]
                unpadded_embeddings.append(valid_emb)
            embeddings = unpadded_embeddings

            '''for i in range(len(seqs)):
                len1 = seq_original_lens[i]
                len2 = embeddings[i].shape[0]
                print(f"Seq {i}: {''.join(seqs[i])}")
                print(f"\t Original len: {len1}; New len: {len2}")
                padding_size = len2 - len1'''
            return embeddings
    
    def extract(self, seqs: List[str]):
        seq_words = [list(seq) for seq in seqs]
        seq_original_lens = [len(seq) for seq in seqs]
        
        # Check for bfloat16 support dynamically (assumes BF16_SUPPORT is defined in your env)
        bf16_supported = torch.cuda.is_available() and torch.cuda.is_bf16_supported()
        dtype = torch.bfloat16 if bf16_supported else torch.float16
        
        with torch.no_grad(), torch.autocast("cuda", dtype=dtype, enabled=True):
            tokenized = self.tokenizer.batch_encode_plus(
                seq_words,
                add_special_tokens=True,
                padding=True,
                is_split_into_words=True,
                return_tensors="pt",
            )
            
            input_ids = tokenized["input_ids"].to(self.device)
            attention_mask = tokenized['attention_mask'].to(self.device)
            
            # Explicitly request attentions for PoolPARTI
            outputs = self.model(
                input_ids=input_ids, 
                attention_mask=attention_mask,
                output_attentions=True 
            )
            
            # Embeddings -> FP32 -> NumPy
            embeddings_np = outputs.last_hidden_state.cpu().to(torch.float32).numpy()
            
            # Calculate the running maximum iteratively to prevent VRAM explosion
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
            pooled_attention_np = pooled_attention.cpu().numpy()
            
            unpadded_embeddings = []
            unpadded_attentions = []
            
            # Move input_ids to CPU list to check against special tokens
            input_ids_list = input_ids.cpu().tolist()
            
            for i in range(len(seqs)):
                # Dynamically construct the biological residue mask!
                # This ignores [CLS], [EOS], and [PAD] automatically.
                residue_idx = [
                    idx for idx, token_id in enumerate(input_ids_list[i]) 
                    if token_id not in self.special_token_ids
                ]
                
                # Sanity check to ensure our mask perfectly matches the input sequence length
                assert len(residue_idx) == seq_original_lens[i], f"Mask length mismatch at seq {i}"
                
                # Slicing the 2D Embedding matrix (L, D)
                valid_emb = embeddings_np[i, residue_idx, :]
                
                # Slicing the 2D Attention matrix (L, L) along both axes
                valid_attn = pooled_attention_np[i][np.ix_(residue_idx, residue_idx)]
                
                unpadded_embeddings.append(valid_emb)
                unpadded_attentions.append(valid_attn)
                
            # Return both for your downstream pipeline
            return unpadded_embeddings, unpadded_attentions

class ESMModel(PLMModel):
    def __init__(self, model_name: str, cache_path):
        super().__init__(model_name, cache_path)
        pass

class DPLMModel(PLMModel):
    def __init__(self, model_name: str, cache_path):
        super().__init__(model_name, cache_path)
        pass

class ProfluentE1Model(PLMModel):
    def __init__(self, model_name: str, cache_path):
        super().__init__(model_name, cache_path)
        
        self.model, self.batch_preparer = self.load_model_and_tokenizer(model_name)
        self.single_special_tokens = set([0])
        self.tuple_special_tokens = set([(1,6), (7,2)])

    def load_model_and_tokenizer(self, model_name):
        from E1.batch_preparer import E1BatchPreparer
        from E1.modeling import E1ForMaskedLM

        #device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        model = E1ForMaskedLM.from_pretrained("Profluent-Bio/E1-150m").to(self.device)
        batch_preparer = E1BatchPreparer()
        model.eval()
        return model, batch_preparer
    
    def unpad(self, seq_original_lens, batch, embeddings_np, pooled_attention_np):

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
                for special_tuple in self.tuple_special_tokens:
                    t_len = len(special_tuple)
                    # Look ahead to see if the next tokens match the tuple exactly
                    if idx + t_len <= len(seq_ids) and tuple(seq_ids[idx:idx+t_len]) == special_tuple:
                        idx += t_len  # Skip over the entire boundary sequence
                        matched_tuple = True
                        break
                        
                if matched_tuple:
                    continue
                    
                # 2. Check for single-token special characters (e.g., 0 for PAD)
                if seq_ids[idx] in self.single_special_tokens:
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

    
    def extract(self, seqs: List[str]):
        bf16_supported = torch.cuda.is_available() and torch.cuda.is_bf16_supported()
        dtype = torch.bfloat16 if bf16_supported else torch.float32
        batch = self.batch_preparer.get_batch_kwargs(seqs, device=self.device)
        
        with torch.autocast(str(self.device), dtype=dtype, enabled=True):
            outputs = self.model(
                input_ids=batch["input_ids"],
                within_seq_position_ids=batch["within_seq_position_ids"],
                global_position_ids=batch["global_position_ids"],
                sequence_ids=batch["sequence_ids"],
                past_key_values=None,
                use_cache=False,
                output_attentions=True,
                output_hidden_states=False,
            )

        embeddings: torch.Tensor = outputs.embeddings  # (B, L, E)
        embeddings = embeddings.detach().cpu().float().numpy()
        pooled_attention = None
    
        # 1. Convert the tuple to a list so we can actively pop elements out of it
        attentions_list = list(outputs.attentions)
        
        # 2. Delete the original tuple reference to allow Python's Garbage Collector to work
        outputs.attentions = None 
        
        while len(attentions_list) > 0:
            # pop(0) extracts the layer and removes it from the list
            layer_attn = attentions_list.pop(0)
            
            # 3. Offload the tensor to CPU memory immediately
            layer_attn_cpu = layer_attn.detach().cpu()
            
            # 4. Explicitly delete the GPU tensor reference to free up VRAM
            del layer_attn 
            
            # 5. Now do the heavy FP32 casting and pooling safely on the CPU
            layer_max = layer_attn_cpu.to(torch.float32).amax(dim=1) 
            
            # 6. Update the global running max (also on the CPU)
            if pooled_attention is None:
                pooled_attention = layer_max
            else:
                pooled_attention = torch.max(pooled_attention, layer_max)
                
        # pooled_attention is already a detached CPU tensor at this point
        pooled_attention_np = pooled_attention.numpy()

        seq_original_lens = [len(seq.replace(',', '')) for seq in seqs]
        unpadded_embeddings, unpadded_attentions = self.unpad(seq_original_lens, batch, embeddings, pooled_attention_np)

        return unpadded_embeddings, unpadded_attentions

    
def plm_master_loader(model_name, cache_path):
    plm_type = define_plm_class(model_name)
    if plm_type == "ANKH":
        return ANKHModel(model_name, cache_path)
    elif plm_type == "PROFLUENT":
        return ProfluentE1Model(model_name, cache_path)
    
    return None