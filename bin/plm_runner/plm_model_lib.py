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
            attn_implementation="eager"
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

    def load_model_and_tokenizer(self, model_name):
        from E1.batch_preparer import E1BatchPreparer
        from E1.modeling import E1ForMaskedLM

        model = E1ForMaskedLM.from_pretrained(model_name)
        model.to(device=self.device)
        model.eval()

        batch_preparer = E1BatchPreparer()
        return model, batch_preparer
    
    def extract(self, seqs: List[str]):
        bf16_supported = torch.cuda.is_available() and torch.cuda.is_bf16_supported()
        batch = self.batch_preparer.get_batch_kwargs(seqs, device="cuda:0")

        dtype = torch.bfloat16 if bf16_supported else torch.float32
        with torch.autocast("cuda", dtype=dtype, enabled=True):
            outputs = self.model(
                input_ids=batch["input_ids"],
                within_seq_position_ids=batch["within_seq_position_ids"],
                global_position_ids=batch["global_position_ids"],
                sequence_ids=batch["sequence_ids"],
                past_key_values=None,
                use_cache=False,
                output_attentions=False,
                output_hidden_states=False,
            )
        
        logits: torch.Tensor = outputs.logits  # (B, L, V)
        embeddings: torch.Tensor = outputs.embeddings  # (B, L, E)

        print(logits)
        print(embeddings)

        embeddings = [emb.cpu().numpy()
            for emb in embeddings]

        last_emb = embeddings[-1]
        print(f"Embeddings desc: shape={embeddings.shape}, dtype={embeddings.dtype}")
        print(f"Last emb desc: shape={last_emb.shape}, dtype={last_emb.dtype}")

        return embeddings, None

    
def plm_master_loader(model_name, cache_path):
    plm_type = define_plm_class(model_name)
    if plm_type == "ANKH":
        return ANKHModel(model_name, cache_path)
    
    return None