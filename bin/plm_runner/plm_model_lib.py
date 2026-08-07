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
            model_name, output_attentions=False, token=self.token
        )
        model.to(device=self.device)
        model.eval()
        return model, tokenizer

    def extract(self, seqs: List[str]):
        shift_left = 0
        shift_right = -1
        seqs = [list(seq) for seq in seqs]
        seq_original_lens = [len(seq) for seq in seqs]
        with torch.no_grad():
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
            embeddings = embeddings.last_hidden_state.cpu().numpy()
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

        return embeddings

    
def plm_master_loader(model_name, cache_path):
    plm_type = define_plm_class(model_name)
    if plm_type == "ANKH":
        return ANKHModel(model_name, cache_path)
    
    return None