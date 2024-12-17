
from transformers import AutoTokenizer, AutoModelForSeq2SeqLM
from transformers import EsmTokenizer, EsmModel

model_names = [
    'facebook/esm2_t6_8M_UR50D',
    'facebook/esm2_t12_35M_UR50D',
    'facebook/esm2_t30_150M_UR50D',
    'facebook/esm2_t33_650M_UR50D',
    'facebook/esm2_t36_3B_UR50D',
    'ElnaggarLab/ankh-base',
    'ElnaggarLab/ankh-large',
    'ElnaggarLab/ankh2-ext1',
    'ElnaggarLab/ankh2-ext2'
]

for mn in model_names:
    if 'esm' in mn:
        model = EsmModel.from_pretrained(mn,).to(device='cpu')
    else:
        model = AutoModelForSeq2SeqLM.from_pretrained(mn, trust_remote_code=True).to(device='cpu')

    print(mn+':', model.num_parameters() / 1000000, 'Million parameters')