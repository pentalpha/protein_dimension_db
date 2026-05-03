include {
    train_interpro_autoencoder_256 ;
    train_interpro_autoencoder_512 ;
    train_interpro_autoencoder_640 ;
    train_interpro_autoencoder_896 ;
    train_interpro_autoencoder_1280
} from './modules/interpro_processes.nf'

params.mode = "release"

params.max_protein_len = 1800

workflow {
    join_interpro_consults_path = params.release_dir + "/raw_data/interpro_consults/concatenated.tsv"
    interpro_vocab_ia_json_path = params.release_dir + "/raw_data/interpro_vocab/interpro_vocab_ia.json"

    train_interpro_autoencoder_256(
        join_interpro_consults_path,
        interpro_vocab_ia_json_path,
    )

    train_interpro_autoencoder_512(
        join_interpro_consults_path,
        interpro_vocab_ia_json_path,
    )

    train_interpro_autoencoder_640(
        join_interpro_consults_path,
        interpro_vocab_ia_json_path,
    )

    train_interpro_autoencoder_896(
        join_interpro_consults_path,
        interpro_vocab_ia_json_path,
    )

    train_interpro_autoencoder_1280(
        join_interpro_consults_path,
        interpro_vocab_ia_json_path,
    )
}
