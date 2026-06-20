import cv2
import os
import glob
import sys
import imageio.v3 as iio
import numpy as np


def generate_pca_video_cv2(directory, output_filename="latent_evolution.mp4", fps=2.5):
    """
    Compila todos os scatter plots do PCA salvos em um vídeo MP4.
    fps: frames por segundo (2.5 fps = 400ms por frame).
    """
    if directory is None:
        return

    print("\nStitching PCA frames into Video...")

    # Encontra todas as imagens e ordena alfabeticamente/cronologicamente
    search_pattern = os.path.join(directory, "pca_epoch_*.png")
    image_paths = sorted(glob.glob(search_pattern))

    if not image_paths:
        print("No PCA images found to make a video.")
        return

    # Lê o primeiro frame para extrair as dimensões necessárias para o vídeo
    first_frame = cv2.imread(image_paths[0])
    height, width, layers = first_frame.shape
    size = (width, height)

    output_path = os.path.join(directory, output_filename)

    # Define o codec mp4v (padrão, sem dependências externas complexas) e cria o VideoWriter
    fourcc = cv2.VideoWriter_fourcc(*"mp4v")
    video = cv2.VideoWriter(output_path, fourcc, fps, size)

    # Adiciona cada imagem como um frame no vídeo
    for img_path in image_paths:
        frame = cv2.imread(img_path)
        video.write(frame)

    # Libera o arquivo e finaliza a gravação
    video.release()
    print(f"Video successfully saved at: {output_path}")


def generate_pca_video(directory, output_filename="latent_evolution.mp4", fps=2.5):
    """
    Compila os frames do PCA em um vídeo MP4 compatível com WhatsApp (codec H.264).
    """
    if directory is None:
        return

    print("\nStitching PCA frames into Video (WhatsApp compatible)...")

    search_pattern = os.path.join(directory, "pca_epoch_*.png")
    image_paths = sorted(glob.glob(search_pattern))

    if not image_paths:
        print("No PCA images found to make a video.")
        return

    print("image_paths", image_paths)
    print("directory", directory)
    print("output_filename", output_filename)
    print("fps", fps)

    output_path = os.path.join(directory, output_filename)

    # Lê todas as imagens da pasta para a memória
    frames = [iio.imread(img) for img in image_paths]

    # Empilha as imagens em um único array 4D (Frames, Altura, Largura, Canais de Cor)
    video_array = np.array(frames)

    # Grava o vídeo já padronizado no formato H.264
    iio.imwrite(
        output_path,
        video_array,
        fps=fps,
        extension=".mp4",
        macro_block_size=2,  # Critical for preserving plot sharpness
        ffmpeg_params=[
            "-crf",
            "17",  # Constant Rate Factor (Quality)
            "-preset",
            "slower",  # Encoding effort
            "-pix_fmt",
            "yuv420p",  # Broadest compatibility
        ],
    )

    assert os.path.exists(output_path)

    print(f"Video successfully saved at: {output_path}")


if __name__ == "__main__":
    dirname = sys.argv[1]
    generate_pca_video(dirname)
