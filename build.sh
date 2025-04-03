mkdir -p singularity_images
singularity build --fakeroot singularity_images/env2.sif conda_envs/env2.def > singularity_images/env2.stdout
singularity build --fakeroot singularity_images/basic_env.sif conda_envs/basic_env.def > singularity_images/basic_env.stdout