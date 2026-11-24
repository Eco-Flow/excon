# Running excon without internet access

Many HPC clusters give compute nodes no internet access. excon runs there, but a few things have to be downloaded beforehand, on a machine that does have access, such as a login node.

## What excon downloads

| What | When | How to avoid it |
|---|---|---|
| Nextflow itself and its plugins (`nf-schema`, `nf-amazon`) | When Nextflow starts | Run Nextflow once with internet access, then set `NXF_OFFLINE=true` |
| Container images | Before each step first runs | Pre-pull them with `utils/pull_singularity_images.sh` |
| Genomes given as NCBI accessions (2-column input rows) | `NCBIGENOMEDOWNLOAD` | Download them yourself and use 3-column rows (`name,genome.fna,annotation.gff`) |
| The eggNOG database (~45 GB) | `EGGNOG_DOWNLOAD`, with `--run_eggnog` | Download it once and pass `--eggnog_data_dir` |
| BUSCO lineage datasets | `BUSCO_BUSCO`, with `--stats` | Download the lineage and pass `--busco_lineage <name>` and `--busco_lineages_path` |

Everything else runs offline.

## The two set-ups

### Nextflow on the login node, steps on compute nodes

This is the simplest set-up. Nextflow's own process stays on the login node (or another node with internet access), and submits each step to the scheduler with an executor configuration (e.g. `process.executor = 'slurm'` in a `-c` config file). With Singularity or Apptainer, Nextflow downloads each image itself, from that node, before the steps that use it are submitted. So compute nodes never need internet access for containers. Only the downloads in the table above run on compute nodes, so supply those inputs locally.

Many clusters don't allow long-running processes on login nodes. Ask your administrators, or use the next set-up.

### Everything offline

1. **On a machine with internet access**, from a clone of excon, download the container images into a shared directory:

   ```bash
   utils/pull_singularity_images.sh /shared/path/excon_images
   ```

   It lists the image of every step in the pipeline with `nextflow inspect`, whether or not your run uses that step. It then pulls each image under the file name Nextflow looks for. Use `--dry-run` to see the list first, and `--engine apptainer` if you have Apptainer and Singularity together.

   The script names no image or version itself: the images and the registry (`quay.io`, from the profile's `registry` setting) come from the pipeline's own config. So it stays correct as excon updates its tools, but **re-run it after you update excon** (`git pull` or a new release): steps whose tool version changed need new images. Re-running only fetches the images that are missing, and old ones can be deleted. Each pull request's CI checks that the script's file names still match the ones Nextflow uses.

   Running Nextflow here also downloads Nextflow's own files and plugins into `~/.nextflow`. That directory needs to be readable from the offline nodes (or copy it there, or set `NXF_HOME` to a shared path for both runs).

2. **Download the other inputs** in the table above that your run needs: genomes, the eggNOG database (`bin/download_eggnog_db.sh`), and BUSCO lineages.

3. **On the offline nodes**, point Nextflow at the images and tell it not to look online:

   ```bash
   export NXF_OFFLINE=true
   export NXF_SINGULARITY_LIBRARYDIR=/shared/path/excon_images   # NXF_APPTAINER_LIBRARYDIR with -profile apptainer
   nextflow run main.nf -profile singularity --input local_genomes.csv ...
   ```

   Use `NXF_SINGULARITY_LIBRARYDIR` rather than `NXF_SINGULARITY_CACHEDIR`. The `singularity` and `apptainer` profiles set their own cache directory (`singularity/cachedir` in the launch directory), and that setting overrides `NXF_SINGULARITY_CACHEDIR`. Nextflow checks the library directory before the cache. Alternatively, copy the images into `singularity/cachedir` in the directory you launch from.

## Starting with internet access and resuming offline

If it is allowed on your cluster, you can also start a run where there is internet access, stop it once the images and downloads are done, and continue on the offline nodes with `-resume`. This only works if you resume from the same launch directory, with the same `work/` directory and options. The steps that already finished, including the downloads, then come from the cache. Pre-pulling the images as above is more reliable, because a resumed run still needs the images for every step it hasn't run yet.

## Docker

With `-profile docker`, each node's Docker daemon pulls the images, not Nextflow. On offline nodes, load the images into each node's Docker beforehand (`docker save` / `docker load`). Use the image names from `nextflow inspect main.nf -profile docker --outdir tmp`, prefixed with `quay.io/` where they name no registry.
