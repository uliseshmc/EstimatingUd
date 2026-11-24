# EstimatingUd

A research pipeline for estimating Ud, the deleterious mutation rate, from
Human/Chimpanzee/Orangutan whole-genome alignments (Ensembl 114, downloaded with
`ensembl_tui`) using `cogent3` substitution models.

The analysis code is an installable package under `src/estimating_ud/` with a Click
command line, `ud`. Visualisation and summary statistics live in `notebooks/`, SLURM
wrappers in `scripts/`, tables and figures in `results/`. The alignments themselves
(terabytes) stay outside the repository; every command finds them through the
`DATA_HUMCHIMPORANGOR114` environment variable, or through a `.env` file in the repo
root when the variable is not exported.

```
pyproject.toml            package metadata and dependency pins
requirements-udenv.txt    exact pip freeze of the project env (install with --no-deps)
src/estimating_ud/        apps.py, models.py, dinuc_models.py, regions.py, config.py, cli.py, stages/
scripts/slurm/            one array job per stage, calling `ud`
scripts/download/         one array job per region, calling `eti`
scripts/submit.sh         creates logs/ then sbatch
notebooks/                01 ... 05 in run order; exploratory/ for the alignment-quality work
results/tables, figures   outputs of stage 5 and of the notebooks
tests/                    pytest suite on synthetic fixtures (no real data needed)
docs/                     the pipeline documents, read them before changing a stage
archive/                  earlier generations of the pipeline, kept for reference only
```

## Install

One micromamba environment, `Udenv`, runs everything, downloads included:

```
micromamba create -f environment.yml
micromamba activate Udenv
pip install --no-deps -r requirements-udenv.txt   # exact versions; see pyproject.toml for why --no-deps
pip install --no-deps -e .
python -m ipykernel install --user --name Udenv
cp .env.example .env   # then edit: DATA_HUMCHIMPORANGOR114=/path/to/HumChimpOranGor_10_114
pytest
```

## Run

```
ud --help
ud align-cds      --chrm 22 --submodel singlent      # stage 1 (cds only)
ud genic-coords                                      # utr5/utr3/gene-span tables, before their download
ud filter-noncds  --chrm 22 --submodel singlent      # stage 3 (there is no stage 2)
ud human-seq      --chrm 22                          # Human sequence per region
ud fit            --chrm 22 --submodel singlent      # stage 4
ud constraint                                        # stage 5, writes results/tables/
ud quality --help                                    # alignment-quality branch
```

Repeat stages 1 to 4 with `--submodel trinuc`. On the cluster, submit the wrappers from the
repo root: `scripts/submit.sh scripts/slurm/04_fit.sh`. The full order, from download to
Ud, is in `docs/pipeline_downloaddata.md`, `docs/pipeline_filtering.md` and
`docs/pipeline_constraint_estimates.md`. The alignment-quality side branch (`ud quality`) is
described in `docs/pipeline_qualityaligned.md`.
