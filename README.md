# R-australiae-Genomics

Scripts and workflow notes for "Population and seascape genomics of the Critically Endangered bottlenose wedgefish (*Rhynchobatus australiae*) from the Indo-West Pacific region" (Groeneveld et al. in prep).

## Pipeline

| Stage | Contents | Where it runs | Document |
|---|---|---|---|
| 1 | Raw data QC, demultiplexing, de novo assembly (Stacks), filtering | HPC, with some checks in R | [`docs/01_assembly_and_filtering.md`](docs/01_assembly_and_filtering.md), [`scripts/01_mMn_optimisation.R`](scripts/01_mMn_optimisation.R) |
| 2 | Environmental data, variable selection, partial RDA, outlier detection, four-dataset preparation | R, with ordistep and partial RDA on HPC | [`docs/02_seascape_outliers.md`](docs/02_seascape_outliers.md), [`scripts/02_seascape_outliers.R`](scripts/02_seascape_outliers.R), [`scripts/set1_RDA_ordistep11.R`](scripts/set1_RDA_ordistep11.R), [`scripts/set1_partialRDA.R`](scripts/set1_partialRDA.R) |
| 3 | Population genomics (diversity, FST, DAPC, sNMF, IBD/IBE) | R | [`docs/03_population_genomics.md`](docs/03_population_genomics.md), [`scripts/03_population_genomics.R`](scripts/03_population_genomics.R) |

## Repository layout

```
.
├── README.md
├── docs/
│   ├── 01_assembly_and_filtering.md
│   ├── 02_seascape_outliers.md
│   └── 03_population_genomics.md
└── scripts/
    ├── 01_mMn_optimisation.R
    ├── 02_seascape_outliers.R
    ├── 03_population_genomics.R
    ├── set1_RDA_ordistep11.R
    └── set1_partialRDA.R
```

Data files (VCFs, popmaps, coordinates, environmental data) are not included. R version 4.4.1 was used locally (packages are loaded at the top of each script); 4.3.2 was used for the cluster R jobs (ordistep and partial RDA).

## Running the pipeline

Run stage 1, then stage 2 (with the two cluster jobs), then stage 3. Each document lists the input files its scripts read under "Input". The R scripts are written to be run section by section, in order, from the top of the file.

## Data availability

Raw sequence data: [<figshare>]. Sample metadata and locality information are not included in this repository.

## Contact

ResearchGate: [https://www.researchgate.net/profile/Mia-Groeneveld](https://www.researchgate.net/profile/Mia-Groeneveld).
