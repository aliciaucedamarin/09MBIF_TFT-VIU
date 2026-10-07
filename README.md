# 09MBIF_TFT-VIU
Final Master's Degree Project at VIU

## Project structure

```text
09MBIF_TFT/
│
├── .gitignore
├── 09MBIF_TFT.Rproj
├── LICENSE
├── README.md
│
├── code/
│   └── # Analysis scripts
│
├── data/
│   └── original/
│       ├── all_24_content_clean_annot_nodup.bed.xlsx
│       └── Panel-Regions.bed
│
├── metadata/
│   ├── patients_metadata.txt
│   └── original/
│       └── INFO-PATIENTS_Methylseq.xlsx
│
└── results/
    ├── figures/
    │   └── *.pdf / *.png
    │
    └── tables/
        └── *.txt / *.csv / *.tsv
```

## Description

This project contains the scripts and files used for data processing and analysis, with the generated results organized into:

* `results/figures/`: figures generated during the analysis.
* `results/tables/`: tables and summary files generated during the analysis.

## Data

The `data/` and `metadata/` directories contain files used as inputs for the analyses.

Some files contain **confidential or sensitive information**. Therefore, certain files are not included in version control and are not synchronized with GitHub. These files are kept only in the local working environment.

The `.gitignore` file is used to prevent confidential or sensitive files from being accidentally uploaded to the repository.

## Code

Analysis scripts are stored in:

```text
code/
```

At this stage, the repository mainly contains the analysis code and project structure. The scripts generate the corresponding tables and figures in the `results/` directory.

## Results

Analysis results are organized as follows:

```text
results/
├── figures/
└── tables/
```

Figures are mainly generated in PDF or PNG format, while tables and summary files are generated in TXT, CSV, or TSV format.

## Reproducibility

Reproducing the analyses requires access to the corresponding input files. Due to the confidential nature of some of the data, certain input files are not included in the repository and are therefore not publicly available.

The project is developed using **R** and an RStudio project file (`09MBIF_TFT.Rproj`).

## License

This project is distributed under the terms specified in the `LICENSE` file.