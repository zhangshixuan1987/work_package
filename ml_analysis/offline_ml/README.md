# Offline ML analysis

This directory contains reproducible offline-ML analysis code. Large input
datasets and generated artifacts are kept outside Git.

## Repository layout

```text
offline_ml/
├── README.md
├── configs/
├── data/       # documentation only; inputs are user supplied
└── scripts/    # Python, notebooks, NCL, and launch scripts
```

The default external workspace is
`/global/cfs/cdirs/m1867/zhan391/paper_work/offline_ml`, organized into
`data/`, `outputs/`, `figures/`, and `logs/`. New scripts should accept input
paths from users and write generated artifacts beneath this external root.

