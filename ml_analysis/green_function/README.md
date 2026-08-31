# Green-function analysis

This directory contains the reproducible analysis code and lightweight
configuration for the green-function work. Large inputs and generated
artifacts are intentionally kept outside the Git repository.

## Repository layout

```text
green_function/
├── README.md
├── configs/
├── data/       # documentation only; actual inputs are user supplied
└── scripts/
```

## Paths

- Input data must be supplied explicitly with `--data-dir`.
- Generated data, figures, and logs are written below
  `/global/cfs/cdirs/m1867/zhan391/paper_work/green_function`.
- The external output root can be overridden with `--output-root` when needed.

Each analysis creates separate `outputs/`, `figures/`, and `logs/`
subdirectories beneath its output root. See `configs/paths.toml` for the
documented default.

For example:

```bash
python green_function/scripts/abrupt4CO2/decode_data_abrupt4CO2_FOM.py \
  --data-dir /path/to/input/data
```

