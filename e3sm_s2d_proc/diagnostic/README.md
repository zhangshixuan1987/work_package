# RTENSO2026 case-study diagnostics

Diagnostics for the RTENSO2026 2015 hindcasts (`../hindcast_run/RTENSO2026/`,
post-processed by `../postprocess/RTENSO2026/`), treated as a case study of the
2015/16 El Nino. The notebooks follow the layout and conventions of the ESP-Lab
suite (`~/code/ESP-Lab/jupyter`): a run-control block at the top, cache-first
products under `S2D_DIAG_ROOT`, and figures named with `figure_filename`. ESP-Lab
measures skill over many start years; here one event is followed through 4 starts.

| Experiment | Starts (2015) | Members | Differs by |
| --- | --- | --- | --- |
| `Reanalysis-CLIM` (`RTENSO2026_clim`) | Feb, May, Aug, Nov | 10 | land/sea-ice/river ICs from the BruteForce year closest to the same-date climatology |
| `Reanalysis-OBS` (`RTENSO2026_obs`) | Feb, May, Aug, Nov | 10 | land/sea-ice/river ICs from the BruteForce year closest to the ERA5 anomaly |
| `RTENSO2026-obs-nsidc` (`RTENSO2026_obs_ice_nsidc`) | placeholder | 10 | `Reanalysis-OBS` with NSIDCv6 sea ice; enabled by setting its `init_months` once post-processed |
| `Reanalysis-FOSIRL` (BruteForce) | Feb, May, Aug, Nov | 10 | the actual 2015 restarts |

All share the ERA5 atmosphere and CFSR ocean ICs. Cases whose `init_months` is empty
are placeholders and are dropped from the notebooks.

## Drift correction

Model anomalies remove the Reanalysis (BruteForce) 1981-2010 lead-time climatology,
the same IC system:

- May and Nov starts use the Reanalysis climatology for that start month and lead.
- Feb and Aug starts have no Reanalysis climatology (only 1980 and 1985 exist). The
  default `CLIM_METHOD = "lead_drift"` uses the observed 1981-2010 climatology of each
  target month plus the Reanalysis drift at the same lead, averaged over the May and
  Nov starts. `"target_month"` instead borrows the May/Nov climatology for the same
  target month at the nearest lead; it is biased at short leads because the
  equatorial cold tongue drifts about 1 C in the first months (Nino3.4 error at the
  Aug start, lead 0: +0.82 C with `target_month`, -0.10 C with `lead_drift`).

Observed anomalies are relative to the 1981-2010 monthly climatology: HadISST2 SST
(land-masked, verified 60S-60N), ERA5 `TREFHT`, `PSL`, `U850`, and GPCP v2.3 `PRECT`,
all conservatively regridded to the model 1x1 grid.

## Files

| File | Purpose |
| --- | --- |
| `rtenso_case.py` | Shared helpers: case registry, hindcast loading, lead-time climatologies, observation anomalies, pattern ACC/RMSE, regions. |
| `build_caches.py` | Build the Reanalysis lead-time climatologies and observation caches ahead of the notebooks (`python build_caches.py [FIELD ...] [--mode rebuild]`). |
| `1_enso_plume.ipynb` | Nino3.4/3/4/1+2 forecast plumes per start against HadISST2, peak-season summary, ensemble-mean error and spread (after 3a). |
| `2_anomaly_maps_acc.ipynb` | DJF 2015/16 anomaly maps from every start, and spatial pattern ACC and RMSE by lead season and region, per member and ensemble mean (after 1a/1b). Set `FIELD` to `TS`, `TREFHT`, `PRECT` or `PSL`. |
| `3_ic_sensitivity.ipynb` | `Reanalysis-CLIM` vs `Reanalysis-OBS` vs `Reanalysis-FOSIRL` differences: land-mean differences with member-based intervals, fraction of land with significant differences by lead, stippled difference maps, and hemispheric sea-ice area against HadISST2. |
| `4_teleconnections.ipynb` | DJF 2015/16 teleconnections of `PRECT`/`TREFHT`/`PSL` per degC of Nino3.4 against the observed event and the 1982-2010 observed regression, an inter-member regression, and equatorial Hovmollers of SST and U850 (after 3b). |

## Running

Use the `e3sm_analysis` kernel (ESP-Lab installed in editable mode) and start Jupyter in
this directory so `rtenso_case` imports. Build the caches once (about 4 minutes per
field and start month; the U850 reads are the slowest):

```bash
conda activate e3sm_analysis
python build_caches.py
```

Each notebook's run-control block holds `derivation_mode` ("auto" reuses valid caches,
"rebuild" recomputes, "require" fails on a missing cache) and `CLIM_METHOD`.

Outputs:

- Caches: `$S2D_DIAG_ROOT/Reanalysis/leadtime_clim/`, `$S2D_DIAG_ROOT/observations/case_study/`
  and `$S2D_DIAG_ROOT/RTENSO2026/case_study/` (`S2D_DIAG_ROOT` defaults to
  `/global/cfs/cdirs/e3sm/S2S2D/s2d_diag`).
- Figures: `/global/cfs/cdirs/e3sm/www/zhan391/e3sm-s2d_rtenso2026_case/`, kept apart from
  the ESP-Lab gallery.

Land and sea-ice diagnostics in `3_ic_sensitivity.ipynb` need the RTENSO2026 `lnd_mon1`
and `ice_mon1` post-processing; cases not yet processed are skipped.
