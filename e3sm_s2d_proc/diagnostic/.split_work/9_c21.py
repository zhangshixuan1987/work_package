DJF_DAYS = np.array([31.0, 31.0, 28.0])          # day weights of Dec, Jan, Feb (no-leap model calendar)


def djf_mean(series, w):
    # day-weighted DJF mean of a monthly series stamped on the 15th; NaN unless all three months exist
    times = [np.datetime64(f"{w + dy}-{m:02d}-15") for dy, m in DJF]
    if not all(t in series["time"].values for t in times):
        return None
    v = series.sel(time=times)
    return (v * xr.DataArray(DJF_DAYS, dims="time", coords={"time": v["time"]})).sum("time") / DJF_DAYS.sum()


def djf_leads(start_month, k):
    # leads of Dec, Jan, Feb of the k-th winter after a start month
    dec = (12 - start_month) % 12 + 12 * k
    return [dec, dec + 1, dec + 2]


def ref_djf_anomalies(m, k):
    # Reanalysis DJF anomalies (year, member) of start month m in its k-th winter, about the mean over years and members
    leads = djf_leads(m, k)
    if leads[-1] >= REF_ALL[m].sizes["lead"]:
        return None
    u = sc.lat_mean(REF_ALL[m]["U"].sel(plev=REF_LEVEL), sc.U60).isel(lead=leads)
    # skipna=False: years without a start (Feb/Aug: all but every fifth year) stay NaN, not zero
    u = (u * xr.DataArray(DJF_DAYS, dims="lead", coords={"lead": u["lead"]})).sum("lead", skipna=False) / DJF_DAYS.sum()
    return (u - u.mean()).values.ravel()


def tercile_bounds(init, k):
    # (lo, hi, kind): matched to the start and lead when the Reanalysis has the start month, else a surrogate
    if init in REF_ALL:
        a, kind = ref_djf_anomalies(init, k), "matched" if init in ref else "matched (Feb/Aug starts every fifth year)"
    else:
        parts = [ref_djf_anomalies(m, k) for m in ref]
        a, kind = np.concatenate([p for p in parts if p is not None]), "surrogate (May+Nov pooled)"
    a = a[np.isfinite(a)]
    return *np.percentile(a, [100 / 3, 200 / 3]), kind, a.size


# observed DJF anomalies: calendar-month climatology subtracted month by month
era5_month_clim = e_u60_month.sel(time=slice(f"{CLIM_YEARS[0]}-01-01", f"{CLIM_YEARS[1]}-12-31")).groupby("time.month").mean("time")
e_anom_month = (e_u60_month.groupby("time.month") - era5_month_clim).drop_vars("month")
assert e_anom_month.dims == ("time",), e_anom_month.dims
e_djf = pd.Series({w: float(v) for w in range(ERA5_FIRST_YEAR, YEAR) if (v := djf_mean(e_anom_month, w)) is not None}).dropna()
e_lo, e_hi = np.percentile(e_djf.values, [100 / 3, 200 / 3])

# drift-correction alternatives for the DJF climatology (Feb/Aug); the May climatology is direct
CLIM_VARIANTS = {"direct": ("direct", list(REF_ALL)), "lead_drift": ("lead_drift", list(ref)),
                 "target_month": ("target_month", list(ref)),
                 "lead_drift, May drift only": ("lead_drift", [5]), "lead_drift, Nov drift only": ("lead_drift", [11])}

m_djf, prob_rows = {}, []
for init in INIT_MONTHS:
    raw = mm_u60[init].assign_coords(time=monthly[init]["time"]).swap_dims(lead="time")
    for k, w in enumerate((YEAR, YEAR + 1)):
        v_raw = djf_mean(raw, w)
        if v_raw is None:
            continue
        lo, hi, kind, n_ref = tercile_bounds(init, k)
        for name, (method, refs) in CLIM_VARIANTS.items():
            if init in ref and name != "direct":
                continue                                           # May/Nov: every method is the direct climatology
            lc = sc.lead_climatology(init, {m: REF_ALL[m] for m in refs}, e_month_clim_m, method=method)
            clim = sc.lat_mean(lc["U"].sel(plev=REF_LEVEL), sc.U60).assign_coords(time=("lead", monthly[init]["time"].values))
            v = v_raw - djf_mean(clim.swap_dims(lead="time"), w)
            if name == CLIM_METHOD:
                m_djf[(init, w)] = v
            prob_rows.append({"start": MONTH_NAMES[init], "winter": f"{w}/{str(w + 1)[2:]}", "climatology": name,
                              "members": int(v.size), "mean anomaly (m/s)": float(v.mean()),
                              "P(weak)": float((v < lo).mean()), "P(normal)": float(((v >= lo) & (v <= hi)).mean()),
                              "P(strong)": float((v > hi).mean()), "bounds": f"{lo:+.1f} / {hi:+.1f} ({kind}, n={n_ref})"})
probs = pd.DataFrame(prob_rows)
probs.to_csv(DATA_DIR / f"djf_vortex_terciles_{YEAR}.csv", index=False)
print(f"[WRITE] {DATA_DIR / f'djf_vortex_terciles_{YEAR}.csv'}; ERA5 tercile bounds {e_lo:+.1f} / {e_hi:+.1f} m/s")
nino_years = [w for w in e_djf.index if obs_djf_nino.get(w, np.nan) >= ENSO_THRESHOLD]
print(f"ERA5 El Niño winters: weak {np.mean(e_djf[nino_years] < e_lo):.0%}, strong {np.mean(e_djf[nino_years] > e_hi):.0%} "
      f"({len(nino_years)} winters)")
print(f"Uncalibrated ensemble tercile probabilities ({CLIM_METHOD} climatology marks the plotted anomalies; Feb/Aug rows compare the climatology methods):")
display(probs.round(2))

fig, ax = plt.subplots(figsize=(13, 5))
pos, labels = [], []
for k, ((init, w), v) in enumerate(m_djf.items()):
    xs = np.full(v.size, k) + rng.uniform(-0.15, 0.15, v.size)
    ax.scatter(xs, v.values, color=INIT_COLORS[init], s=18, alpha=0.7)
    ax.plot([k - 0.3, k + 0.3], [float(v.mean())] * 2, color=INIT_COLORS[init], lw=3)
    lo, hi, kind, _ = tercile_bounds(init, w - YEAR)
    ax.plot([k - 0.35, k + 0.35], [lo, lo], color="0.4", lw=1, ls="--" if kind.startswith("surrogate") else "-")
    ax.plot([k - 0.35, k + 0.35], [hi, hi], color="0.4", lw=1, ls="--" if kind.startswith("surrogate") else "-")
    pos.append(k)
    labels.append(f"{MONTH_NAMES[init]} start\n{w}/{str(w + 1)[2:]}")
k0 = len(pos) + 0.5
xs = np.full(e_djf.size, k0) + rng.uniform(-0.15, 0.15, e_djf.size)
nino_mask = np.array([w in nino_years for w in e_djf.index])
ax.scatter(xs[~nino_mask], e_djf.values[~nino_mask], color="0.5", s=18, label="ERA5 winters")
ax.scatter(xs[nino_mask], e_djf.values[nino_mask], color="tab:orange", edgecolor="k", s=40, label="ERA5 El Niño winters")
ax.plot([k0 - 0.35, k0 + 0.35], [e_lo, e_lo], color="0.4", lw=1)
ax.plot([k0 - 0.35, k0 + 0.35], [e_hi, e_hi], color="0.4", lw=1)
pos.append(k0)
labels.append(f"ERA5\n{int(e_djf.index[0])}-{int(e_djf.index[-1]) + 1}")
ax.axhline(0, color="k", lw=0.8)
ax.set_xticks(pos, labels, fontsize=fontz * 0.75)
ax.set_ylabel("DJF U(60°N, 10 hPa) anomaly (m s$^{-1}$)")
ax.set_title("DJF polar-vortex strength (drift-corrected; grey: tercile bounds, dashed = surrogate)", loc="left")
ax.legend(frameon=False, loc="lower left", fontsize=fontz * 0.7)
ax.grid(alpha=0.3, axis="y")
save(fig, "djf", "vortex", "terciles", YEAR)
plt.show()
