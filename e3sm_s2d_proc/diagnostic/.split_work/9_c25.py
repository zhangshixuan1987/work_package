COMPOSITE_WINDOW = (-30, 60)          # days that must all be in the record
PRE_WINDOW, PCH_WINDOWS = (-20, -1), [(-20, -1), (0, 14), (15, 29), (30, 59)]
PCH_LEVELS = [10.0, 50.0, 100.0, 300.0, 500.0]
NBOOT = 2000

# polar-cap SLP: ERA5 and model (drift-corrected PSL climatology), standardized by the ERA5 daily sigma
e_slp = sc.era5_polar_cap_slp(ERA5_DAILY_DIR, years=(ERA5_FIRST_YEAR, YEAR))
slp_clim, slp_std = sc.daily_climatology(e_slp, CLIM_YEARS)
e_slp_s = (e_slp - on_days(slp_clim, e_slp)) / on_days(slp_std, e_slp)
m_slp_s = {}
for init in INIT_MONTHS:
    d = daily[init]
    pc = sc.lat_mean(rc.matched_climatology("PSL", init, method=CLIM_METHOD, mode=derivation_mode).mean("lon"), sc.PCAP)
    clim = sc.daily_from_lead(pc, YEAR, init, d["time"].values).assign_coords(time=d["time"])
    raw = sc.lat_mean(d["PSL"], sc.PCAP)
    m_slp_s[init] = (raw - clim) / on_days(slp_std, raw)

lags = np.arange(COMPOSITE_WINDOW[0], COMPOSITE_WINDOW[1] + 1)


def windows(series, dates):
    return sc.lag_composite(series, dates, lags)


# ERA5 composites, clusters (winters) and season-matching weights
e_dates = list(obs_events["central_date"])
e_comp = {"u": windows(e_u60, e_dates), "vt": windows(e_vt, e_dates), "pch": windows(e_pch, e_dates), "slp": windows(e_slp_s, e_dates)}
e_full = np.isfinite(e_comp["u"].values).all(axis=1)
e_months = np.array([d.month for d in e_dates])
e_clusters = np.array([sc.winter_of(d) for d in e_dates])

# model composites (raw-wind events), clusters (start x member trajectory)
m_parts, m_months, m_clusters = {k: [] for k in e_comp}, [], []
for init in INIT_MONTHS:
    for _, ev in model_ssw[init].iterrows():
        mem, d = ev["series"], ev["central_date"]
        for key, src in (("u", m_u60[init]), ("vt", m_vt[init]), ("pch", m_pch[init]), ("slp", m_slp_s[init])):
            m_parts[key].append(windows(src.sel(member=mem), [d]))
        m_months.append(d.month)
        m_clusters.append(f"{init}:{mem}")
m_comp = {k: xr.concat(v, dim="event") for k, v in m_parts.items()} if m_months else None
m_months, m_clusters = np.array(m_months), np.array(m_clusters)
m_full = np.isfinite(m_comp["u"].values).all(axis=1) if m_comp else np.array([], bool)
print(f"events with the whole {COMPOSITE_WINDOW} window: ERA5 {e_full.sum()} of {e_full.size}, "
      f"E3SM {m_full.sum()} of {m_full.size}")

# season-matching weights of the ERA5 events (model onset-month frequency / ERA5 onset-month frequency)
p_mod = pd.Series(m_months[m_full]).value_counts(normalize=True)
p_obs = pd.Series(e_months[e_full]).value_counts(normalize=True)
e_w = np.array([p_mod.get(m, 0.0) / p_obs[m] if f else 0.0 for m, f in zip(e_months, e_full)])
missing = sorted(set(p_mod.index) - set(p_obs.index))
if missing:
    print(f"[note] model onset months with no ERA5 event: {missing}; the matched ERA5 composite cannot represent them")
print("onset months, model / ERA5 (fractions):")
display(pd.DataFrame({"E3SM": p_mod, "ERA5": p_obs}).fillna(0).sort_index().round(2).T)


def comp_stats(comp, full, weights, clusters, boot=True):
    # weighted mean, cluster-bootstrap 90% CI of the mean, and the 10-90% event spread
    x = comp.values[full]
    w = np.asarray(weights, float)[full]
    mean = np.nansum(x * w.reshape((-1,) + (1,) * (x.ndim - 1)), 0) / w.sum()
    if not boot:
        return mean
    cl = np.asarray(clusters)[full]
    labs = np.unique(cl)
    idx = {l: np.flatnonzero(cl == l) for l in labs}
    boot = []
    for _ in range(NBOOT):
        ev = np.concatenate([idx[l] for l in rng.choice(labs, labs.size)])
        boot.append(np.nansum(x[ev] * w[ev].reshape((-1,) + (1,) * (x.ndim - 1)), 0) / w[ev].sum())
    boot = np.asarray(boot)
    return mean, np.nanpercentile(boot, 5, 0), np.nanpercentile(boot, 95, 0), np.nanpercentile(x, 10, 0), np.nanpercentile(x, 90, 0)


ones_m = np.ones(m_full.size)
ones_e = np.where(e_full, 1.0, 0.0)
fig = plt.figure(figsize=(17, 10.5))
gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.18)
for k, (key, title, units) in enumerate([("u", f"(a) U(60°N, {REF_LEVEL:g} hPa)", "m s$^{-1}$"),
                                         ("vt", "(b) [v*T*] 100 hPa, 45-75°N (absolute)", "K m s$^{-1}$")]):
    ax = fig.add_subplot(gs[0, k])
    for comp, full, w, cl, color, label in ((e_comp[key], e_full, e_w, e_clusters, "k", "ERA5, season-matched"),
                                            (m_comp[key], m_full, ones_m, m_clusters, "tab:red", "E3SM")):
        mu, lo, hi, p10, p90 = comp_stats(comp, full, w, cl)
        ax.fill_between(lags, lo, hi, color=color, alpha=0.2, lw=0)
        ax.plot(lags, mu, color=color, lw=2.5, label=f"{label} ({int(full.sum())} events)")
        ax.plot(lags, p10, color=color, lw=0.8, ls=":")
        ax.plot(lags, p90, color=color, lw=0.8, ls=":")
    mu_raw = comp_stats(e_comp[key], e_full, ones_e, e_clusters)[0]
    ax.plot(lags, mu_raw, color="k", lw=1, alpha=0.6, label="ERA5, all events (unweighted)")
    ax.axvline(0, color="0.4", lw=1)
    ax.axhline(0, color="0.4", lw=0.8)
    ax.set_xlabel("days from central date")
    ax.set_ylabel(units)
    ax.set_title(title, loc="left")
    ax.legend(frameon=False, fontsize=fontz * 0.7)
    ax.grid(alpha=0.3)
pch_lev = np.array([-2.5, -2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2, 2.5])
for k, (comp, full, w, label) in enumerate(((e_comp["pch"], e_full, e_w, "ERA5, season-matched"), (m_comp["pch"], m_full, ones_m, "E3SM"))):
    ax = fig.add_subplot(gs[1, k])
    mu = xr.DataArray(comp_stats(comp, full, w, None, boot=False), dims=("lag", "plev"),
                      coords={"lag": lags, "plev": comp["plev"].values})
    c = mu.sel(plev=slice(1, 1000)).dropna("plev", how="all")
    cf = ax.contourf(c["lag"], c["plev"], c.T, levels=pch_lev, cmap="RdBu_r", extend="both")
    ax.set_yscale("log")
    ax.set_ylim(1000, 1 if label.startswith("ERA5") else 10)
    ax.set_yticks([1000, 300, 100, 30, 10, 3, 1], ["1000", "300", "100", "30", "10", "3", "1"])
    ax.axvline(0, color="k", lw=1)
    ax.set_xlabel("days from central date")
    ax.set_ylabel("hPa")
    ax.set_title(f"({'cd'[k]}) {label}: polar-cap height anomaly / σ ({int(full.sum())} events)", loc="left", fontsize=fontz * 0.9)
fig.colorbar(cf, ax=fig.axes[2:], label="σ (ERA5 daily)", shrink=0.8, pad=0.01)
save(fig, "ssw", "composite", YEAR)
plt.show()
