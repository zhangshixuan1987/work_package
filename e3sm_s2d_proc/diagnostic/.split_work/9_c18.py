rng = np.random.default_rng(0)


def freq_row(label, counts, kind, **extra):
    counts = np.asarray(counts, dtype=float)
    boots = rng.choice(counts, size=(5000, counts.size)) if counts.size else np.full((1, 1), np.nan)
    return {"set": label, **extra, "wind": kind, "winters": counts.size, "events": int(np.nansum(counts)),
            "events per winter": counts.mean() if counts.size else np.nan,
            "5%": np.percentile(boots.mean(1), 5), "95%": np.percentile(boots.mean(1), 95),
            "P(≥1 SSW)": (counts > 0).mean() if counts.size else np.nan,
            "P 5%": np.percentile((boots > 0).mean(1), 5), "P 95%": np.percentile((boots > 0).mean(1), 95)}


def mid_month(times):
    # the 15th of each time's month (datetime64[ns])
    return sc.to_datetime64(times).astype("datetime64[M]").astype("datetime64[ns]") + np.timedelta64(14, "D")


def case_of(member):
    return str(member).split(":")[0]


# observed
obs_events = sc.ssw_table(e_u60.values, e_u60["time"].values, "ERA5", **SSW)
obs_winters = sc.covered_winters(e_u60["time"].values)
obs_counts = pd.Series(0, index=obs_winters).add(obs_events.groupby("winter").size(), fill_value=0).astype(int)
nino_obs = rc.box_mean(rc.obs_anomaly("TS", years=(ERA5_FIRST_YEAR, YEAR), mode=derivation_mode, product=NINO_PRODUCT_OBS),
                       rc.NINO_REGIONS["Nino3.4"]).load()
nino_obs = nino_obs.assign_coords(time=mid_month(nino_obs["time"].values))


def winter_mean(series, w, months=DJF):
    times = [np.datetime64(f"{w + dy}-{m:02d}-15") for dy, m in months]
    try:
        v = series.sel(time=times)
    except KeyError:
        return np.nan
    return float(v.mean()) if np.isfinite(v).all() else np.nan


obs_djf_nino = pd.Series({w: winter_mean(nino_obs, w) for w in obs_winters})
rows = [freq_row(f"ERA5 {obs_winters[0]}/{str(obs_winters[0] + 1)[2:]}-{obs_winters[-1]}/{str(obs_winters[-1] + 1)[2:]}",
                 obs_counts.values, "ERA5"),
        freq_row("ERA5 El Niño winters", obs_counts[obs_djf_nino >= ENSO_THRESHOLD].values, "ERA5"),
        freq_row("ERA5 La Niña winters", obs_counts[obs_djf_nino <= -ENSO_THRESHOLD].values, "ERA5")]
print("ERA5 El Niño winters:", ", ".join(f"{w}/{str(w + 1)[2:]}" for w in obs_counts[obs_djf_nino >= ENSO_THRESHOLD].index))

# model: per start x case x winter, then per start (cases pooled, member-count weights), raw and bias-adjusted
model_counts = {}                      # (init, winter, kind) -> per-member counts (all cases)
for init in INIT_MONTHS:
    winters = sc.covered_winters(m_u60[init]["time"].values)
    members = [str(m) for m in m_u60[init]["member"].values]
    for kind, tab in (("raw", model_ssw[init]), ("bias-adjusted", model_ssw_adj[init])):
        for w in winters:
            c = tab[tab["winter"] == w].groupby("series").size().reindex(members, fill_value=0)
            model_counts[(init, w, kind)] = c
            label = f"E3SM {MONTH_NAMES[init]} start, {w}/{str(w + 1)[2:]}"
            for case in dict.fromkeys(case_of(m) for m in members):
                cc = c[[case_of(m) == case for m in c.index]]
                rows.append(freq_row(label, cc.values, kind, case=case))
            weights = ", ".join(f"{k} {v}" for k, v in pd.Series([case_of(m) for m in members]).value_counts().items())
            rows.append(freq_row(label, c.values, kind, case=f"pooled ({weights})"))
for kind in ("raw", "bias-adjusted"):
    for w in (YEAR, YEAR + 1):
        pooled = [c.values for (i, ww, k), c in model_counts.items() if ww == w and k == kind]
        if pooled:
            rows.append(freq_row(f"E3SM all starts, {w}/{str(w + 1)[2:]}", np.concatenate(pooled), kind, case="pooled"))
freq = pd.DataFrame(rows)
freq.to_csv(DATA_DIR / f"ssw_frequency_{YEAR}.csv", index=False)
print(f"[WRITE] {DATA_DIR / f'ssw_frequency_{YEAR}.csv'}")
display(freq[freq["wind"] != "bias-adjusted"].round(2))
print("Threshold sensitivity (bias-adjusted wind), pooled per start:")
display(freq[(freq["wind"] == "bias-adjusted") & freq["case"].astype(str).str.startswith("pooled")].round(2))

# event catalogue with per-event properties (raw wind primary)
cat = [sc.ssw_event_details(e_u60.values, e_u60["time"].values, obs_events["central_date"], SSW["recover_days"])
       .assign(source="ERA5", series="ERA5", wind="ERA5")]
for i in INIT_MONTHS:
    for kind, tab, src in (("raw", model_ssw[i], m_u60[i]), ("bias-adjusted", model_ssw_adj[i], m_u60_adj[i])):
        for mem, ev in tab.groupby("series"):
            cat.append(sc.ssw_event_details(src.sel(member=mem).values, src["time"].values, ev["central_date"],
                                            SSW["recover_days"]).assign(source=f"E3SM init{i:02d}", series=mem, wind=kind))
events_all = pd.concat(cat, ignore_index=True)
events_all["winter"] = [sc.winter_of(d) for d in events_all["central_date"]]
events_all["case"] = [case_of(s) for s in events_all["series"]]
events_all.to_csv(DATA_DIR / f"ssw_central_dates_{YEAR}.csv", index=False)
print(f"[WRITE] {DATA_DIR / f'ssw_central_dates_{YEAR}.csv'}")
display(events_all[events_all["wind"] != "bias-adjusted"].groupby(["source", "winter"])
        [["min_u", "easterly_days", "days_to_recovery"]].median().round(1))
