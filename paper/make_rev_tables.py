#!/usr/bin/env python3
"""Minor revision tables, read from the stamped result files (no typed numbers):
  tables/tab_derived.tex     main text: headline arms at the derived charges 0.025 and 0.114 SAR
  tables/tab_s_derived.tex   supplement: paired statistics at the derived charges and under the
                             per-window admission criterion
Also prints every value quoted in the text.
"""
import csv, os, sys
RES = os.environ.get("RBRL_RESULTS", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results"))
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tables")
CITIES = ("riyadh", "jeddah", "taif")
CN = dict(riyadh="Riyadh", jeddah="Jeddah", taif="Taif")
DER = (0.025, 0.114)


def rd(stem):
    with open(os.path.join(RES, stem + ".csv")) as f:
        return [r for r in csv.DictReader(f) if not r.get(next(iter(r)), "").startswith("#")]


def sw(r):
    return 100.0 * (1.0 - float(r["J_arm"]) / float(r["J_eq"]))


def pooled(rows):
    return 100.0 * (1.0 - sum(float(r["J_arm"]) for r in rows) / sum(float(r["J_eq"]) for r in rows))


def fm(x, d=1, sign=False):
    s = f"{x:+.{d}f}" if sign else f"{x:.{d}f}"
    return s.replace("-", "-") if not s.startswith("-0.0") or float(s) != 0 else s.replace("-", "")


MARK = {"pooled_0.03": "", "pooled_0.10": "^{b}"}


def rule(arm, city, c, adm="pooled"):
    rs = [r for r in rd("T38_derived_rules") if r["city"] == city and r["arm"] == arm
          and r["admission"] == adm and abs(float(r["c_sw"]) - c) < 1e-12]
    assert len(rs) == 1, (arm, city, c, adm, len(rs))
    return rs[0]


def mpc(fc, city, c, adm="pooled"):
    sfx = ("" if abs(c - 0.15) < 1e-12 else f"_csw{c:g}") + ("_pwadm" if adm == "per_window" else "")
    stem = ("T21v3_pers" if fc == "persistence" else "T21v3") + sfx + f"_{city}"
    rs = [r for r in rd(stem) if abs(float(r["horizon_h"]) - 24) < 1e-9]
    assert len(rs) == 1, (stem, len(rs))
    return rs[0]


def bound(city, c):
    rs = [r for r in rd("T5v3_bound") if r["city"] == city and abs(float(r["c_sw"]) - c) < 1e-12]
    assert len(rs) == 1
    return float(rs[0]["bound_reported"])


def stats(c, adm, fam, city):
    rs = [r for r in rd("T39_derived_stats") if abs(float(r["c_sw"]) - c) < 1e-12 and r["admission"] == adm
          and r["family"] == fam and r["city"] == city]
    assert len(rs) == 1, (c, adm, fam, city, len(rs))
    return rs[0]


def sg(x, d=1):
    """signed value at d decimals, with zero printed unsigned."""
    t = f"{float(x):+.{d}f}"
    return t[1:] if float(t) == 0 else t


def cell(r):
    return f"${float(r['s_pooled_eq']):.1f}{MARK.get(r['fit_status'], '^{?}')}$"


def hlcell(r):
    res = "" if r["resolution"] == "resolved" else r"^{\mathrm{nr}}"
    v = float(r['hl'])
    t = "0.0" if abs(v) < 0.05 else f"{v:+.1f}"
    return f"${t}{res}$"


def tab_derived():
    L = [r"\begin{table*}[t]", r"\centering",
         r"\caption{Headline arms at the charges the wear estimate supports, $\csw=0.025$ and $0.114\,\SAR$ per "
         r"on-command (tariff in force, occupied band $[22,24]\,^\circ$C, $\Delta t=15\,$min, seven held-out windows). "
         r"Upper block: pooled equalised saving $s$ (\%), every arm refitted at each charge by Rule~A. Lower block: "
         r"Hodges--Lehmann estimate of the per-window difference (points) for the pre-declared families. "
         r"$^{b}$fitted at the pooled $0.10\,$K fallback; $^{\mathrm{nr}}$pooled difference under $3$ points, "
         r"not resolved. Intervals and $p$-values: Supplementary Table~S9; the $0.15\,\SAR$ sensitivity case: "
         r"Table~\ref{tab:main}.}",
         r"\label{tab:derived}", r"\small",
         r"\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}lrrrrrr@{}}", r"\toprule",
         r" & \multicolumn{2}{c}{Riyadh} & \multicolumn{2}{c}{Jeddah} & \multicolumn{2}{c}{Taif} \\",
         r"\cmidrule(lr){2-3}\cmidrule(lr){4-5}\cmidrule(l){6-7}",
         r"$\csw$ ($\SAR$) & $0.025$ & $0.114$ & $0.025$ & $0.114$ & $0.025$ & $0.114$ \\", r"\midrule"]
    rows = [("TUNED", lambda ci, c: cell(rule("TUNED", ci, c))),
            ("Guard hysteresis", lambda ci, c: cell(rule("guard", ci, c))),
            (r"MPC $24\,$h, perfect", lambda ci, c: cell(mpc("perfect", ci, c))),
            (r"MPC $24\,$h, persistence", lambda ci, c: cell(mpc("persistence", ci, c))),
            (r"$B^{\ast}$ ($V_N$)", lambda ci, c: f"$\\le{bound(ci, c):.1f}$")]
    for name, fn in rows:
        L.append(name + " & " + " & ".join(fn(ci, c) for ci in CITIES for c in DER) + r" \\")
    L.append(r"\midrule")
    for fam, lab in (("F1", "F1: TUNED $-$ THERM"), ("F2", r"F2: MPC $24\,$h, perfect $-$ TUNED"),
                     ("F3", r"F3: MPC $24\,$h, persistence $-$ TUNED"), ("F4", "F4: guard $-$ TUNED")):
        L.append(lab + " & " + " & ".join(hlcell(stats(c, "pooled", fam, ci)) for ci in CITIES for c in DER) + r" \\")
    L += [r"\bottomrule", r"\end{tabular*}", r"\end{table*}"]
    return "\n".join(L) + "\n"


def tab_s_derived():
    L = [r"\begin{table}[!ht]", r"\centering",
         r"\caption{Paired comparisons at the derived charges ($\csw=0.025$ and $0.114\,\SAR$, every arm refitted by "
         r"Rule~A) and under the per-window $0.03\,$K admission criterion (TUNED and both MPC arms refitted under it, "
         r"at $0.025$, $0.114$ and $0.15\,\SAR$). Per-window difference of the equalised savings, $n=7$ held-out "
         r"windows; Hodges--Lehmann estimate with its exact $95.3\%$ interval, windows with positive/negative "
         r"difference, Holm-adjusted exact Wilcoxon $p$ within each three-city family, pooled difference (points). "
         r"$^{b}$an arm fitted at the $0.10\,$K fallback. Pooled differences under $3$ points are not resolved.}",
         r"\label{tab:s-derived}", r"\scriptsize", r"\setlength{\tabcolsep}{3pt}",
         r"\begin{tabularx}{\textwidth}{@{}llXrrrr@{}}", r"\toprule",
         r"Criterion, $\csw$ & Family & City & HL [95.3\% CI] & $+/-$ & $p_{\mathrm{adj}}$ & Pooled $\Delta s$ \\",
         r"\midrule"]
    blocks = [("pooled", c, f) for c in DER for f in ("F1", "F2", "F3", "F4")] + \
             [("per_window", c, f) for c in DER + (0.15,) for f in ("F2", "F3")]
    lastkey = None
    for adm, c, fam in blocks:
        for i, ci in enumerate(CITIES):
            r = stats(c, adm, fam, ci)
            key = (adm, c)
            lab = ((r"Pooled, " if adm == "pooled" else r"Per-window, ") + f"${c:g}$") if key != lastkey else ""
            lastkey = key
            if r.get("admitted") == "False":
                L.append(f"{lab} & {fam if i == 0 else ''} & {CN[ci]} & \\multicolumn{{4}}{{l}}{{no admissible "
                         f"MPC fit under this criterion}} \\\\")
                continue
            b = "^{b}" if "pooled_0.10" in (r["fit_status_a"] + r["fit_status_b"]) else ""
            L.append(f"{lab} & {fam if i == 0 else ''} & {CN[ci]} & ${sg(r['hl'])}$ "
                     f"[${sg(r['hl_lo'])}$, ${sg(r['hl_hi'])}$] & {r['n_pos']}/{r['n_neg']} & "
                     f"${float(r['p_adj']):.3f}$ & ${sg(r['d_pooled'])}{b}$ "
                     f"({'res.' if r['resolution'] == 'resolved' else 'n.r.'}) \\\\")
        if (adm, c, fam) != blocks[-1]:
            pass
    L += [r"\bottomrule", r"\end{tabularx}", r"\end{table}"]
    return "\n".join(L) + "\n"


def report():
    print("== derived range")
    for c in DER:
        for ci in CITIES:
            t, g = rule("TUNED", ci, c), rule("guard", ci, c)
            mp, ms = mpc("perfect", ci, c), mpc("persistence", ci, c)
            print(f"c={c:g} {ci:7s} TUNED {float(t['s_pooled_eq']):.2f} on {float(t['start_sav_eq']):.1f} dT {float(t['dT_pooled']):+.3f} | "
                  f"guard {float(g['s_pooled_eq']):.2f} | MPCp {float(mp['s_pooled_eq']):.2f} {mp['fit_status']} dT {float(mp['dT_pooled']):+.3f} on {float(mp['start_sav_eq']):.1f} | "
                  f"MPCs {float(ms['s_pooled_eq']):.2f} {ms['fit_status']} dT {float(ms['dT_pooled']):+.3f} | B* {bound(ci, c):.1f}")
    print("== per-window admission")
    for c in DER + (0.15,):
        for ci in CITIES:
            t = (rule("TUNED", ci, c, "per_window") if c != 0.15 else
                 [r for r in rd("T1v3_sens") if r["city"] == ci and r["arm"] == "TUNED" and r["variant"] == "per_window_0.03"][0])
            out = f"c={c:g} {ci:7s} TUNEDpw {float(t['s_pooled_eq']):.2f} dT {float(t['dT_pooled']):+.3f}"
            for fc in ("perfect", "persistence"):
                try:
                    m = mpc(fc, ci, c, "per_window")
                    out += (f" | {fc} {float(m['s_pooled_eq']):.2f} {m['fit_status']} {m.get('knob')} {m.get('backoff')} "
                            f"dT {float(m['dT_pooled']):+.3f} on {float(m['start_sav_eq']):.1f} lam/l* {float(m['lam_over_lam_star']):.3f}")
                except FileNotFoundError:
                    out += f" | {fc} NOT ADMITTED"
            print(out)



def sweep(tag, city, c, tariff="inforce", db=1.0):
    rs = [r for r in rd("T40_derived_sweep" + tag) if r["city"] == city and abs(float(r["c_sw"]) - c) < 1e-12
          and r["tariff"] == tariff and abs(float(r["deadband"]) - db) < 1e-12]
    assert len(rs) == 1, (tag, city, c, tariff, db, len(rs))
    return rs[0]


def rep(arm, city, c):
    rs = [r for r in rd("T41_ppo_repriced") if r["arm"] == arm and r["city"] == city and r["seed"] == "all"
          and abs(float(r["c_sw"]) - c) < 1e-12]
    return rs[0] if rs else None


def _sv(r):
    v = float(r["s_pooled_eq"])
    t = f"{v:.1f}"
    return "$" + ("0.0" if float(t) == 0 else t) + "$"


def _rng(r):
    if r is None:
        return "---"
    lo, hi = float(r["s_min"]), float(r["s_max"])
    a, b = f"{lo:.1f}", f"{hi:.1f}"
    return f"${a}$" if a == b else f"${a}$--${b}$"


def tab_s_sweeps():
    L = [r"\begin{table}[!ht]", r"\centering",
         r"\caption{The tuned thermostat and the learned and shielded arms at the derived charges, $\csw=0.025$ and "
         r"$0.114\,\SAR$: pooled equalised saving $s$ (\%) on the held-out windows. (a) TUNED refitted by Rule~A at each "
         r"charge in each scenario of main Tables~8 and~9, one setting changed at a time from the headline case "
         r"(tariff in force, THERM deadband $1\,^\circ$C, $\Cth=4{,}800\,$kJ/K, $\Delta t=15\,$min). (b) The PPO policies, "
         r"their ceilings and the shield margins are those trained and fitted at $0.15\,\SAR$; their held-out "
         r"trajectories are re-priced at each charge, which is exact for a fixed trajectory (TUNED re-priced in this way "
         r"reproduces its refits). Seeded arms: range over the five seeds; PPO without rules has one admissible seed, "
         r"in Riyadh.}",
         r"\label{tab:s-sweeps}", r"\scriptsize", r"\setlength{\tabcolsep}{4pt}",
         r"\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}lrrrrrr@{}}", r"\toprule",
         r" & \multicolumn{2}{c}{Riyadh} & \multicolumn{2}{c}{Jeddah} & \multicolumn{2}{c}{Taif} \\",
         r"\cmidrule(lr){2-3}\cmidrule(lr){4-5}\cmidrule(l){6-7}",
         r"$\csw$ ($\SAR$) & $0.025$ & $0.114$ & $0.025$ & $0.114$ & $0.025$ & $0.114$ \\", r"\midrule",
         r"\multicolumn{7}{@{}l}{\emph{(a) TUNED, refitted}}\\"]
    scen = [("Headline case", "", "inforce", 1.0), ("Tariff $0.21\\,\\SAR$/kWh (2027)", "", "y2027", 1.0),
            ("THERM deadband $0.50\\,^\\circ$C", "", "inforce", 0.50), ("THERM deadband $0.75\\,^\\circ$C", "", "inforce", 0.75),
            ("THERM deadband $1.25\\,^\\circ$C", "", "inforce", 1.25), ("$\\Cth=2{,}400\\,$kJ/K", "_cth2400", "inforce", 1.0),
            ("$\\Cth=7{,}200\\,$kJ/K", "_cth7200", "inforce", 1.0), ("$\\Cth=9{,}600\\,$kJ/K", "_cth9600", "inforce", 1.0),
            ("$\\Delta t=5\\,$min", "_dt5", "inforce", 1.0), ("$\\Delta t=30\\,$min", "_dt30", "inforce", 1.0)]
    for lab, tag, tar, db in scen:
        L.append(lab + " & " + " & ".join(_sv(sweep(tag, ci, c, tar, db)) for ci in CITIES for c in DER) + r" \\")
    L.append(r"\midrule")
    L.append(r"\multicolumn{7}{@{}l}{\emph{(b) Learned and shielded arms, re-priced}}\\")
    for lab, arm in (("Guard hysteresis", "guard"), ("RBRL, 5 seeds", "RBRL"), ("D1 shield on latch", "D1_latch"),
                     ("D1 on PPO, 5 seeds", "D1_ppo"), ("PPO without rules", "NORULES")):
        L.append(lab + " & " + " & ".join(_rng(rep(arm, ci, c)) for ci in CITIES for c in DER) + r" \\")
    L += [r"\bottomrule", r"\end{tabular*}", r"\end{table}"]
    return "\n".join(L) + "\n"

if __name__ == "__main__":
    report()
    if "--write" in sys.argv:
        os.makedirs(OUT, exist_ok=True)
        for fn, f in (("tab_derived.tex", tab_derived), ("tab_s_derived.tex", tab_s_derived), ("tab_s_sweeps.tex", tab_s_sweeps)):
            open(os.path.join(OUT, fn), "w").write(f())
            print("wrote", fn)
