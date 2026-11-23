#!/usr/bin/env python3
"""Minor revision: the derived charge range as the primary instance, and the per-window
admission criterion for the MPC comparisons.

  python3 rev_minor.py rules   -> T38_derived_rules / _perwindow
        TUNED and the guard refitted at c_sw = 0.025 and 0.114 SAR (in-force tariff, the three
        cities) under Rule A4 with the pooled 0.03 K criterion (as T1), and TUNED refitted under
        the per-window 0.03 K criterion of t1_sens at the same charges; scored by Rule C on the
        seven held-out windows (t1_v2.evaluate). The 0.15 SAR rows are those of T1v3 / T1v3_sens.
  python3 rev_minor.py stats   -> T39_derived_stats
        The pre-declared families F1-F4 (stats_v3.contrast, Holm within each three-city family)
        at 0.025 and 0.114 SAR, and F2/F3 with every arm admitted under the per-window criterion
        at 0.025, 0.114 and 0.15 SAR. The 0.15 SAR pooled-criterion families are recomputed from
        the headline files and checked against T20v3_paired_stats.
  python3 rev_minor.py reprice -> T41_ppo_repriced
        The learned and shielded arms re-priced at c_sw = 0.025 and 0.114 SAR: their held-out
        trajectories (PPO policies trained, ceilings and shield margins fitted at 0.15 SAR) are kept
        and each window's cost is re-priced, J(c) = J + (c - 0.15) f_s, for the arm and for its
        equalised THERM. Exact for a fixed trajectory (TUNED reproduces T38 to 1e-6); a window whose
        comparator has several crossings with different J keeps the crossing chosen at 0.15 SAR and
        is counted in n_multi.
  RBRL_TAG=<tag> [RBRL_CTH=.. | RBRL_DT_MIN=..] python3 rev_minor.py sweep -> T40_derived_sweep<tag>
        TUNED refitted by Rule A4 at c_sw = 0.025 and 0.114 SAR in one scenario of the sensitivity
        tables (thermal capacity or control interval from the environment); with no tag also the
        2027 tariff and the THERM deadbands 0.50, 0.75 and 1.25 C (comparator refitted as in T7).
"""
from __future__ import annotations
import os, sys, json, time
sys.path.insert(0, "/opt/rbrl")
import numpy as np
import config as CFG
import final_tables as F
import t1_v2 as T1V2
import t1_sens as SENS
import stats_v3 as SV
import v3_names as V
import rbrl_io as IO

PRICE = CFG.PRICE
DERIVED = (0.025, 0.114)
REF = CFG.C_SW                      # 0.15
RULES = "T38_derived_rules"
STATS = "T39_derived_stats"


def run_rules():
    summ, perw = [], []
    c0 = F.C_SW
    try:
        for c in DERIVED:
            F.C_SW = c                                   # read at call time by every run and fit
            for city in F.CITIES:
                t0 = time.time()
                th = {w: F.therm_run(city, w, PRICE) for w in F.WINDOWS}
                arms = []
                bt = T1V2.fit_tuned(city, PRICE, th)
                if bt:
                    arms.append(("TUNED", "pooled", dict(on=round(float(bt[1]), 3), off=round(float(bt[2]), 3)),
                                 bt.fit_status))
                bg, gst = T1V2.fit_guard(city, PRICE, th)
                if bg:
                    cap = bg[1]
                    arms.append(("guard", "pooled", dict(cap=cap, on=round(cap - T1V2.GUARD_HOT, 3),
                                                         off=T1V2.GUARD_OFF), gst))
                p, st, _ = SENS.tuned_per_window(city, PRICE)
                if p is not None:
                    arms.append(("TUNED", "per_window", dict(on=p[0], off=p[1]), st))
                for arm, adm, prm, st in arms:
                    runs = {w: F.rule_run(city, w, PRICE, prm["on"], prm["off"]) for w in F.EVAL_W}
                    s, pw = T1V2.evaluate(city, PRICE, runs, th)
                    summ.append(dict(c_sw=c, city=city, arm=arm, admission=adm, params=json.dumps(prm),
                                     fit_status=st, **s))
                    perw += [dict(c_sw=c, city=city, arm=arm, admission=adm, fit_status=st, **r) for r in pw]
                    print(f"c_sw={c:g} {city:7s} {arm:5s} {adm:10s} {prm} {st} s_eq={s['s_pooled_eq']:+.2f}% "
                          f"dT={s['dT_pooled']:+.4f} starts={s['start_sav_eq']:+.1f}% "
                          f"E={s['energy_sav_eq']:+.2f}% over={s['deg_h_over']:.4f}", flush=True)
                print(f"  ({time.time() - t0:.0f}s)", flush=True)
    finally:
        F.C_SW = c0
    T1V2.write(RULES, summ)
    T1V2.write(RULES + "_perwindow", perw)


# ------------------------------------------------------------------ per-window rows
def _by_window(rows, what):
    out = {}
    for r in rows:
        w = int(float(r["window"]))
        if w in out:
            raise ValueError(f"{what}: window {w} twice")
        out[w] = r
    if len(out) != len(F.EVAL_W):
        raise ValueError(f"{what}: {len(out)} windows")
    return out


def rule_pw(arm, city, c, adm):
    """{w: row} of TUNED / guard at charge c under admission adm."""
    if abs(c - REF) < 1e-12:
        if adm == "pooled":
            return V.perwindow(arm, city)
        rows = [r for r in IO.read_csv(V.name("T1_SENS") + "_perwindow")
                if r["city"] == city and r["arm"] == arm and r["variant"] == "per_window_0.03"]
        return _by_window(rows, f"{arm} {city} per-window 0.15")
    rows = [r for r in IO.read_csv(RULES + "_perwindow")
            if r["city"] == city and r["arm"] == arm and r["admission"] == adm and abs(float(r["c_sw"]) - c) < 1e-12]
    return _by_window(rows, f"{arm} {city} {adm} {c:g}")


def mpc_pw(forecast, city, c, adm):
    sfx = ("" if abs(c - REF) < 1e-12 else f"_csw{c:g}") + ("_pwadm" if adm == "per_window" else "")
    stem = V.mpc_stem(forecast, False, False, sfx)
    rows = [r for r in V.read_stem(f"{stem}_{city}_perwindow") if abs(float(r["horizon_h"]) - 24.0) < 1e-9
            and r.get("forecast", "perfect") == forecast]
    return _by_window(rows, f"{stem}_{city}")


def arm(name, city, c, adm):
    pw = (rule_pw(name, city, c, adm) if name in ("TUNED", "guard")
          else mpc_pw("persistence" if name == "MPCp24" else "perfect", city, c, adm))
    fit = "/".join(sorted({r.get("fit_status", "") for r in pw.values()}))
    return {w: V.s_w(r) for w, r in pw.items()}, V.pooled(pw.values()), fit, pw


def _admitted(name, city, c, adm):
    """False when the arm has no admissible fit at this charge and criterion (mpc_lambda writes no
    file then; its log records 'no admissible lam')."""
    try:
        arm(name, city, c, adm)
        return True
    except FileNotFoundError:
        return False


def family(fam, a, b, c, adm):
    block = []
    for city in F.CITIES:
        miss = [x for x in (a, b) if x is not None and not _admitted(x, city, c, adm)]
        if miss:
            block.append(dict(c_sw=c, admission=adm, family=fam, contrast=f"{SV.LABEL[a]} vs {SV.LABEL[b]}",
                              arm_a=a, arm_b=b or "THERM", city=city, admitted=False,
                              not_admitted="+".join(miss)))
            continue
        sa, pa, fa, pwa = arm(a, city, c, adm)
        if b is None:
            sb, pb, fb, pwb = None, 0.0, "-", None
        else:
            sb, pb, fb, pwb = arm(b, city, c, adm)
        row = dict(c_sw=c, admission=adm, family=fam, contrast=f"{SV.LABEL[a]} vs {SV.LABEL[b]}",
                   arm_a=a, arm_b=b or "THERM", city=city, fit_status_a=fa, fit_status_b=fb,
                   dT_a=float(np.mean([float(r["dT"]) for r in pwa.values()])) if "dT" in next(iter(pwa.values())) else "",
                   admitted=True, not_admitted="", **SV.contrast(sa, sb, pa, pb))
        block.append(row)
    ok = [r for r in block if r["admitted"]]          # Holm over the cities with an admitted pair
    adj = SV.holm([r["p_wilcoxon"] for r in ok])
    for r, p in zip(ok, adj):
        r.update(holm=True, family_size=len(ok), p_adj=float(p), p_floor=SV.p_floor(r["n"], len(ok)))
    return block


FAMS = (("F1", "TUNED", None), ("F2", "MPC24", "TUNED"), ("F3", "MPCp24", "TUNED"), ("F4", "guard", "TUNED"))


def run_stats():
    rows = []
    for c in DERIVED + (REF,):
        for fam, a, b in FAMS:
            rows += family(fam, a, b, c, "pooled")
    for c in DERIVED + (REF,):
        for fam, a, b in FAMS[1:3]:
            rows += family(fam, a, b, c, "per_window")
    # the 0.15 SAR pooled families must reproduce the published T20v3_paired_stats
    ref = {(r["family"], r["city"]): r for r in V.read("STATS")}
    for r in rows:
        if r["admission"] == "pooled" and abs(r["c_sw"] - REF) < 1e-12 and r["admitted"]:
            q = ref[(r["family"], r["city"])]
            for k in ("hl", "hl_lo", "hl_hi", "p_wilcoxon", "p_adj", "d_pooled"):
                if abs(float(q[k]) - float(r[k])) > 1e-4:
                    raise SystemExit(f"0.15 check failed: {r['family']} {r['city']} {k} {q[k]} vs {r[k]}")
    print("0.15 SAR pooled families reproduce T20v3_paired_stats", flush=True)
    IO.write_csv(STATS, rows, ndigits=6)
    for r in rows:
        if not r["admitted"]:
            print(f"c_sw={r['c_sw']:<5g} {r['admission']:10s} {r['family']} {r['city']:7s} NOT ADMITTED: "
                  f"{r['not_admitted']}", flush=True)
            continue
        print(f"c_sw={r['c_sw']:<5g} {r['admission']:10s} {r['family']} {r['contrast']:30s} {r['city']:7s} "
              f"HL={r['hl']:+.2f} [{r['hl_lo']:+.2f},{r['hl_hi']:+.2f}] +{r['n_pos']}/-{r['n_neg']} "
              f"p_adj={r['p_adj']:.3f} pooled {r['s_pooled_a']:.1f} vs {r['s_pooled_b']:.1f} "
              f"d={r['d_pooled']:+.2f} {r['resolution']} fit {r['fit_status_a']}|{r['fit_status_b']}", flush=True)



# ------------------------------------------------------------------ re-pricing of fixed trajectories
REPR = "T41_ppo_repriced"
SWEEP = "T40_derived_sweep"


def _jc(r, c):
    k = c - REF
    return (float(r["J_arm"]) + k * float(r["fs_arm"]), float(r["J_eq"]) + k * float(r["fs_eq"]))


def pooled_at(rows, c):
    ja = je = 0.0
    for r in rows:
        a, e = _jc(r, c)
        ja += a
        je += e
    return 100.0 * (1.0 - ja / je)


def _multi(rows):
    return sum(1 for r in rows if r.get("J_eq_max") not in (None, "")
               and abs(float(r["J_eq_max"]) - float(r["J_eq"])) > 1e-9)


def _d1_ppo(city):
    """{seed: {w: row}} of the D1 shield on each saved PPO model, the sigma chosen by A4 on the
    tuning weeks (lowest tune_obj among fitted rows), as for the latch (v3_names.d1_selected)."""
    summ = [r for r in IO.read_csv(f"T_D1_ppo_{city}") if r.get("fit_status") not in ("", "none")]
    pw = IO.read_csv(f"T_D1_ppo_{city}_perwindow")
    out = {}
    for arm in sorted({r["arm"] for r in summ}):
        best = min((r for r in summ if r["arm"] == arm), key=lambda r: float(r["tune_obj"]))
        out[arm] = _by_window([r for r in pw if r["arm"] == arm and abs(float(r["sigma"]) - float(best["sigma"])) < 1e-12],
                              f"D1-PPO {city} {arm}")
    return out


def _norules(city):
    rows = [r for r in IO.read_csv(V.name("NORULES_PW")) if r.get("city", city) == city]
    out = {}
    for r in rows:
        out.setdefault(str(r.get("seed", "")), []).append(r)
    return {k: _by_window(v, f"NORULES {city} {k}") for k, v in out.items()}


def run_reprice():
    rows = []
    for city in F.CITIES:
        arms = {"TUNED": {"-": V.perwindow("TUNED", city)}, "guard": {"-": V.perwindow("guard", city)},
                "RBRL": V.perwindow_seeds("RBRL", city), "D1_latch": {"-": V.d1_selected(city)[1]},
                "D1_ppo": _d1_ppo(city)}
        try:
            arms["NORULES"] = _norules(city)
        except (FileNotFoundError, ValueError, LookupError) as e:
            print(f"NORULES {city}: {e}", flush=True)
        for arm, seeds in arms.items():
            if not seeds:                  # no admissible seed in this city (PPO without rules)
                continue
            for c in DERIVED + (REF,):
                per = {sd: pooled_at(pw.values(), c) for sd, pw in seeds.items()}
                nm = sum(_multi(pw.values()) for pw in seeds.values())
                for sd, v in per.items():
                    rows.append(dict(c_sw=c, city=city, arm=arm, seed=sd, s_pooled_eq=v, n_multi=_multi(seeds[sd].values())))
                vals = list(per.values())
                rows.append(dict(c_sw=c, city=city, arm=arm, seed="all", s_pooled_eq=float(np.mean(vals)),
                                 s_min=min(vals), s_max=max(vals), n_seeds=len(vals), n_multi=nm))
                print(f"c={c:<5g} {city:7s} {arm:8s} seeds={len(vals)} mean={np.mean(vals):6.2f} "
                      f"range={min(vals):6.2f}..{max(vals):6.2f} n_multi={nm}", flush=True)
    # TUNED must reproduce the refits of T38 exactly
    t38 = {(r["city"], float(r["c_sw"])): float(r["s_pooled_eq"]) for r in IO.read_csv(RULES)
           if r["arm"] == "TUNED" and r["admission"] == "pooled"}
    for r in rows:
        if r["arm"] == "TUNED" and r["seed"] == "all" and r["c_sw"] in DERIVED:
            if abs(t38[(r["city"], r["c_sw"])] - r["s_pooled_eq"]) > 1e-3:
                raise SystemExit(f"re-pricing check failed: TUNED {r['city']} {r['c_sw']}")
    print("TUNED re-priced reproduces T38", flush=True)
    IO.write_csv(REPR, rows, ndigits=6)


def run_sweep():
    tag = os.environ.get("RBRL_TAG", "")
    settings = [("inforce", 1.0)] + ([("y2027", 1.0), ("inforce", 0.50), ("inforce", 0.75), ("inforce", 1.25)]
                                     if not tag else [])
    rows = []
    c0 = F.C_SW
    try:
        for c in DERIVED:
            F.C_SW = c
            for tname, db in settings:
                pv = F.TARIFFS[tname]
                for city in F.CITIES:
                    t0 = time.time()
                    th = {w: T1V2.therm_d(city, w, pv, db) for w in F.WINDOWS} if db != 1.0 else \
                         {w: F.therm_run(city, w, pv) for w in F.WINDOWS}
                    bt = T1V2.fit_tuned(city, pv, th, db=db)
                    base = dict(c_sw=c, city=city, tariff=tname, deadband=db, cth=CFG.C_TH, dt_min=CFG.DT_MIN,
                                therm_fs=float(np.mean([th[w]["f_s"] for w in F.EVAL_W])))
                    if not bt:
                        rows.append(dict(base, fit_status="none"))
                        print(f"c={c:g} {city} {tname} db={db} cth={CFG.C_TH:g} dt={CFG.DT_MIN:g}: no admissible rule", flush=True)
                        continue
                    runs = {w: F.rule_run(city, w, pv, bt[1], bt[2]) for w in F.EVAL_W}
                    s, _ = T1V2.evaluate(city, pv, runs, th, db=db)
                    rows.append(dict(base, on=round(float(bt[1]), 3), off=round(float(bt[2]), 3),
                                     fit_status=bt.fit_status, **s))
                    print(f"c={c:g} {city:7s} {tname:7s} db={db:.2f} cth={CFG.C_TH:g} dt={CFG.DT_MIN:g} "
                          f"rule {bt[1]:.2f}/{bt[2]:.2f} {bt.fit_status} s_eq={s['s_pooled_eq']:+.2f}% "
                          f"({time.time() - t0:.0f}s)", flush=True)
    finally:
        F.C_SW = c0
    T1V2.write(SWEEP, rows)

if __name__ == "__main__":
    {"rules": run_rules, "stats": run_stats, "reprice": run_reprice, "sweep": run_sweep}[sys.argv[1]]()
