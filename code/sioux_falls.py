"""Exploratory Sioux Falls illustration: link activations of static deterministic
user equilibrium (BPR costs) along a demand multiplier theta.

Usage:
  python3 -I sioux_falls.py NET_TNTP TRIPS_TNTP OUT_JSON OUT_PDF [--quick]

Deterministic: fixed tolerances, fixed iteration caps, no randomness.
"""
import hashlib
import json
import math
import re
import sys
import time

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra

GAP_TOL = 1e-12
MAX_SWEEPS = 4000
INNER_ROUNDS = 3
USED_FRAC = 1e-9          # link counted as used if v_a > USED_FRAC * total demand
THETA_GRID = [round(0.05 + 0.01 * i, 2) for i in range(146)]
N_ACT = 3                 # activations analysed in detail (nearest to theta=1)
BISECT_WIDTH = 1e-8
EXP_OFFSETS = [1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2]
FD_STEPS = [1e-4, 3e-5, 1e-5]
FD_CONTROL_OFFSET = 5e-4
HIDDEN_OFFSETS = [-1e-3, -3e-4, -1e-4, 1e-4, 3e-4, 1e-3, 3e-3]
LOOSE_TOLS = [1e-4, 1e-5]
FW_MAX_ITERS = 30000


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        h.update(f.read())
    return h.hexdigest()


def read_net(path):
    with open(path) as f:
        lines = f.read().splitlines()
    start = next(i for i, l in enumerate(lines) if "<END OF METADATA>" in l)
    rows = []
    for l in lines[start + 1:]:
        s = l.strip()
        if not s or s.startswith("~"):
            continue
        parts = s.replace(";", " ").split()
        if len(parts) < 7:
            continue
        rows.append(parts)
    net = {
        "init": [int(r[0]) - 1 for r in rows],
        "term": [int(r[1]) - 1 for r in rows],
        "cap": [float(r[2]) for r in rows],
        "t0": [float(r[4]) for r in rows],
        "B": [float(r[5]) for r in rows],
        "power": [float(r[6]) for r in rows],
    }
    assert len(rows) == 76, len(rows)
    return net


def read_trips(path):
    with open(path) as f:
        txt = f.read()
    body = txt.split("<END OF METADATA>", 1)[1]
    od = {}
    blocks = re.split(r"Origin\s+(\d+)", body)
    for i in range(1, len(blocks), 2):
        o = int(blocks[i]) - 1
        for d, q in re.findall(r"(\d+)\s*:\s*([0-9.eE+-]+)\s*;", blocks[i + 1]):
            q = float(q)
            d = int(d) - 1
            if q > 0 and d != o:
                od[(o, d)] = q
    return od


class Solver:
    def __init__(self, net, od):
        self.init = net["init"]
        self.term = net["term"]
        self.m = len(self.init)
        self.n = max(max(self.init), max(self.term)) + 1
        self.t0 = net["t0"]
        self.p = net["power"]
        self.k = [net["t0"][a] * net["B"][a] / net["cap"][a] ** net["power"][a]
                  for a in range(self.m)]
        self.lid = {(i, j): a for a, (i, j) in enumerate(zip(self.init, self.term))}
        G = csr_matrix((np.arange(1, self.m + 1, dtype=float),
                        (self.init, self.term)), shape=(self.n, self.n))
        self.perm = (G.data.astype(int) - 1)
        self.G = G
        self.ods = sorted(od)
        self.q0 = [od[x] for x in self.ods]
        self.total0 = sum(self.q0)
        self.origins = sorted({o for o, _ in self.ods})
        self.by_origin = {o: [i for i, (oo, _) in enumerate(self.ods) if oo == o]
                          for o in self.origins}
        self.theta = None
        self.reset()

    def reset(self):
        self.paths = [[] for _ in self.ods]
        self.keys = [dict() for _ in self.ods]
        self.h = [[] for _ in self.ods]
        self.theta = None

    def tc(self, a, x):
        return self.t0[a] + self.k[a] * x ** self.p[a]

    def dtc(self, a, x):
        return self.p[a] * self.k[a] * x ** (self.p[a] - 1.0)

    def costs_of(self, v):
        return [self.tc(a, max(v[a], 0.0)) for a in range(self.m)]

    def tree(self, tlist, o):
        self.G.data = np.asarray(tlist)[self.perm]
        dist, pred = dijkstra(self.G, directed=True, indices=o, return_predecessors=True)
        return dist, pred

    def sp_path(self, pred, o, d):
        links = []
        j = d
        while j != o:
            i = int(pred[j])
            links.append(self.lid[(i, j)])
            j = i
        links.reverse()
        return links

    def link_flows(self):
        v = [0.0] * self.m
        for P, H in zip(self.paths, self.h):
            for p, x in zip(P, H):
                if x != 0.0:
                    for a in p:
                        v[a] += x
        return v

    def set_theta(self, theta):
        if self.theta is None:
            t = list(self.t0)
            for o in self.origins:
                _, pred = self.tree(t, o)
                for i in self.by_origin[o]:
                    p = self.sp_path(pred, o, self.ods[i][1])
                    self.paths[i] = [p]
                    self.keys[i] = {tuple(p): 0}
                    self.h[i] = [self.q0[i] * theta]
        else:
            for i in range(len(self.ods)):
                tot = sum(self.h[i])
                target = self.q0[i] * theta
                self.h[i] = [x * target / tot for x in self.h[i]]
        self.theta = theta
        self.v = self.link_flows()
        self.t = self.costs_of(self.v)

    def gap(self):
        self.v = self.link_flows()
        self.t = self.costs_of(self.v)
        self.G.data = np.asarray(self.t)[self.perm]
        dist = dijkstra(self.G, directed=True, indices=self.origins)
        row = {o: r for r, o in enumerate(self.origins)}
        sptt = math.fsum(self.q0[i] * self.theta * dist[row[o], d]
                         for i, (o, d) in enumerate(self.ods))
        tstt = math.fsum(x * c for x, c in zip(self.v, self.t))
        return (tstt - sptt) / tstt, tstt, dist, row

    def _shift(self, lp, ls, hp):
        v, tc, dtc = self.v, self.tc, self.dtc

        def g(D):
            return (sum(tc(a, max(v[a] - D, 0.0)) for a in lp)
                    - sum(tc(a, v[a] + D) for a in ls))

        def dg(D):
            return -(sum(dtc(a, max(v[a] - D, 0.0)) for a in lp)
                     + sum(dtc(a, v[a] + D) for a in ls))

        if g(hp) >= 0.0:
            return hp
        lo, hi = 0.0, hp
        D = 0.0
        gD = g(D)
        for _ in range(200):
            d = dg(D)
            Dn = D - gD / d if d < 0.0 else 0.5 * (lo + hi)
            if not (lo < Dn < hi):
                Dn = 0.5 * (lo + hi)
            D = Dn
            gD = g(D)
            if gD > 0.0:
                lo = D
            else:
                hi = D
            if gD == 0.0 or hi - lo <= 1e-16 * max(hp, 1.0):
                break
        return D

    def equalize(self, i):
        P, H = self.paths[i], self.h[i]
        t, v = self.t, self.v
        for _ in range(INNER_ROUNDS):
            c = [sum(t[a] for a in p) for p in P]
            s = min(range(len(P)), key=lambda j: c[j])
            moved = False
            ss = set(P[s])
            for j in range(len(P)):
                if j == s or H[j] <= 0.0 or c[j] <= c[s]:
                    continue
                sj = set(P[j])
                lp = [a for a in P[j] if a not in ss]
                ls = [a for a in P[s] if a not in sj]
                D = self._shift(lp, ls, H[j])
                if D <= 0.0:
                    continue
                moved = True
                if D >= H[j]:
                    D = H[j]
                    H[j] = 0.0
                else:
                    H[j] -= D
                H[s] += D
                for a in lp:
                    v[a] = max(v[a] - D, 0.0)
                    t[a] = self.tc(a, v[a])
                for a in ls:
                    v[a] += D
                    t[a] = self.tc(a, v[a])
            if not moved:
                break

    def sweep(self):
        for o in self.origins:
            _, pred = self.tree(self.t, o)
            for i in self.by_origin[o]:
                p = self.sp_path(pred, o, self.ods[i][1])
                key = tuple(p)
                if key not in self.keys[i]:
                    self.keys[i][key] = len(self.paths[i])
                    self.paths[i].append(p)
                    self.h[i].append(0.0)
                self.equalize(i)

    def solve(self, theta, tol=GAP_TOL, max_sweeps=MAX_SWEEPS):
        self.set_theta(theta)
        sweeps = 0
        while True:
            gp, tstt, dist, row = self.gap()
            if gp <= tol or sweeps >= max_sweeps:
                break
            self.sweep()
            sweeps += 1
        return {"theta": theta, "gap": gp, "tstt": tstt, "sweeps": sweeps,
                "v": list(self.v), "dist": dist, "row": row,
                "total": self.total0 * theta}

    def margin(self, a, res):
        """min over OD of reduced cost of best route through link a (>=0)."""
        dist, row = res["dist"], res["row"]
        i, j = self.init[a], self.term[a]
        ta = self.t[a]
        best = math.inf
        for (o, d) in self.ods:
            r = row[o]
            if i == d or j == o:
                continue
            di = dist[r, i]
            dj = dist[self.origins.index(j), d] if j in row else math.inf
            val = di + ta + dj - dist[r, d]
            best = min(best, val)
        return best


def fw_solve(S, theta, tol, max_iters=FW_MAX_ITERS):
    """Frank-Wolfe with exact (bisection) line search; cold start from AON."""
    q = np.array(S.q0) * theta
    init = np.array(S.init)
    t0 = np.array(S.t0)
    k = np.array(S.k)
    pw = np.array(S.p)

    def cost(v):
        return t0 + k * np.maximum(v, 0.0) ** pw

    def aon(t):
        y = np.zeros(S.m)
        for o in S.origins:
            _, pred = S.tree(list(t), o)
            for i in S.by_origin[o]:
                for a in S.sp_path(pred, o, S.ods[i][1]):
                    y[a] += q[i]
        return y

    v = aon(t0)
    it = 0
    while True:
        t = cost(v)
        y = aon(t)
        tstt = float(v @ t)
        gp = (tstt - float(y @ t)) / tstt
        if gp <= tol or it >= max_iters:
            break
        dvec = y - v
        lo, hi = 0.0, 1.0
        if dvec @ cost(v + dvec) <= 0:
            al = 1.0
        else:
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                if dvec @ cost(v + mid * dvec) > 0:
                    hi = mid
                else:
                    lo = mid
            al = 0.5 * (lo + hi)
        v = v + al * dvec
        it += 1
    return {"v": v.tolist(), "gap": gp, "iters": it, "tstt": tstt}


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    quick = "--quick" in sys.argv
    net_path, trips_path, out_json, out_pdf = args
    net = read_net(net_path)
    od = read_trips(trips_path)
    assert set(net["power"]) == {4.0}
    S = Solver(net, od)
    T0 = time.time()
    log = lambda *x: print(f"[{time.time() - T0:7.1f}s]", *x, flush=True)

    def used(res):
        return [x > USED_FRAC * res["total"] for x in res["v"]]

    if quick:
        r = S.solve(1.0)
        log("theta=1 gap", r["gap"], "sweeps", r["sweeps"], "TSTT", r["tstt"])
        log("unused links", [a for a, u in enumerate(used(r)) if not u])
        return

    # --- grid scan ---
    grid = []
    for th in THETA_GRID:
        r = S.solve(th)
        u = used(r)
        grid.append({"theta": th, "gap": r["gap"], "sweeps": r["sweeps"],
                     "tstt": r["tstt"], "unused_links": [a + 1 for a, x in enumerate(u) if not x],
                     "min_used_link_flow_frac": min([x for x in r["v"] if x > USED_FRAC * r["total"]]) / r["total"],
                     "v": r["v"]})
        log(f"theta={th:.2f} gap={r['gap']:.2e} sweeps={r['sweeps']} unused={len(grid[-1]['unused_links'])}")
    assert all(g["gap"] <= GAP_TOL for g in grid)

    events = []
    for g1, g2 in zip(grid, grid[1:]):
        s1, s2 = set(g1["unused_links"]), set(g2["unused_links"])
        for a in sorted(s1 - s2):
            events.append({"link": a, "type": "activation", "lo": g1["theta"], "hi": g2["theta"]})
        for a in sorted(s2 - s1):
            events.append({"link": a, "type": "deactivation", "lo": g1["theta"], "hi": g2["theta"]})
    log("events", events)
    acts = sorted([e for e in events if e["type"] == "activation"],
                  key=lambda e: (abs(0.5 * (e["lo"] + e["hi"]) - 1.0), e["link"]))[:N_ACT]

    details = []
    for e in acts:
        a = e["link"] - 1
        lo, hi = e["lo"], e["hi"]
        S.solve(lo)
        while hi - lo > BISECT_WIDTH:
            mid = 0.5 * (lo + hi)
            r = S.solve(mid)
            if r["v"][a] > USED_FRAC * r["total"]:
                hi = mid
            else:
                lo = mid
        th_b = 0.5 * (lo + hi)
        # margin (reduced cost) on inactive side, linear extrapolation
        marg = []
        for dlt in [1e-3, 2e-3]:
            r = S.solve(th_b - dlt)
            marg.append((th_b - dlt, S.margin(a, r)))
        (x1, m1), (x2, m2) = marg
        th_m = x1 + m1 * (x1 - x2) / (m2 - m1)
        # flow on the active side
        flows = []
        for dlt in EXP_OFFSETS:
            r = S.solve(th_b + dlt)
            flows.append({"delta": dlt, "flow": r["v"][a], "frac_of_total": r["v"][a] / r["total"],
                          "gap": r["gap"]})
        # linear extrapolation of flow from the two smallest offsets
        f1, f2 = flows[0]["flow"], flows[1]["flow"]
        d1, d2 = flows[0]["delta"], flows[1]["delta"]
        th_f = th_b + d1 - f1 * (d2 - d1) / (f2 - f1)
        th_star = th_f  # flow extrapolation; agrees with bisection to <1e-8
        slopes = []
        for A, Bq in zip(flows, flows[1:]):
            row_ = {"delta_range_above_bisection": [A["delta"], Bq["delta"]]}
            for ref_name, ref in [("ref_margin", th_m), ("ref_flow_extrap", th_f), ("ref_bisection", th_b)]:
                x1 = A["delta"] + th_b - ref
                x2 = Bq["delta"] + th_b - ref
                row_[ref_name] = ((math.log(Bq["flow"]) - math.log(A["flow"])) / (math.log(x2) - math.log(x1))
                                  if x1 > 0 and x2 > 0 else None)
            slopes.append(row_)
        # flow at theta*+0.01, +0.05
        frac = {}
        exact = {}
        for dlt in [0.01, 0.05]:
            r = S.solve(th_star + dlt)
            exact[dlt] = r
            frac[f"{dlt:.2f}"] = {"flow": r["v"][a], "frac_of_total": r["v"][a] / r["total"]}
        # loose-tolerance solves
        loose = {}
        for dlt in [0.01, 0.05]:
            ex = np.array(exact[dlt]["v"])
            tot = exact[dlt]["total"]
            truly_unused = ex <= USED_FRAC * tot
            for tol in LOOSE_TOLS:
                S2 = Solver(net, od)
                rp = S2.solve(th_star + dlt, tol=tol)
                rf = fw_solve(S, th_star + dlt, tol)
                for name, rr, extra in [("path_based_cold", rp, {"sweeps": rp["sweeps"]}),
                                        ("frank_wolfe", rf, {"iters": rf["iters"]})]:
                    vv = np.array(rr["v"])
                    err = np.abs(vv - ex)
                    loose[f"d{dlt:.2f}_tol{tol:.0e}_{name}"] = dict(
                        extra, gap=rr["gap"], new_link_flow=float(vv[a]),
                        new_link_reported_used=bool(vv[a] > USED_FRAC * tot),
                        new_link_abs_error=float(err[a]),
                        max_abs_link_error=float(err.max()),
                        max_abs_link_error_frac_total=float(err.max() / tot),
                        false_used_links=[int(b) + 1 for b in np.where(truly_unused & (vv > USED_FRAC * tot))[0]],
                        false_unused_links=[int(b) + 1 for b in np.where((~truly_unused) & (vv <= USED_FRAC * tot))[0]])
            log("loose done", e["link"], dlt)
        # hidden-activation window: loose solvers close to theta*
        hidden = []
        for dlt in HIDDEN_OFFSETS:
            th = th_star + dlt
            rex = S.solve(th)
            ex = rex["v"][a]
            tot = rex["total"]
            row_ = {"delta": dlt, "exact_flow": ex, "exact_frac_of_total": ex / tot,
                    "exact_used": bool(ex > USED_FRAC * tot)}
            for tol in LOOSE_TOLS:
                S2 = Solver(net, od)
                rp = S2.solve(th, tol=tol)
                rf = fw_solve(S, th, tol)
                exv = np.array(rex["v"])
                for name, rr in [("path_based_cold", rp), ("frank_wolfe", rf)]:
                    vv = np.array(rr["v"])
                    row_[f"tol{tol:.0e}_{name}"] = {
                        "gap": rr["gap"], "reported_flow": float(vv[a]),
                        "reported_used": bool(vv[a] > USED_FRAC * tot),
                        "max_abs_link_error_all_links": float(np.abs(vv - exv).max())}
            hidden.append(row_)
        # TSTT one-sided slopes
        def tstt(th):
            return S.solve(th)["tstt"]
        fd = {}
        for center, label in [(th_star, "at_theta_star"), (th_star + FD_CONTROL_OFFSET, "control_theta_star_plus_5e-4")]:
            f0 = tstt(center)
            for hstep in FD_STEPS:
                fm1, fm2 = tstt(center - hstep), tstt(center - 2 * hstep)
                fp1, fp2 = tstt(center + hstep), tstt(center + 2 * hstep)
                sl = (3 * f0 - 4 * fm1 + fm2) / (2 * hstep)
                sr = (-3 * f0 + 4 * fp1 - fp2) / (2 * hstep)
                fd[f"{label}_h{hstep:.0e}"] = {"left_slope": sl, "right_slope": sr,
                                               "relative_jump": (sr - sl) / abs(sl)}
        # local curve for figure
        curve = []
        for i in range(-20, 61):
            th = th_star + 0.001 * i
            r = S.solve(th)
            curve.append([th, r["v"][a]])
        details.append({"link": e["link"], "link_nodes": [S.init[a] + 1, S.term[a] + 1],
                        "grid_bracket": [e["lo"], e["hi"]],
                        "theta_star_bisection_1e-9_threshold": th_b,
                        "theta_star_margin_extrapolation": th_m,
                        "theta_star_flow_extrapolation": th_f,
                        "theta_star_used": th_star,
                        "margin_samples": [list(m) for m in marg],
                        "flows_above_threshold": flows, "local_log_log_slopes": slopes,
                        "flow_at_offsets": frac, "loose_tolerance": loose,
                        "tstt_one_sided_slopes": fd, "hidden_window": hidden, "curve": curve})
        log("activation done", e["link"], th_star)

    out = {
        "label": "EXPLORATORY illustration; static deterministic UE on the TNTP Sioux Falls test network",
        "data": {
            "net_url": "https://raw.githubusercontent.com/bstabler/TransportationNetworks/master/SiouxFalls/SiouxFalls_net.tntp",
            "trips_url": "https://raw.githubusercontent.com/bstabler/TransportationNetworks/master/SiouxFalls/SiouxFalls_trips.tntp",
            "net_sha256": sha256(net_path), "trips_sha256": sha256(trips_path),
            "n_links": S.m, "n_nodes": S.n, "n_od_pairs_positive": len(S.ods),
            "total_demand_theta1": S.total0},
        "settings": {"cost": "t0*(1+B*(v/C)^4) from file columns free_flow_time, b, capacity, power",
                     "solver": "path-based Gauss-Seidel equalisation (origin-ordered column generation, "
                               "pairwise exact 1-D safeguarded-Newton shifts to shortest path)",
                     "gap_definition": "(TSTT - sum_od D_od*mincost_od)/TSTT",
                     "gap_tol": GAP_TOL, "max_sweeps": MAX_SWEEPS, "inner_rounds": INNER_ROUNDS,
                     "used_threshold_frac_of_total_demand": USED_FRAC,
                     "theta_grid": [THETA_GRID[0], THETA_GRID[-1], 0.01],
                     "bisection_width": BISECT_WIDTH, "fd_steps": FD_STEPS,
                     "loose_tols": LOOSE_TOLS, "fw_max_iters": FW_MAX_ITERS,
                     "fd_control_offset": FD_CONTROL_OFFSET, "hidden_offsets": HIDDEN_OFFSETS,
                     "theta_star_estimator": "linear extrapolation to zero of new-link flow at bisection+1e-5, +3e-5",
                     "warm_start": "path flows rescaled from previous solve"},
        "grid": [{k: g[k] for k in ("theta", "gap", "sweeps", "tstt", "unused_links", "min_used_link_flow_frac")} for g in grid],
        "events": events, "activations": details,
    }
    with open(out_json, "w") as f:
        json.dump(out, f, sort_keys=True, indent=2)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 3.8))
    for dct in details:
        c = np.array(dct["curve"])
        ax.plot(c[:, 0] - dct["theta_star_used"], c[:, 1],
                label=f"link {dct['link']} ({dct['link_nodes'][0]}→{dct['link_nodes'][1]}), θ*={dct['theta_star_used']:.5f}")
    ax.set_xlabel("θ − θ*")
    ax.set_ylabel("flow on newly used link")
    ax.set_title("Sioux Falls (exploratory): link activation along demand multiplier", fontsize=9)
    ax.legend(fontsize=7, loc="upper left")
    ins = ax.inset_axes([0.62, 0.12, 0.34, 0.38])
    for dct in details:
        fl = dct["flows_above_threshold"]
        x = [r["delta"] + dct["theta_star_bisection_1e-9_threshold"] - dct["theta_star_used"] for r in fl]
        ins.loglog(x, [r["flow"] for r in fl], "o-", ms=2)
    ins.set_xlabel("θ − θ*", fontsize=6)
    ins.tick_params(labelsize=5)
    fig.tight_layout()
    fig.savefig(out_pdf, metadata={"CreationDate": None, "ModDate": None})
    log("done")


if __name__ == "__main__":
    main()
