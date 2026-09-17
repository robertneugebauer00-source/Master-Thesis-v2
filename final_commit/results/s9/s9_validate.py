## s9_validate.py -- Python mirror of analysis/S9_simulation.Rmd
## Purpose: numerical validation of the S9 design (no R in this environment).
## Mirrors the R chunk logic 1:1; calibration from s8_backbone/prep.npz (Elbe).
## MDL floors are NOT in the npz, so they are approximated: a compound whose
## column minimum is a point mass (>10% of values at the min) is treated as
## censored at MDL = 10^min (i.e. floor = min in log10 space); all other
## compounds are treated as uncensored. This is a validation-only proxy --
## the R version uses the real MDL/2 vector from prep_unit()$pseudo.
import numpy as np, json, time
from scipy import stats

rng = np.random.default_rng(42)
t0 = time.time()

## ---- calibration -------------------------------------------------------------
d = np.load('/mnt/agents/output/s8_backbone/prep.npz', allow_pickle=True)
Xcal = d['Xlog__Basin_Elbe']
elbe_sites = [s.decode() if isinstance(s, bytes) else str(s) for s in d['sites__Basin_Elbe']]
tu_sites = [s.decode() if isinstance(s, bytes) else str(s) for s in d['sum_tu_sites']]
tu_all = d['sum_tu']
idx = [tu_sites.index(s) for s in elbe_sites]
TU_POOL = tu_all[idx]

n_cal, p = Xcal.shape
MU = Xcal.mean(0); SDV = Xcal.std(0, ddof=1)
SIGMA = np.corrcoef(Xcal, rowvar=False)

## MDL proxy (see header)
FLOORLOG = Xcal.min(0)                       # log10 floor per compound
frac_at_min = (Xcal == Xcal.min(0)).mean(0)
CENSORED = frac_at_min > 0.10                # point mass at min => censored compound

## PC1 source direction
w, V = np.linalg.eigh(SIGMA)
G_PC = V[:, -1] * np.sign(V[:, -1].sum())

## broken structure: shrink 50% of off-diagonal entries, repair to nearest cor
BROKEN_FRAC = 0.5
ut = np.triu_indices(p, 1)
sel = rng.random(len(ut[0])) < BROKEN_FRAC
SIGMA_B = SIGMA.copy()
SIGMA_B[ut[0][sel], ut[1][sel]] *= 0.1
SIGMA_B[ut[1][sel], ut[0][sel]] *= 0.1
def near_cor(M, eps=1e-4):
    M = (M + M.T) / 2
    w_, V_ = np.linalg.eigh(M)
    w_ = np.maximum(w_, eps)
    M2 = (V_ * w_) @ V_.T
    dd = np.sqrt(np.diag(M2))
    M2 = M2 / np.outer(dd, dd)
    return (M2 + M2.T) / 2
SIGMA_B = near_cor(SIGMA_B)

def factor(Sigma, sdv):
    C = np.outer(sdv, sdv) * Sigma
    w_, V_ = np.linalg.eigh(C)
    return V_ * np.sqrt(np.maximum(w_, 0))
L_INTACT = factor(SIGMA, SDV)
L_BROKEN = factor(SIGMA_B, SDV)

## skeleton: top-5% |r| fallback (R uses the cached glasso skeleton when present)
thr = np.quantile(np.abs(SIGMA[ut]), 0.95)
SKEL = np.abs(SIGMA) > thr
np.fill_diagonal(SKEL, False)

## ---- params (mirror S9) -------------------------------------------------------
BETA_KAPPA, KAPPA_NOISE, ALPHA_MEAN = 1.5, 1.0, 0.5
NTAXA, DEPTH, DM_PHI = 300, 10000, 10
BETA1 = 0.8
GRID_N = [50, 100, 200, 500]
GRID_BETA2 = [0, 0.3, 0.6, 1.0]
NREP = 10

## ---- chemistry simulator -------------------------------------------------------
def sim_chemistry(n, rng):
    s = rng.choice(TU_POOL, n, replace=True)
    zs = (np.log10(s) - np.log10(s).mean()) / np.log10(s).std(ddof=1)
    kappa = 1 / (1 + np.exp(-(BETA_KAPPA * zs + rng.normal(0, KAPPA_NOISE, n))))
    Z1 = rng.normal(size=(n, p)) @ L_INTACT.T
    Z2 = rng.normal(size=(n, p)) @ L_BROKEN.T
    X = Z1 * np.sqrt(1 - kappa)[:, None] + Z2 * np.sqrt(kappa)[:, None]
    X = X + MU + ALPHA_MEAN * np.outer(zs, G_PC)
    ## censoring: values below the floor (log10 MDL/2 proxy) are set to the floor
    X[:, CENSORED] = np.maximum(X[:, CENSORED], FLOORLOG[CENSORED])
    return X, kappa, s

## ---- estimators (validated ports, vectorised LOO downdates) --------------------
def loo_cors(X):
    n = X.shape[0]
    S1 = X.sum(0); S2 = (X**2).sum(0); Sxy = X.T @ X
    n1 = n - 1
    out = np.empty((n, p, p))
    for q in range(n):
        x = X[q]
        s1 = S1 - x; s2 = S2 - x**2; sxy = Sxy - np.outer(x, x)
        num = sxy - np.outer(s1, s1) / n1
        den = s2 - s1**2 / n1
        den[den <= 0] = np.nan
        r = num / np.sqrt(np.outer(den, den))
        r[np.isnan(r)] = 0
        out[q] = (r + r.T) / 2
    return out

def lioness_strength(X):
    n = X.shape[0]
    G_all = np.corrcoef(X, rowvar=False)
    loos = loo_cors(X)
    e = n * G_all[None] - (n - 1) * loos
    return np.abs(e[:, SKEL]).sum(1)

def bonobo_strength(X):
    n = X.shape[0]
    S1 = X.sum(0); S2 = (X**2).sum(0); Sxy = X.T @ X
    xbar = X.mean(0)
    n1 = n - 1
    out = np.empty(n)
    for q in range(n):
        x = X[q]
        s1 = S1 - x; sxy = Sxy - np.outer(x, x)
        S_q = (sxy - np.outer(s1, s1) / n1) / (n1 - 1)     # LOO covariance, ddof=1
        dg = np.diag(S_q)
        v_pop = ((dg - dg.mean())**2).mean()               # population variance
        delta = 1 / (3 + 2 * np.sqrt(dg).mean() / v_pop) if v_pop > 0 else 1/3
        dx = x - xbar
        Sq = delta * np.outer(dx, dx) + (1 - delta) * S_q
        sd = np.sqrt(np.diag(Sq)); sd[sd == 0] = 1
        R = Sq / np.outer(sd, sd)
        out[q] = np.abs(R[SKEL]).sum()
    return out

## ---- microbiome simulator -------------------------------------------------------
def sim_microbiome(kappa, tu, beta1, beta2, rng):
    n = len(kappa)
    lt = np.log10(tu); zt = (lt - lt.mean()) / lt.std(ddof=1)
    zk = (kappa - kappa.mean()) / kappa.std(ddof=1)
    n_sens = round(0.2 * NTAXA); n_tol = round(0.1 * NTAXA)
    sens = np.array([-1]*n_sens + [1]*n_tol + [0]*(NTAXA - n_sens - n_tol))
    base = rng.lognormal(0, 1.2, NTAXA)
    eta = np.outer(zt, beta1 * sens) + np.outer(zk, beta2 * sens) + np.log(base)
    alpha = np.exp(eta) * DM_PHI
    Gm = rng.gamma(alpha, 1)
    P = Gm / Gm.sum(1, keepdims=True)
    counts = np.array([rng.multinomial(DEPTH, P[q]) for q in range(n)])
    return counts[:, sens == -1].sum(1) / DEPTH

## ---- biomarker F-test (nested OLS, mirrors anova(lm,lm)) ------------------------
def biomarker_p(metric, tu, strength):
    n = len(metric)
    lt = np.log10(tu); zt = (lt - lt.mean()) / lt.std(ddof=1)
    zs = (strength - strength.mean()) / strength.std(ddof=1)
    X0 = np.column_stack([np.ones(n), zt])
    X1 = np.column_stack([np.ones(n), zt, zs])
    b0 = np.linalg.lstsq(X0, metric, rcond=None)[0]
    b1 = np.linalg.lstsq(X1, metric, rcond=None)[0]
    rss0 = ((metric - X0 @ b0)**2).sum(); rss1 = ((metric - X1 @ b1)**2).sum()
    F = (rss0 - rss1) / (rss1 / (n - 3))
    return stats.f.sf(F, 1, n - 3)

def auc(score, truth):
    r = stats.rankdata(score)
    m = truth.sum(); nn = len(truth) - m
    return (r[truth == 1].sum() - m * (m + 1) / 2) / (m * nn)

## ---- recovery grid ----------------------------------------------------------------
recovery = []
for n in GRID_N:
    for rep in range(NREP):
        r_ = np.random.default_rng(42 + 100 * n + rep)
        X, kappa, tu = sim_chemistry(n, r_)
        sb = bonobo_strength(X); sl = lioness_strength(X)
        atyp = (kappa > np.quantile(kappa, 0.75)).astype(int)
        recovery.append(dict(n=n, rep=rep,
            rho_bonobo=stats.spearmanr(sb, kappa).statistic,
            rho_lioness=stats.spearmanr(sl, kappa).statistic,
            auc_bonobo=auc(-sb, atyp), auc_lioness=auc(-sl, atyp)))
    print(f"[recovery] n={n} done  ({time.time()-t0:.0f}s)", flush=True)

## ---- power grid (BONOBO only, beta1 = 0) -------------------------------------------
power = []
for i, n in enumerate(GRID_N):
    for b2 in GRID_BETA2:
        pv = []
        for rep in range(NREP):
            r_ = np.random.default_rng(42 + 100000 * (i * 4 + GRID_BETA2.index(b2)) + rep)
            X, kappa, tu = sim_chemistry(n, r_)
            sb = bonobo_strength(X)
            metric = sim_microbiome(kappa, tu, 0, b2, r_)
            pv.append(biomarker_p(metric, tu, sb))
        power.append(dict(n=n, beta2=b2, power=float(np.mean(np.array(pv) < 0.05)),
                          median_p=float(np.median(pv))))
        print(f"[power] n={n} beta2={b2} power={np.mean(np.array(pv)<0.05):.2f}  ({time.time()-t0:.0f}s)", flush=True)

## ---- controls ------------------------------------------------------------------------
pv_amount, pv_both = [], []
for rep in range(NREP):
    r_ = np.random.default_rng(42 + 900000 + rep)
    X, kappa, tu = sim_chemistry(200, r_)
    sb = bonobo_strength(X)
    pv_amount.append(biomarker_p(sim_microbiome(kappa, tu, BETA1, 0, r_), tu, sb))
    r_ = np.random.default_rng(42 + 950000 + rep)
    X, kappa, tu = sim_chemistry(200, r_)
    sb = bonobo_strength(X)
    pv_both.append(biomarker_p(sim_microbiome(kappa, tu, BETA1, 0.6, r_), tu, sb))
controls = dict(H_amount_rej=float(np.mean(np.array(pv_amount) < 0.05)),
                H_both_rej=float(np.mean(np.array(pv_both) < 0.05)))

res = dict(recovery=recovery, power=power, controls=controls,
           note="MDL floors approximated from point masses (validation-only proxy)",
           runtime_s=time.time() - t0)
with open('/mnt/agents/output/s9_sim/s9_validation_results.json', 'w') as f:
    json.dump(res, f, indent=1)
print(f"DONE in {time.time()-t0:.0f}s")
