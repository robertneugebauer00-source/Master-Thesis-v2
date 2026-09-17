## s8_run.py -- S8 backbone-sensitivity run (Python port of lioness_glasso).
## Per unit (conc, 17 units):
##   1. StARS lambda selection (rep.num = 20, nlambda = 30, huge 2.0.1 rules)
##   2. aggregate glasso at lambda* -> pcor A_agg, skeleton
##   3. lioness_glasso: LOO glasso per site, e_q = n*A_agg - (n-1)*A_loo,
##      restricted to the aggregate skeleton -> strength per site
##   4. BONOBO (exact R port) on the SAME python skeleton -> strength per site
##      (same-skeleton comparison + validation vs S7 mean_weight)
## Saves one CSV per unit + a run summary JSON, progressively.
import numpy as np, pandas as pd, json, time, os, sys
from sklearn.covariance import graphical_lasso
from joblib import Parallel, delayed

OUT = '/mnt/agents/output/s8_backbone'
os.makedirs(f'{OUT}/per_unit', exist_ok=True)
logf = open(f'{OUT}/s8_run.log', 'a', buffering=1)
def log(msg):
    line = f'[{time.strftime("%H:%M:%S")}] {msg}'
    print(line); logf.write(line + '\n')

d = np.load(f'{OUT}/prep.npz', allow_pickle=True)
units = [k.split('__', 1)[1] for k in d.files if k.startswith('Xlog__')]
sum_tu = pd.Series(d['sum_tu'], index=d['sum_tu_sites'])
st = pd.read_csv('/mnt/agents/output/s7_site_level/site_table.csv')

NLAMBDA, REPNUM, THRESH = 30, 20, 0.05
TOL, MAXITER = 1e-3, 100

def cor_scale(X):
    S = np.corrcoef(X, rowvar=False)
    return S

def glasso_admm(S, lam, rho=1.0, max_iter=1000, abstol=1e-5, reltol=1e-4):
    ## ADMM glasso (Boyd et al. 2011 sec 6.5) -- fallback for the LOO fits
    ## where sklearn's coordinate descent diverges (near-perfectly correlated
    ## compound pairs at p > n). Validated against sklearn where CD converges:
    ## 99.98% edge agreement, weight r = 0.9999996.
    p = S.shape[0]
    Theta = np.eye(p); Z = np.eye(p); U = np.zeros((p, p))
    for it in range(max_iter):
        M = rho * (Z - U) - S
        M = (M + M.T) / 2
        w, Q = np.linalg.eigh(M)
        eig = (w + np.sqrt(w**2 + 4 * rho)) / (2 * rho)
        Theta = (Q * eig) @ Q.T
        Z_old = Z.copy()
        A = Theta + U
        Z = np.sign(A) * np.maximum(np.abs(A) - lam / rho, 0)
        np.fill_diagonal(Z, np.diag(A))
        U = U + Theta - Z
        r = np.linalg.norm(Theta - Z, 'fro')
        s = rho * np.linalg.norm(Z - Z_old, 'fro')
        eps_p = p * abstol + reltol * max(np.linalg.norm(Theta, 'fro'), np.linalg.norm(Z, 'fro'))
        eps_d = p * abstol + reltol * rho * np.linalg.norm(U, 'fro')
        if r < eps_p and s < eps_d:
            break
    return Theta

def glasso_pcor(S, lam):
    ## retry ladder: plain -> more iterations -> looser tol -> tiny diagonal
    ## jitter (p > n-1 makes S singular) -> ADMM fallback (documented in the
    ## S8 methods note)
    p = S.shape[0]
    for tol, mi, jit in ((TOL, MAXITER, 0.0), (TOL, 500, 0.0), (1e-2, 500, 0.0),
                         (TOL, 500, 1e-6), (1e-2, 1000, 1e-5)):
        try:
            Sj = S + jit * np.eye(p) if jit else S
            _, prec = graphical_lasso(emp_cov=Sj, alpha=lam, tol=tol, max_iter=mi)
            dg = np.sqrt(np.diag(prec))
            pc = -prec / np.outer(dg, dg)
            np.fill_diagonal(pc, 0)
            return (pc + pc.T) / 2
        except Exception:
            continue
    try:
        prec = glasso_admm(S, lam)
        dg = np.sqrt(np.diag(prec))
        pc = -prec / np.outer(dg, dg)
        np.fill_diagonal(pc, 0)
        return (pc + pc.T) / 2
    except Exception:
        return None

def stars_select(X, seed):
    n, p = X.shape
    S_full = cor_scale(X)
    np.fill_diagonal(S_full, 0)
    lmax = np.abs(S_full).max()
    grid = lmax * 0.1 ** (np.arange(NLAMBDA) / (NLAMBDA - 1))
    b = int(np.floor(10 * np.sqrt(n))) if n > 144 else int(np.floor(0.8 * n))
    rng = np.random.default_rng(seed)
    idxs = [rng.choice(n, size=b, replace=False) for _ in range(REPNUM)]
    def one(args):
        r, k = args
        Sb = cor_scale(X[idxs[r]])
        pc = glasso_pcor(Sb, grid[k])
        if pc is None: return k, None
        return k, (np.abs(pc) > 1e-8).astype(float)
    res = Parallel(n_jobs=2)(delayed(one)((r, k)) for r in range(REPNUM) for k in range(NLAMBDA))
    adj = np.zeros((NLAMBDA, p, p))
    cnt = np.zeros(NLAMBDA)
    for k, a in res:
        if a is not None: adj[k] += a; cnt[k] += 1
    ut = np.triu_indices(p, 1)
    var = np.full(NLAMBDA, np.nan)
    for k in range(NLAMBDA):
        if cnt[k] == 0: continue
        xi = adj[k][ut] / cnt[k]
        var[k] = 4 * (xi.sum() - (xi ** 2).sum()) / (p * (p - 1))
    cross = np.where(np.nan_to_num(var, nan=-1) >= THRESH)[0]
    opt = int(max(cross[0] - 1, 0)) if len(cross) else NLAMBDA - 1
    return grid, var, opt

def bonobo_nets(X):
    n, p = X.shape
    xbar = X.mean(axis=0)
    nets = np.empty((n, p, p))
    for q in range(n):
        Xq = np.delete(X, q, axis=0)
        S_q = np.cov(Xq, rowvar=False)          # ddof=1, as R cov()
        dg = np.diag(S_q)
        v_pop = np.mean((dg - dg.mean()) ** 2)  # population variance, ddof=0
        delta = 1 / (3 + 2 * np.mean(np.sqrt(dg)) / v_pop) if v_pop > 0 else 1/3
        dx = X[q] - xbar
        Sig = delta * np.outer(dx, dx) + (1 - delta) * S_q
        sd = np.sqrt(np.diag(Sig)); sd[(sd == 0) | ~np.isfinite(sd)] = 1
        R = Sig / np.outer(sd, sd)
        nets[q] = (R + R.T) / 2
    return nets

sizes = {u: d[f'Xlog__{u}'].shape[0] for u in units}
order = sorted(units, key=lambda u: sizes[u])
log(f'units (asc n): {[(u, sizes[u]) for u in order]}')

summary = {}
if os.path.exists(f'{OUT}/run_summary.json'):     # resume: keep earlier entries
    summary = json.load(open(f'{OUT}/run_summary.json'))
for ui, u in enumerate(order):
    if os.path.exists(f'{OUT}/per_unit/{u}.csv'):   # resume mode: skip finished units
        log(f'{u}: already done, skipping')
        continue
    t0 = time.time()
    X = d[f'Xlog__{u}']; sites = list(d[f'sites__{u}'])
    n, p = X.shape
    try:
        grid, var, opt = stars_select(X, seed=42 + ui)
        lam = float(grid[opt])
        S_full = cor_scale(X)
        A_agg = glasso_pcor(S_full, lam)
        skel_full = np.abs(A_agg) > 1e-8
        ut = np.triu_indices(p, 1)
        skel_ut = skel_full[ut]
        n_edges = int(skel_ut.sum())

        ## lioness_glasso: LOO refits at fixed lambda
        def loo(q):
            Xq = np.delete(X, q, axis=0)
            return glasso_pcor(cor_scale(Xq), lam)
        A_loo = Parallel(n_jobs=2)(delayed(loo)(q) for q in range(n))
        if any(a is None for a in A_loo):
            raise RuntimeError(f'{sum(a is None for a in A_loo)} LOO glasso fits did not converge')
        e = n * A_agg[None] - (n - 1) * np.stack(A_loo)
        e = e * skel_full[None, :, :]          # restrict to aggregate skeleton
        for q in range(n): np.fill_diagonal(e[q], 0)
        strength_lg = np.abs(e[:, ut[0], ut[1]][:, skel_ut]).sum(axis=1)

        ## BONOBO on the same skeleton (validation + same-skeleton comparison)
        nets = bonobo_nets(X)
        strength_bo = np.abs(nets[:, ut[0], ut[1]][:, skel_ut]).sum(axis=1)
        mean_w = nets[:, ut[0], ut[1]].mean(axis=1)

        out = pd.DataFrame({'unit': u, 'site': sites,
                            'strength_lionessglasso': strength_lg,
                            'strength_bonobo_pyskel': strength_bo,
                            'mean_weight_bonobo': mean_w,
                            'sum_tu': sum_tu.reindex(sites).to_numpy()})
        out.to_csv(f'{OUT}/per_unit/{u}.csv', index=False)
        s7 = st[st.unit == u].set_index('site')['mean_weight']
        mw_diff = float(np.max(np.abs(pd.Series(mean_w, index=sites).reindex(s7.index) - s7)))
        summary[u] = dict(n=n, p=p, lambda_star=lam, opt_idx=opt, n_edges=n_edges,
                          bonobo_mean_weight_maxdiff_vs_s7=mw_diff,
                          variability=[None if np.isnan(v) else float(v) for v in var],
                          secs=round(time.time() - t0, 1))
        with open(f'{OUT}/run_summary.json', 'w') as f: json.dump(summary, f, indent=1)
        log(f'{u}: n={n} p={p} lam={lam:.4f} (idx {opt}) edges={n_edges} '
            f'mean_weight_maxdiff={mw_diff:.2e} [{time.time()-t0:.0f}s]')
    except Exception as ex:
        log(f'{u}: FAILED -- {ex!r}')
        summary[u] = dict(error=repr(ex))
        with open(f'{OUT}/run_summary.json', 'w') as f: json.dump(summary, f, indent=1)
log('ALL DONE')
