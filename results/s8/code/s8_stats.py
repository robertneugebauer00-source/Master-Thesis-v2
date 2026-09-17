## s8_stats.py -- S8 robustness statistics from the per-unit lioness_glasso run.
## Mirrors analysis/S8_backbone_sensitivity.Rmd: pooled permutation tests, UDF
## external proxy, trivial-score benchmark, biomarker F-test, S7 comparison,
## 2-panel figure. Covariates (UDF, basin, n_det, mean_log) come from the S7
## site_table so definitions match the delivered S7 numbers exactly.
import numpy as np, pandas as pd, json, glob
from scipy import stats as sst
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = '/mnt/agents/output/s8_backbone'
NPERM, SEED, MIN_UDF = 5000, 42, 10
rng = np.random.default_rng(SEED)

rows = []
for f in sorted(glob.glob(f'{OUT}/per_unit/*.csv')):
    rows.append(pd.read_csv(f))
s8 = pd.concat(rows, ignore_index=True)
st = pd.read_csv('/mnt/agents/output/s7_site_level/site_table.csv')
s8 = s8.merge(st[['unit', 'site', 'UDF', 'stream_order', 'basin', 'n_det', 'mean_log']],
              on=['unit', 'site'], how='left')
print('S8 site table:', s8.shape, '| units:', s8.unit.nunique())
s8.to_csv(f'{OUT}/tableS8_site_table.csv', index=False)

units = sorted(s8.unit.unique())
def fast_rho(ry, x):
    return np.corrcoef(ry, pd.Series(x).rank())[0, 1]
def stouffer(rhos, n):
    z = np.arctanh(np.clip(rhos, -.999, .999)) * np.sqrt(np.array(n) - 3)
    return z.sum() / np.sqrt(len(rhos))

ranks_lg = [s8.loc[s8.unit == u, 'strength_lionessglasso'].rank().to_numpy() for u in units]
ranks_bo = [s8.loc[s8.unit == u, 'strength_bonobo_pyskel'].rank().to_numpy() for u in units]
x_vals   = [s8.loc[s8.unit == u, 'sum_tu'].to_numpy() for u in units]
ns       = [len(x) for x in x_vals]

obs_lg = np.array([fast_rho(r, x) for r, x in zip(ranks_lg, x_vals)])
obs_bo = np.array([fast_rho(r, x) for r, x in zip(ranks_bo, x_vals)])
null = np.empty((len(units), NPERM))
for b in range(NPERM):
    null[:, b] = [fast_rho(r, rng.permutation(x)) for r, x in zip(ranks_lg, x_vals)]
p_emp = (np.abs(null) >= np.abs(obs_lg[:, None])).mean(axis=1)

obs = dict(neg=int((obs_lg < 0).sum()), med=float(np.median(obs_lg)), Z=float(stouffer(obs_lg, ns)))
nstats = dict(neg=(null < 0).sum(axis=0), med=np.median(null, axis=0),
              Z=np.array([stouffer(null[:, b], ns) for b in range(NPERM)]))
p_neg = float((nstats['neg'] >= obs['neg']).mean())
p_med = float((nstats['med'] <= obs['med']).mean())
p_Z   = float((nstats['Z'] <= obs['Z']).mean())
Z_bo  = float(stouffer(obs_bo, ns))

per_unit = pd.DataFrame({'unit': units, 'n_sites': ns,
                         'rho_lionessglasso': np.round(obs_lg, 3),
                         'p_emp_lionessglasso': np.round(p_emp, 4),
                         'rho_bonobo_same_skeleton': np.round(obs_bo, 3)})
per_unit.to_csv(f'{OUT}/tableS8_per_unit_permutation.csv', index=False)

gt = pd.DataFrame({
    'test': ['sign count negative (lioness_glasso)', 'median rho (lioness_glasso)',
             'Stouffer Z (lioness_glasso)', 'Stouffer Z (lioness_glasso)',
             'Stouffer Z (bonobo, same skeleton)'],
    'statistic': np.round([obs['neg'], obs['med'], obs['Z'], obs['Z'], Z_bo], 4),
    'p_value': [p_neg, p_med, 2*sst.norm.cdf(-abs(obs['Z'])), p_Z, 2*sst.norm.cdf(-abs(Z_bo))],
    'method': ['permutation', 'permutation', 'normal theory', 'permutation', 'normal theory']})
gt.to_csv(f'{OUT}/tableS8_global_tests.csv', index=False)
print(gt.to_string(index=False))

## UDF external proxy (basin units)
B = s8[s8.unit.str.startswith('Basin_')].copy()
U = B.dropna(subset=['UDF'])
udf_rows = []
for bn, d in U.groupby('basin'):
    if len(d) < MIN_UDF: continue
    r_lg, p_lg = sst.spearmanr(d.strength_lionessglasso, d.UDF)
    r_bo, p_bo = sst.spearmanr(d.strength_bonobo_pyskel, d.UDF)
    udf_rows.append(dict(basin=bn, n=len(d),
                         rho_lionessglasso_strength_UDF=round(r_lg, 3), p_lionessglasso=round(p_lg, 4),
                         rho_bonobo_strength_UDF=round(r_bo, 3), p_bonobo=round(p_bo, 4)))
udf_tab = pd.DataFrame(udf_rows)
udf_tab.to_csv(f'{OUT}/tableS8_external_proxy_UDF.csv', index=False)
print(udf_tab.to_string(index=False))

## pooled FE OLS + Huber (lioness_glasso)
import statsmodels.api as sm
B['z_strength'] = B.groupby('basin')['strength_lionessglasso'].transform(lambda s: (s - s.mean()) / s.std())
B['log_tu'] = np.log1p(B['sum_tu'])
X = pd.get_dummies(B[['log_tu', 'basin']], drop_first=True).astype(float)
X = sm.add_constant(X)
ols = sm.OLS(B['z_strength'], X).fit(cov_type='cluster', cov_kwds={'groups': B['basin']})
hub = sm.RLM(B['z_strength'], X, M=sm.robust.norms.HuberT()).fit()
print(f"FE OLS slope {ols.params['log_tu']:+.4f} (p {ols.pvalues['log_tu']:.3f}) | "
      f"Huber {hub.params['log_tu']:+.4f} (p {hub.pvalues['log_tu']:.3f})")

## biomarker F-test
Ub = U[U.basin.isin(udf_tab.loc[udf_tab.n >= 30, 'basin'])].copy()
Ub['z_udf'] = Ub.groupby('basin')['UDF'].transform(lambda s: (s - s.mean()) / s.std())
Ub['z_strength'] = Ub.groupby('basin')['strength_lionessglasso'].transform(lambda s: (s - s.mean()) / s.std())
Ub['log_tu'] = np.log1p(Ub['sum_tu'])
X0 = sm.add_constant(pd.get_dummies(Ub[['log_tu', 'basin']], drop_first=True).astype(float))
X1 = X0.copy(); X1['z_strength'] = Ub['z_strength']
m0, m1 = sm.OLS(Ub['z_udf'], X0).fit(), sm.OLS(Ub['z_udf'], X1).fit()
F = ((m0.ssr - m1.ssr) / 1) / m1.mse_resid
pF = 1 - sst.f.cdf(F, 1, m1.df_resid)
print(f'biomarker F = {F:.2f}, p = {pF:.3f}')

## trivial-score benchmark (partial Spearman per basin)
def partial_spearman(d, y, x, covs):
    def rres(v):
        Xm = d[covs].rank().to_numpy()
        Xm = sm.add_constant(Xm)
        return sm.OLS(pd.Series(v).rank().to_numpy(), Xm).fit().resid
    return sst.spearmanr(rres(d[y]), rres(d[x]))[0]
bench = []
for bn, d in B.groupby('basin'):
    bench.append(dict(basin=bn, n=len(d),
                      raw_rho=round(sst.spearmanr(d.strength_lionessglasso, d.sum_tu)[0], 3),
                      partial_rho_given_trivial=round(
                          partial_spearman(d, 'strength_lionessglasso', 'sum_tu', ['n_det', 'mean_log']), 3)))
bench = pd.DataFrame(bench)
bench.to_csv(f'{OUT}/tableS8_benchmark.csv', index=False)
print(bench.to_string(index=False))

## comparison vs S7
s7g = pd.read_csv('/mnt/agents/output/s7_site_level/tableS7_global_tests.csv')
cmp = pd.concat([pd.concat([pd.DataFrame({'backbone': ['correlation (S7)'] * len(s7g)}), s7g], axis=1),
                 pd.concat([pd.DataFrame({'backbone': ['glasso (S8)'] * len(gt)}), gt], axis=1)])
cmp.to_csv(f'{OUT}/tableS8_vs_S7_comparison.csv', index=False)
print(cmp.to_string(index=False))

## figure
fig, axes = plt.subplots(1, 2, figsize=(14.5, 7.2))
ax = axes[0]
ord_ = np.argsort(obs_lg)
s7u = pd.read_csv('/mnt/agents/output/s7_site_level/tableS7_per_unit_permutation.csv')
rho_lio_s7 = s7u.set_index('unit').reindex(units)['rho_lioness'].to_numpy()
y = np.arange(len(units))
ax.axvline(0, color='0.6', lw=1)
ax.scatter(obs_lg[ord_], y + 0.18, marker='^', color='#9B2226', label='lioness_glasso (S8, glasso backbone)', zorder=3)
ax.scatter(obs_bo[ord_], y - 0.18, marker='o', color='#2A9D8F', label='BONOBO (same skeleton)', zorder=3)
ax.scatter(rho_lio_s7[ord_], y, marker='s', color='#E9C46A', label='LIONESS correlation (S7)', zorder=3)
ax.set_yticks(y); ax.set_yticklabels([units[i] for i in ord_], fontsize=7)
ax.set_xlabel('Spearman rho: per-site strength vs site sum-TU')
ax.set_title('(a) per-unit rho: three estimators, one question')
ax.legend(fontsize=8, loc='lower right')
ax = axes[1]
ax.scatter(s8.strength_bonobo_pyskel, s8.strength_lionessglasso, s=8, alpha=0.35, color='#264653')
r = sst.spearmanr(s8.strength_bonobo_pyskel, s8.strength_lionessglasso)[0]
ax.set_xlabel('BONOBO strength (correlation backbone)')
ax.set_ylabel('lioness_glasso strength (glasso backbone)')
ax.set_title(f'(b) per-site agreement (Spearman {r:.2f})')
lims = [0, max(s8.strength_bonobo_pyskel.max(), s8.strength_lionessglasso.max())]
ax.plot(lims, lims, ls=':', color='0.6')
fig.suptitle('S8 -- backbone sensitivity of the site-level tests', fontweight='bold')
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(f'{OUT}/figS8_backbone_sensitivity.png', dpi=150)
print('figure saved')

json.dump(dict(obs=obs, p_neg=p_neg, p_med=p_med, p_Z=p_Z, Z_bo=Z_bo,
               biomarker_F=float(F), biomarker_p=float(pF),
               fe_ols_slope=float(ols.params['log_tu']), fe_ols_p=float(ols.pvalues['log_tu']),
               huber_slope=float(hub.params['log_tu']), huber_p=float(hub.pvalues['log_tu'])),
          open(f'{OUT}/s8_stats_summary.json', 'w'), indent=1)
print('DONE')
