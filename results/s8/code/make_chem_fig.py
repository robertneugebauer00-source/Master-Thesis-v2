## make_chem_fig.py -- "the chemical side of atypical sites" figure (Elbe).
## 4 panels: (a) strength~stress decoupling, (b) chemical fingerprints,
## (c) global edge-weight shift, (d) TU drivers: same toxicant, different context.
import numpy as np, pandas as pd, json, warnings
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.covariance import graphical_lasso
from scipy import stats as sst

OUT = '/mnt/agents/output'
d = np.load(f'{OUT}/s8_backbone/prep.npz', allow_pickle=True)
X = d['Xlog__Basin_Elbe']; sites = list(d['sites__Basin_Elbe']); comps = list(d['comps__Basin_Elbe'])
n, p = X.shape
sum_tu = pd.Series(d['sum_tu'], index=d['sum_tu_sites']).reindex(sites)
lam = json.load(open(f'{OUT}/s8_backbone/run_summary.json'))['Basin_Elbe']['lambda_star']

## BONOBO nets (exact R port)
xbar = X.mean(axis=0)
nets = np.empty((n, p, p))
for q in range(n):
    Xq = np.delete(X, q, axis=0)
    S_q = np.cov(Xq, rowvar=False)
    dg = np.diag(S_q)
    v_pop = np.mean((dg - dg.mean())**2)
    delta = 1/(3 + 2*np.mean(np.sqrt(dg))/v_pop) if v_pop > 0 else 1/3
    dx = X[q] - xbar
    Sig = delta*np.outer(dx, dx) + (1-delta)*S_q
    sd = np.sqrt(np.diag(Sig)); sd[(sd == 0) | ~np.isfinite(sd)] = 1
    R = Sig/np.outer(sd, sd)
    nets[q] = (R + R.T)/2

## aggregate skeleton at lambda*
S_full = np.corrcoef(X, rowvar=False)
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    _, prec = graphical_lasso(emp_cov=S_full, alpha=lam, tol=1e-3, max_iter=500)
dg = np.sqrt(np.diag(prec)); A_agg = -prec/np.outer(dg, dg); np.fill_diagonal(A_agg, 0)
ut = np.triu_indices(p, 1)
skel = np.abs(A_agg[ut]) > 1e-8
ei, ej = ut[0][skel], ut[1][skel]

## detection + floors
mdl_tab = pd.read_csv('/mnt/agents/upload/MDL_CECs_summary.csv')
xl = '/mnt/agents/upload/Finckh_Carmona_2023_PANGAEA_R1.xlsx'
ab = pd.read_excel(xl, sheet_name='A_B_results', header=None)
name_row = ab.index[ab[0] == 'Name'][0]
first_eus = next(i for i, v in enumerate(ab.iloc[name_row]) if isinstance(v, str) and v.startswith('EUS_'))
compounds_all = ab.iloc[name_row+1:, 0].astype(str).tolist()
ik_all = ab.iloc[name_row+1:, 5].astype(str).tolist()
site_ids_all = [str(v) for v in ab.iloc[name_row, first_eus:]]
conc_arr = np.nan_to_num(ab.iloc[name_row+1:, first_eus:].apply(pd.to_numeric, errors='coerce').to_numpy(float))
cm = pd.DataFrame(conc_arr, index=compounds_all, columns=site_ids_all)
mdl_by_ik = dict(zip(mdl_tab['inchikey'].astype(str).str.strip().str.upper(), mdl_tab['mdl_value_ng_L']))
mdl_by_name = dict(zip(mdl_tab['chemical_name'].astype(str).str.strip().str.lower(), mdl_tab['mdl_value_ng_L']))
mv = np.array([mdl_by_ik.get(i.strip().upper(), np.nan) for i in ik_all])
un = np.isnan(mv); mv[un] = [mdl_by_name.get(c.strip().lower(), np.nan) for c in np.array(compounds_all)[un]]
un = np.isnan(mv); mh = mv/2; mh[un] = 1.0
log_floor = np.log10(pd.Series(mh, index=compounds_all).reindex(comps).to_numpy())
det_all = X > log_floor[None, :] + 1e-12

## TU drivers
pred = pd.read_csv('/mnt/agents/upload/TRIDENT_prediction_results_EC10_GRO.csv')
mapf = pd.read_csv('/mnt/agents/output/lioness_bonobo_addition/data-raw/Danube_TRIDENT_mapping_reconstructed.csv')
mg = mapf.merge(pred[['SMILES', 'predictions (mg/L)']], on='SMILES', how='left')
ec10 = pd.Series((mg['predictions (mg/L)']*1e6).to_numpy(), index=mg['cmpdname'])
ec10 = ec10[np.isfinite(ec10)]
have = [c for c in cm.index if c in ec10.index]
tu_df = cm.loc[have].div(ec10.loc[have], axis=0)

strength = np.abs(nets[:, ei, ej]).sum(axis=1)
n_det = det_all.sum(axis=1)
med = np.median(X, axis=0)
rho, pv = sst.spearmanr(strength, sum_tu)

NAMED = ['EUS_434', 'EUS_053', 'EUS_027', 'EUS_244', 'EUS_254', 'EUS_029']
LBL = {'EUS_434': 'EUS_434\natypical (rank 1/153)', 'EUS_053': 'EUS_053\natypical',
       'EUS_027': 'EUS_027\natypical', 'EUS_244': 'EUS_244\ntypical (clean)',
       'EUS_254': 'EUS_254\ntypical', 'EUS_029': 'EUS_029\nextreme but typical'}
TEAL, RED, SAND, SLATE, GREY = '#2A9D8F', '#9B2226', '#E9C46A', '#264653', '0.75'

fig, axes = plt.subplots(2, 2, figsize=(16.5, 13))
fig.subplots_adjust(hspace=0.34, wspace=0.30, top=0.93, bottom=0.07, left=0.09, right=0.97)

## (a) decoupling
ax = axes[0, 0]
sc = ax.scatter(sum_tu, strength, c=n_det, cmap='viridis', s=26, alpha=0.75, edgecolor='none')
ax.set_xscale('log')
OFF = {'EUS_434': (-112, -4), 'EUS_053': (-58, -16), 'EUS_027': (8, 6),
       'EUS_244': (-88, 6), 'EUS_254': (-76, -20), 'EUS_029': (-118, -4)}
for s in NAMED:
    q = sites.index(s)
    ax.annotate(LBL[s], (sum_tu[s], strength[q]), textcoords='offset points',
                xytext=OFF[s], fontsize=7.5,
                color=RED if 'atypical' in LBL[s] else SLATE)
    ax.scatter([sum_tu[s]], [strength[q]], s=60, facecolor='none',
               edgecolor=RED if 'atypical' in LBL[s] else SLATE, linewidth=1.6, zorder=5)
ax.set_xlabel('site stress  sum-TU (log scale)')
ax.set_ylabel('BONOBO strength on skeleton')
ax.set_title(f'(a) the decoupling: stress rises, network participation falls\n'
             f'Elbe, Spearman rho = {rho:.2f}, p = {pv:.1e}   (colour = compounds detected)', fontsize=10.5)
cb = fig.colorbar(sc, ax=ax, pad=0.02); cb.set_label('n detects', fontsize=8)

## (b) fingerprints
ax = axes[0, 1]
def fp(site, kpos=6, kneg=6):
    q = sites.index(site)
    dev = pd.Series(X[q] - med, index=comps)
    return pd.concat([dev.nlargest(kpos), dev.nsmallest(kneg)]), det_all[q]
t434, d434 = fp('EUS_434')
dev244 = pd.Series(X[sites.index('EUS_244')] - med, index=comps)
t244 = dev244.nlargest(6)                      # the typical site's signature is all-positive
d244 = det_all[sites.index('EUS_244')]
ylab = list(dict.fromkeys(list(t434.index) + list(t244.index)))
# order: positives of both sites on top (by max dev), then negatives
maxdev = {c: max(t434.get(c, np.nan) if not np.isnan(t434.get(c, np.nan)) else -99,
                 t244.get(c, np.nan) if not np.isnan(t244.get(c, np.nan)) else -99) for c in ylab}
ylab = sorted(ylab, key=lambda c: -maxdev[c])
y = np.arange(len(ylab))[::-1]
med244 = pd.Series(X[sites.index('EUS_244')] - med, index=comps)
v434 = [t434.get(c, np.nan) for c in ylab]; v244 = [med244.get(c, np.nan) for c in ylab]
c434 = [TEAL if d434[comps.index(c)] else GREY for c in ylab]
c244 = [TEAL if d244[comps.index(c)] else GREY for c in ylab]
ax.barh(y + 0.21, [0 if np.isnan(v) else v for v in v434], height=0.38, color=c434, alpha=0.9, label='EUS_434 (atypical)')
ax.barh(y - 0.21, [0 if np.isnan(v) else v for v in v244], height=0.38, color=c244, alpha=0.45,
        hatch='//', edgecolor='white', label='EUS_244 (typical)')
ax.set_yticks(y); ax.set_yticklabels(ylab, fontsize=8)
ax.axvline(0, color='0.4', lw=0.8)
ax.set_xlabel('log10 concentration vs Elbe median  (deviation)')
ax.set_title('(b) chemical fingerprints: deviation from the basin median\n'
             'solid colour = detected, grey = below MDL here; hatched = typical site', fontsize=10.5)
ax.legend(fontsize=8, loc='lower right')

## (c) global weight shift
ax = axes[1, 0]
data = [nets[sites.index(s), ei, ej] for s in NAMED]
bp = ax.boxplot(data, vert=True, showfliers=False, widths=0.55, patch_artist=True,
                medianprops=dict(color='black', lw=1.4))
for patch, s in zip(bp['boxes'], NAMED):
    patch.set_facecolor(RED if 'atypical' in LBL[s] else SLATE)
    patch.set_alpha(0.55 if 'atypical' in LBL[s] else 0.35)
ax.axhline(np.median(nets[:, ei, ej]), color=SAND, lw=1.6, ls='--', label='basin median weight')
ax.set_xticklabels([s.replace('EUS_', '') for s in NAMED], fontsize=8.5)
ax.set_xlabel('site (EUS_)')
ax.set_ylabel('per-site BONOBO weight over 1478 skeleton edges')
ax.set_title('(c) atypicality is a GLOBAL down-shift of all edge weights,\nnot a few broken links (red = atypical, blue = typical)', fontsize=10.5)
ax.legend(fontsize=8)

## (d) TU drivers
ax = axes[1, 1]
for k, (s, col, off) in enumerate([('EUS_029', SLATE, 0.21), ('EUS_434', RED, -0.21)]):
    top = tu_df[s].sort_values(ascending=False).head(8)
    names = list(dict.fromkeys(list(top.index)))
    if k == 0:
        allnames = list(dict.fromkeys(list(tu_df['EUS_029'].sort_values(ascending=False).head(8).index) +
                                      list(tu_df['EUS_434'].sort_values(ascending=False).head(8).index)))
        yy = np.arange(len(allnames))[::-1]
    vals = [tu_df.loc[c, s] if c in tu_df.index else 0 for c in allnames]
    ax.barh(yy + off, vals, height=0.38, color=col, alpha=0.75 if k == 0 else 0.6,
            label=f'{s} (sum-TU {sum_tu[s]:.0f})')
ax.set_yticks(yy); ax.set_yticklabels(allnames, fontsize=8)
ax.set_xscale('log'); ax.set_xlabel('toxic units TU = MEC / EC10 (log scale)')
ax.set_title('(d) where the stress comes from: ONE dominant toxicant\n'
             '(Methocarbamol) -- but EUS_029 stays typical, EUS_434 is atypical', fontsize=10.5)
ax.legend(fontsize=8, loc='lower right')

fig.suptitle('The chemical side of single-site networks: what makes a site "atypical" (Elbe, BONOBO)',
             fontsize=13, fontweight='bold')
fig.savefig(f'{OUT}/fig_chemical_atypicality_elbe.png', dpi=150)
print('saved', f'{OUT}/fig_chemical_atypicality_elbe.png')
