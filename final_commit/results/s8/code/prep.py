## prep.py -- rebuild the SGH pipeline data prep in Python (R env lost).
## Mirrors R/01_data.R + prep_unit() from R/02_functions.R exactly, kind = "conc".
## Validates against /mnt/agents/output/s7_site_level/site_table.csv and saves
## prep.npz for the StARS + lioness_glasso background job.
import pandas as pd, numpy as np, json

XL   = '/mnt/agents/upload/Finckh_Carmona_2023_PANGAEA_R1.xlsx'
MDL  = '/mnt/agents/upload/MDL_CECs_summary.csv'
PRED = '/mnt/agents/upload/TRIDENT_prediction_results_EC10_GRO.csv'
MAP  = '/mnt/agents/output/lioness_bonobo_addition/data-raw/Danube_TRIDENT_mapping_reconstructed.csv'
OUT  = '/mnt/agents/output/s8_backbone'

PREV_THRESH, BASIN_MIN_N, SECTION_MIN_N = 0.10, 20, 30
MIN_DETECT_UNIT, DEFAULT_PSEUDO_CONC = 2, 1.0

## ---- sites_meta (sheet B) ----
b = pd.read_excel(XL, sheet_name='B_parameters_info', header=None)
hdr = b.index[b[0] == 'Sample code'][0]
meta = b.iloc[hdr+1:].copy(); meta.columns = list(b.iloc[hdr])
meta = meta.rename(columns={'Sample code':'site_id','Latitude':'lat','Longitude':'lon',
                            'Country':'country','River basin':'river_basin'})
meta['lat'] = pd.to_numeric(meta['lat'], errors='coerce')
meta['lon'] = pd.to_numeric(meta['lon'], errors='coerce')
meta['UDF'] = pd.to_numeric(meta['Urban discharge fraction (UDF)'], errors='coerce')
meta['stream_order'] = pd.to_numeric(meta['Stream order by Strahler'], errors='coerce')
meta = meta.reset_index(drop=True)

## ---- conc_mat (sheet A) ----
ab = pd.read_excel(XL, sheet_name='A_B_results', header=None)
name_row = ab.index[ab[0] == 'Name'][0]
first_eus = next(i for i, v in enumerate(ab.iloc[name_row]) if isinstance(v, str) and v.startswith('EUS_'))
data_rows = ab.index[name_row+1:]
compounds = ab.iloc[name_row+1:, 0].astype(str).tolist()
inchikeys = ab.iloc[name_row+1:, 5].astype(str).tolist()
site_cols = list(range(first_eus, ab.shape[1]))
site_ids  = [str(v) for v in ab.iloc[name_row, first_eus:]]
conc = ab.iloc[name_row+1:, first_eus:].apply(pd.to_numeric, errors='coerce').to_numpy(dtype=float)
conc = np.nan_to_num(conc, nan=0.0)
conc_mat = pd.DataFrame(conc, index=compounds, columns=site_ids)

## ---- MDL/2 floor ----
mdl_tab = pd.read_csv(MDL)
mdl_by_ik   = dict(zip(mdl_tab['inchikey'].astype(str).str.strip().str.upper(), mdl_tab['mdl_value_ng_L']))
mdl_by_name = dict(zip(mdl_tab['chemical_name'].astype(str).str.strip().str.lower(), mdl_tab['mdl_value_ng_L']))
mdl_ngL = np.array([mdl_by_ik.get(ik.strip().upper(), np.nan) for ik in inchikeys])
unm = np.isnan(mdl_ngL)
mdl_ngL[unm] = [mdl_by_name.get(c.strip().lower(), np.nan) for c in np.array(compounds)[unm]]
unm = np.isnan(mdl_ngL)
mdl_half = mdl_ngL / 2.0
mdl_half[unm] = DEFAULT_PSEUDO_CONC
mdl_half_ngL = pd.Series(mdl_half, index=compounds)
print(f'MDL match: {(~unm).sum()}/{len(compounds)} matched; {unm.sum()} fallback +1 ng/L')

## ---- EC10 + tu_mat ----
pred = pd.read_csv(PRED); mapf = pd.read_csv(MAP)
mg = mapf.merge(pred[['SMILES', 'predictions (mg/L)']], on='SMILES', how='left')
ec10_ngL = pd.Series((mg['predictions (mg/L)'] * 1e6).to_numpy(), index=mg['cmpdname'])
ec10_ngL = ec10_ngL[np.isfinite(ec10_ngL)]
have_ec10 = [c for c in conc_mat.index if c in ec10_ngL.index]
tu = conc_mat.loc[have_ec10].to_numpy() / ec10_ngL.loc[have_ec10].to_numpy()[:, None]
tu[~np.isfinite(tu)] = 0
tu_mat = pd.DataFrame(tu, index=have_ec10, columns=conc_mat.columns)
sum_tu = tu_mat.sum(axis=0)   ## colSums(tu_mat) -- the S6/S7 stress axis
print(f'EC10: {len(ec10_ngL)} compounds; tu_mat {tu_mat.shape}')

## ---- units ----
basin_counts = meta['river_basin'].value_counts()
big_basins = basin_counts[basin_counts >= BASIN_MIN_N].index.tolist()
REF_SITES = meta.loc[meta['river_basin'].isin(big_basins), 'site_id'].tolist()
ref = [s for s in REF_SITES if s in conc_mat.columns]
sub = conc_mat[ref]
keep_glob = sub.index[(sub.to_numpy() > 0).sum(axis=1) >= max(2, int(np.ceil(PREV_THRESH * sub.shape[1])))]
KEEP_GLOBAL = list(keep_glob)
print(f'global node set: conc {len(KEEP_GLOBAL)}/{conc_mat.shape[0]} over {len(ref)} ref sites')

groups = [
    dict(label='Upper Danube',  safe='Danube_Upper',  basin='Danube', countries=['Germany','Austria','Czechia']),
    dict(label='Middle Danube', safe='Danube_Middle', basin='Danube', countries=['Slovenia','Slovakia','Hungary','Croatia']),
    dict(label='Lower Danube',  safe='Danube_Lower',  basin='Danube', countries=['Serbia','Bulgaria','Romania','Moldova','Ukraine']),
    dict(label='Upper Rhine',   safe='Rhine_Upper',   basin='Rhine',  countries=['Switzerland']),
    dict(label='Lower Rhine',   safe='Rhine_Lower',   basin='Rhine',  countries=['Germany','Netherlands'])]
for g in groups:
    g['site_ids'] = meta.loc[(meta['river_basin'] == g['basin']) & (meta['country'].isin(g['countries'])), 'site_id'].tolist()
units = [g for g in groups if len(g['site_ids']) >= BASIN_MIN_N]

basin_units = [dict(label=f'{bn} (whole basin)', safe='Basin_' + ''.join(c if c.isalnum() else '_' for c in bn),
                    basin=bn, site_ids=meta.loc[meta['river_basin'] == bn, 'site_id'].tolist()) for bn in big_basins]

SECTION_SPEC = {'Elbe': ('lat', 'low'), 'Rhine': ('lat', 'low'), 'Danube': ('lon', 'low')}
section_units = []
for bn, (ax, up) in SECTION_SPEC.items():
    if bn not in big_basins: continue
    sm = meta[(meta['river_basin'] == bn) & np.isfinite(meta[ax])]
    k = min(3, int(np.floor(len(sm) / SECTION_MIN_N)))
    if k < 2: continue
    v = sm[ax].to_numpy()
    qs = np.quantile(v, [i / k for i in range(1, k)])
    tile = pd.cut(v, bins=[-np.inf] + list(qs) + [np.inf], labels=[f't{i+1}' for i in range(k)])
    pos_names = ['Upper', 'Lower'] if k == 2 else ['Upper', 'Middle', 'Lower']
    posmap = dict(zip([f't{i+1}' for i in range(k)], pos_names if up == 'low' else pos_names[::-1]))
    for pos in pos_names:
        key = [t for t, p in posmap.items() if p == pos][0]
        section_units.append(dict(label=f'{pos} {bn} (section)',
                                  safe='Sec_' + ''.join(c if c.isalnum() else '_' for c in bn) + '_' + pos,
                                  basin=bn, site_ids=sm['site_id'][tile == key].tolist()))
section_units = [u for u in section_units if len(u['site_ids']) >= BASIN_MIN_N]
all_units = {u['safe']: u for u in units + basin_units + section_units}
print('units:', {k: len(v['site_ids']) for k, v in all_units.items()})

## ---- prep_unit (conc) ----
def prep_unit(site_ids):
    sites = [s for s in site_ids if s in conc_mat.columns]
    sub = conc_mat[sites]
    n_det = (sub.to_numpy() > 0).sum(axis=1)
    keep = sub.index.isin(KEEP_GLOBAL) & (n_det >= MIN_DETECT_UNIT)
    sub = sub.loc[keep]
    fl = mdl_half_ngL.loc[sub.index].to_numpy()
    arr = sub.to_numpy()
    arr[arr == 0] = np.repeat(fl[:, None], arr.shape[1], axis=1)[arr == 0]
    Xlog = np.log10(arr).T   ## sites x compounds
    return pd.DataFrame(Xlog, index=sites, columns=sub.index)

prepped = {safe: prep_unit(u['site_ids']) for safe, u in all_units.items()}
print('prepped:', {k: v.shape for k, v in prepped.items()})

## ---- validate against S7 site_table ----
st = pd.read_csv('/mnt/agents/output/s7_site_level/site_table.csv')
report = {}
ok = True
for safe, df in prepped.items():
    sub_st = st[st['unit'] == safe]
    row = {'n_sites_py': df.shape[0], 'n_sites_s7': len(sub_st), 'n_comp': df.shape[1]}
    if len(sub_st) != df.shape[0]: ok = False; row['MISMATCH'] = 'site count'
    stu = sum_tu.reindex(df.index).to_numpy()
    s7 = sub_st.set_index('site').reindex(df.index)
    row['sum_tu_maxdiff'] = float(np.nanmax(np.abs(stu - s7['sum_tu'].to_numpy())))
    row['udf_match'] = bool(np.allclose(meta.set_index('site_id').reindex(df.index)['UDF'].to_numpy(dtype=float),
                                        s7['UDF'].to_numpy(dtype=float), equal_nan=True))
    row['n_det_match'] = bool((((df.to_numpy() > np.log10(mdl_half_ngL.loc[df.columns].to_numpy())[None,:] + 1e-12)).sum(axis=1) == s7['n_det'].to_numpy()).all()) \
        if 'n_det' in s7 else None
    row['mean_log_maxdiff'] = float(np.max(np.abs(df.to_numpy().mean(axis=1) - s7['mean_log'].to_numpy())))
    if row['sum_tu_maxdiff'] > 1e-6 or not row['udf_match'] or row['mean_log_maxdiff'] > 1e-9: ok = False
    report[safe] = row
print(json.dumps(report, indent=1))
print('PREP VALIDATION:', 'PASS' if ok else 'FAIL')

## ---- save ----
np.savez(f'{OUT}/prep.npz',
         **{f'Xlog__{k}': v.to_numpy() for k, v in prepped.items()},
         **{f'sites__{k}': v.index.to_numpy() for k, v in prepped.items()},
         **{f'comps__{k}': v.columns.to_numpy() for k, v in prepped.items()},
         sum_tu_sites=sum_tu.index.to_numpy(), sum_tu=sum_tu.to_numpy(),
         meta_sites=meta['site_id'].to_numpy(),
         meta_basin=meta['river_basin'].astype(str).to_numpy(),
         meta_udf=meta['UDF'].to_numpy(dtype=float),
         meta_so=meta['stream_order'].to_numpy(dtype=float),
         units_json=np.array([json.dumps({k: v['site_ids'] for k, v in all_units.items()})]))
print('saved', f'{OUT}/prep.npz')
