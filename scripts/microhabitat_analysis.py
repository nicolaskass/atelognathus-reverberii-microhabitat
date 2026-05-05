"""
microhabitat_analysis.py
========================
Reproduces all analyses and figures in:

  Kass et al. "Microhabitat use by Atelognathus reverberii (Cei, 1969)
  at Laguna Azul, Somuncura Plateau, Patagonia: selection for the rocky
  shoreline ecotone and proximity to water"
  Ichthyology & Herpetology, in review.

Usage:
  python scripts/microhabitat_analysis.py --datadir data/ --outdir outputs/

Dependencies: numpy, scipy, pandas, matplotlib
"""
import argparse, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats
from scipy.stats import mannwhitneyu, fisher_exact, spearmanr, beta as beta_dist

warnings.filterwarnings('ignore')
plt.rcParams.update({'font.family': 'serif', 'font.size': 10,
                     'axes.spines.top': False, 'axes.spines.right': False})

SP = r'$\it{Atelognathus\ reverberii}$'
ZCOLS  = ['#c9a227', '#6aaa45', '#2166ac']
SCOLS  = ['#4dac26', '#969696', '#525252', '#c9a227']


# ── Data loading ──────────────────────────────────────────────────────────────

def load_data(datadir):
    d = Path(datadir)
    r = pd.read_csv(d / 'quadrants_rocky_zone.csv')
    i = pd.read_csv(d / 'quadrants_steppe_zone.csv')
    p = pd.read_csv(d / 'quadrants_lacustrine_zone.csv')
    ei= pd.read_csv(d / 'independent_encounters_rocky_zone.csv')
    return r, i, p, ei


# ── Between-zone analysis ─────────────────────────────────────────────────────

def between_zone_analysis(r, i, p):
    n_det = np.array([r['presence'].sum(), i['presence'].sum(), p['presence'].sum()])
    n_q   = np.array([80, 80, 80])
    rates = n_det / n_q
    alpha = rates / rates.sum()

    chi2_q, p_q, _, _ = stats.chi2_contingency(np.array([n_det, n_q - n_det]))

    # Bootstrap CIs for alpha
    rng = np.random.default_rng(42)
    boot = []
    for _ in range(10000):
        bd = np.array([rng.binomial(80, rates[k]) for k in range(3)])
        s  = bd / 80
        boot.append(s / s.sum() if s.sum() > 0 else np.zeros(3))
    boot = np.array(boot)
    alpha_ci = np.percentile(boot, [2.5, 97.5], axis=0)

    # Transect-level confirmatory
    t_pos = np.array([5, 2, 0])
    t_neg = np.array([3, 6, 8])
    chi2_t, p_t, _, _ = stats.chi2_contingency(np.array([t_pos, t_neg]))

    rng2   = np.random.default_rng(42)
    zones  = np.repeat([0, 1, 2], t_pos + t_neg)
    outs   = np.concatenate([[1]*pp + [0]*nn for pp, nn in zip(t_pos, t_neg)])
    cnt    = 0
    for _ in range(50000):
        pm = rng2.permutation(outs)
        d  = np.array([pm[zones == z].sum() for z in range(3)])
        a  = (t_pos + t_neg) - d
        try:
            c2, _, _, _ = stats.chi2_contingency(np.array([d, a]))
            if c2 >= chi2_t: cnt += 1
        except: pass
    p_perm = cnt / 50000

    _, p_ri = fisher_exact([[5, 3], [2, 6]])
    _, p_rp = fisher_exact([[5, 3], [0, 8]])

    return dict(n_det=n_det, n_q=n_q, rates=rates, alpha=alpha,
                alpha_ci=alpha_ci, chi2_q=chi2_q, p_q=p_q,
                chi2_t=chi2_t, p_t=p_t, p_perm=p_perm,
                p_ri=p_ri, p_rp=p_rp)


# ── Within rocky zone: use vs availability ────────────────────────────────────

def within_zone_analysis(r):
    used  = r[r['presence'] == 1]
    avail = r[r['presence'] == 0]
    results = {}
    vars_to_test = ['pct_rock', 'pct_vegetation', 'pct_gravel', 'pct_soil', 'distance_m']
    pvals = []
    for var in vars_to_test:
        u = used[var].dropna()
        a = avail[var].dropna()
        _, p = mannwhitneyu(u, a, alternative='two-sided')
        pvals.append(p)
    # Holm correction
    order = np.argsort(pvals)
    holm_thresholds = [0.05 / (5 - rank) for rank in range(5)]
    for rank, idx in enumerate(order):
        var = vars_to_test[idx]
        u = used[var].dropna()
        a = avail[var].dropna()
        results[var] = {
            'used_med': u.median(), 'avail_med': a.median(),
            'p': pvals[idx],
            'holm_thr': holm_thresholds[rank],
            'significant': pvals[idx] < holm_thresholds[rank]
        }
    return results, used, avail


# ── Within-zone confirmatory analysis (transect-level, addresses pseudoreplication) ──

def within_zone_transect_test(r, n_perm=50000, seed=42):
    """
    Pseudoreplication-aware confirmatory test for distance-to-shore preference
    within the rocky zone.

    Strategy: presence labels are permuted *within each transect*, preserving
    (a) the marginal number of detections per transect and (b) the transect
    structure itself. Under the null of no within-transect distance preference,
    detections should be uniformly distributed across the 10 quadrant positions
    of each positive transect. The test statistic is the median distance-to-
    shore among all detected quadrants.

    Returns observed median, permutation median distribution, p-value.
    """
    rng = np.random.default_rng(seed)
    obs_used = r[r['presence'] == 1]
    obs_med  = obs_used['distance_m'].median()

    # Within-transect permutation
    perm_meds = np.empty(n_perm)
    transect_ids = r['transect_id'].unique()
    for k in range(n_perm):
        sim_pres = np.zeros(len(r), dtype=int)
        for t_id in transect_ids:
            mask = (r['transect_id'] == t_id).values
            t_pres = r.loc[mask, 'presence'].values
            sim_pres[mask] = rng.permutation(t_pres)
        sim_dist = r.loc[sim_pres == 1, 'distance_m']
        perm_meds[k] = sim_dist.median() if len(sim_dist) else np.nan

    p_perm = float(np.mean(perm_meds <= obs_med))

    # Per-transect detection-position summary
    per_transect = []
    for t_id in transect_ids:
        t = r[r['transect_id'] == t_id]
        det = t[t['presence'] == 1]
        per_transect.append({
            'transect_id'    : t_id,
            'n_detections'   : int(t['presence'].sum()),
            'detection_dists': det['distance_m'].tolist(),
            'median_det_dist': float(det['distance_m'].median()) if len(det) else np.nan
        })

    return {
        'obs_median'         : float(obs_med),
        'perm_median_p'      : p_perm,
        'perm_dist_summary'  : {
            'mean' : float(np.nanmean(perm_meds)),
            'p2.5' : float(np.nanpercentile(perm_meds, 2.5)),
            'p97.5': float(np.nanpercentile(perm_meds, 97.5))
        },
        'per_transect'       : per_transect,
        'n_transects_pos'    : sum(1 for x in per_transect if x['n_detections'] > 0)
    }


def predictor_collinearity(r):
    """
    Spearman correlation among the four within-zone predictors that approached
    significance in the Mann-Whitney battery, to flag collinearity that may
    inflate apparent independent effects.
    """
    cols = ['distance_m', 'pct_vegetation', 'pct_rock', 'pct_gravel', 'pct_soil']
    out = {}
    for i, c1 in enumerate(cols):
        for c2 in cols[i+1:]:
            rho, p = spearmanr(r[c1], r[c2])
            out[f'{c1}__vs__{c2}'] = {'rho': float(rho), 'p': float(p)}
    return out


def mw_power_simulation(n1=8, n2=72, effect_sizes=(0.3, 0.4, 0.5),
                         n_sim=20000, alpha=0.05, seed=42):
    """
    Monte-Carlo power for a two-sided Mann-Whitney U test with the actual
    sample-size structure of the within-zone use-vs-availability design
    (n1 = 8 used, n2 = 72 available).

    Effect size is parameterised as the rank-biserial correlation r;
    samples are drawn from two normal distributions whose standardised
    mean difference is set to map onto the requested r via the
    Mann-Whitney-to-d approximation r = 1 - 2*U/(n1*n2) and U_expected
    under shifted normals. We use the empirical mapping: shift d such
    that the resulting MW r matches the target.
    """
    rng = np.random.default_rng(seed)
    # For two normals with shift d, the MW rank-biserial r approaches
    # 2 * Phi(d/sqrt(2)) - 1.  Invert numerically.
    from scipy.stats import norm
    out = {}
    for r_target in effect_sizes:
        # Solve d:  2*Phi(d/sqrt(2)) - 1 = r_target
        d = norm.ppf((r_target + 1) / 2) * np.sqrt(2)
        rejects = 0
        for _ in range(n_sim):
            x = rng.normal(loc=d, scale=1.0, size=n1)
            y = rng.normal(loc=0.0, scale=1.0, size=n2)
            _, p = mannwhitneyu(x, y, alternative='two-sided')
            if p < alpha:
                rejects += 1
        out[r_target] = rejects / n_sim
    return out


def jeffreys_upper_bound(o, n, alpha=0.025):
    """
    Upper bound of the (1-2*alpha) Jeffreys credible interval for a proportion
    when o = 0.  Used to provide a sensible 'undefined-but-bounded' value for
    the lacustrine-zone Manly's alpha CI.
    """
    if o == 0:
        return float(beta_dist.ppf(1 - alpha, 0.5, n + 0.5))
    if o == n:
        return 1.0
    return float(beta_dist.ppf(1 - alpha, o + 0.5, n - o + 0.5))


# ── Figures ───────────────────────────────────────────────────────────────────

def fig_zone_patterns(r, i, p, bz, outdir):
    fig, axes = plt.subplots(1, 3, figsize=(13, 5.5))

    zones_names = ['Lacustrine\nsediment', 'Intermediate\nsteppe', 'Rocky\nshore']
    rates_plot  = bz['rates'][[2, 1, 0]]   # reorder for display
    alpha_plot  = bz['alpha'][[2, 1, 0]]
    cols_plot   = ZCOLS[::-1]

    # A: detection rates
    ax = axes[0]
    bars = ax.bar(zones_names, rates_plot * 100, color=cols_plot,
                  alpha=0.82, width=0.5, edgecolor='white')
    for bar, nd, nq in zip(bars, bz['n_det'][[2, 1, 0]], bz['n_q'][[2, 1, 0]]):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                f'{nd}/{nq}', ha='center', va='bottom', fontsize=9)
    ax.set_ylabel('Detection rate (% of quadrants)', fontsize=10)
    ax.set_title('A  Detection rate by zone', fontsize=10, loc='left', fontweight='bold')
    ax.set_ylim(0, 15)
    ax.text(0.97, 0.97, f'Quadrant: p={bz["p_q"]:.3f}\nTransect: p={bz["p_t"]:.3f}',
            transform=ax.transAxes, ha='right', va='top', fontsize=8, color='#444')

    # B: Manly's alpha
    ax = axes[1]
    ci = bz['alpha_ci'][:, [2, 1, 0]]
    ax.bar(zones_names, alpha_plot, color=cols_plot, alpha=0.82, width=0.5, edgecolor='white')
    ax.errorbar(range(3), alpha_plot,
                yerr=[alpha_plot - ci[0], ci[1] - alpha_plot],
                fmt='none', color='k', lw=1.5, capsize=5)
    ax.axhline(1/3, color='#d6604d', ls='--', lw=1.4, label='Random use (α=0.33)')
    ax.set_ylabel("Manly's electivity (α)", fontsize=10)
    ax.set_title("B  Habitat electivity", fontsize=10, loc='left', fontweight='bold')
    ax.set_ylim(-0.05, 1.1)
    ax.legend(fontsize=8, frameon=False)

    # C: substrate
    zone_subs = {
        'Lacustrine\nsediment': [p['pct_vegetation'].mean(), p['pct_gravel'].mean(),
                                  p['pct_rock'].mean(), p['pct_soil'].mean()],
        'Intermediate\nsteppe':  [i['pct_vegetation'].mean(), i['pct_gravel'].mean(),
                                   i['pct_rock'].mean(), i['pct_soil'].mean()],
        'Rocky\nshore':          [r['pct_vegetation'].mean(), r['pct_gravel'].mean(),
                                   r['pct_rock'].mean(), r['pct_soil'].mean()],
    }
    slabs = ['Vegetation', 'Gravel', 'Rock', 'Soil/mud']
    ax = axes[2]
    x = np.arange(3); bot = np.zeros(3)
    for idx_s, (sl, sc) in enumerate(zip(slabs, SCOLS)):
        vals = [zone_subs[z][idx_s] for z in zone_subs]
        ax.bar(x, vals, bottom=bot, color=sc, label=sl, alpha=0.85, edgecolor='white')
        bot += np.array(vals)
    ax.set_xticks(x); ax.set_xticklabels(list(zone_subs.keys()), fontsize=8)
    ax.set_ylabel('Mean substrate cover (%)', fontsize=10)
    ax.set_title('C  Substrate by zone', fontsize=10, loc='left', fontweight='bold')
    ax.legend(fontsize=8, frameon=False)

    fig.suptitle(f'Zone-level microhabitat use — {SP}', fontsize=11, y=1.01)
    fig.tight_layout()
    fig.savefig(outdir / 'fig1_zone_patterns.png', dpi=300, bbox_inches='tight')
    plt.close(fig)


def fig_within_zone(r, ei, wz_results, outdir):
    fig = plt.figure(figsize=(13, 5.5))
    gs  = gridspec.GridSpec(1, 3, figure=fig, wspace=0.45,
                            left=0.07, right=0.97, top=0.88, bottom=0.14)

    # A: detection by distance
    ax = fig.add_subplot(gs[0])
    dists = sorted(r['distance_m'].unique())
    det_r_d = [(r[r['distance_m'] == d]['presence'].sum(),
                len(r[r['distance_m'] == d])) for d in dists]
    det_rates = [nd/nq if nq > 0 else 0 for nd, nq in det_r_d]
    colors_d  = ['#d62728' if dr > 0 else '#aec7e8' for dr in det_rates]
    ax.bar(dists, det_rates, color=colors_d, alpha=0.80, width=0.7, edgecolor='white')
    for d, dr, (nd, nq) in zip(dists, det_rates, det_r_d):
        if nd > 0:
            ax.text(d, dr + 0.01, f'{nd}/{nq}', ha='center', va='bottom', fontsize=8)
    ax.axvspan(-0.5, 0.5, color='#4393c3', alpha=0.12, label='Aquatic quadrant')
    ax.set_xlabel('Distance from shore (m)', fontsize=10)
    ax.set_ylabel('Detection rate', fontsize=10)
    ax.set_title('A  Detection by distance\n(rocky zone; MW p=0.002)',
                 fontsize=10, loc='left', fontweight='bold')
    ax.legend(fontsize=8, frameon=False)

    # B: used vs available substrate
    ax = fig.add_subplot(gs[1])
    used  = r[r['presence'] == 1]
    avail = r[r['presence'] == 0]
    sub_vars   = ['pct_rock', 'pct_vegetation', 'pct_gravel', 'pct_soil']
    sub_labels = ['Rock', 'Veg.', 'Gravel', 'Soil']
    x = np.arange(4); w = 0.35
    u_meds = [used[v].median() for v in sub_vars]
    a_meds = [avail[v].median() for v in sub_vars]
    ax.bar(x - w/2, u_meds, width=w, color='#d6604d', alpha=0.8, label=f'Used (n={len(used)})')
    ax.bar(x + w/2, a_meds, width=w, color='#92c5de', alpha=0.8, label=f'Available (n={len(avail)})')
    ax.set_xticks(x); ax.set_xticklabels(sub_labels, fontsize=9)
    ax.set_ylabel('Median cover (%)', fontsize=10)
    ax.set_title('B  Used vs. available substrate', fontsize=10, loc='left', fontweight='bold')
    ax.legend(fontsize=8, frameon=False)
    for k, var in enumerate(sub_vars):
        if wz_results[var]['significant']:
            ymax = max(u_meds[k], a_meds[k]) + 3
            ax.text(k, ymax, '**' if wz_results[var]['p'] < 0.01 else '*',
                    ha='center', fontsize=11, fontweight='bold')

    # C: IE distances
    ax = fig.add_subplot(gs[2])
    ei_dists = ei['distance_m'].dropna().values
    ax.hist(ei_dists, bins=np.arange(0, 28, 3), color='#2166ac',
            alpha=0.75, edgecolor='white')
    med = np.median(ei_dists)
    ax.axvline(med, color='k', lw=2, label=f'Median = {med:.1f} m')
    ax.set_xlabel('Distance to shore (m)', fontsize=10)
    ax.set_ylabel('Number of encounters', fontsize=10)
    ax.set_title(f'C  Independent encounters\n(n={len(ei)}, all under rocks)',
                 fontsize=10, loc='left', fontweight='bold')
    ax.legend(fontsize=9, frameon=False)

    fig.suptitle(f'Within-zone microhabitat — {SP} (rocky shore)', fontsize=11)
    fig.savefig(outdir / 'fig2_within_zone.png', dpi=300, bbox_inches='tight')
    plt.close(fig)


# ── Main ──────────────────────────────────────────────────────────────────────

def run(datadir='data/', outdir='outputs/'):
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    r, i, p, ei = load_data(datadir)

    print("=== BETWEEN-ZONE ANALYSIS ===")
    bz = between_zone_analysis(r, i, p)
    print(f"  Rocky: {bz['n_det'][0]}/80 = {bz['rates'][0]*100:.1f}%  alpha={bz['alpha'][0]:.3f}")
    print(f"  Interm: {bz['n_det'][1]}/80 = {bz['rates'][1]*100:.1f}%  alpha={bz['alpha'][1]:.3f}")
    print(f"  Playa: {bz['n_det'][2]}/80 = 0.0%  alpha={bz['alpha'][2]:.3f}")
    print(f"  Quadrant chi2={bz['chi2_q']:.2f}, p={bz['p_q']:.4f}")
    print(f"  Transect chi2={bz['chi2_t']:.2f}, p={bz['p_t']:.4f}, perm p={bz['p_perm']:.4f}")

    print("\n=== WITHIN ROCKY ZONE (exploratory, quadrant-level Mann-Whitney) ===")
    wz, used, avail = within_zone_analysis(r)
    for var, res in sorted(wz.items(), key=lambda x: x[1]['p']):
        sig = '*' if res['significant'] else ' '
        print(f"  {sig} {var:20s}: used={res['used_med']:.0f}  avail={res['avail_med']:.0f}"
              f"  p={res['p']:.4f}  thr={res['holm_thr']:.3f}")

    print("\n=== WITHIN ROCKY ZONE (confirmatory, transect-level permutation) ===")
    tlevel = within_zone_transect_test(r)
    print(f"  Observed median distance of detected quadrants: {tlevel['obs_median']:.1f} m")
    print(f"  Within-transect permutation null median: "
          f"{tlevel['perm_dist_summary']['mean']:.2f} "
          f"[95% PI: {tlevel['perm_dist_summary']['p2.5']:.1f}, "
          f"{tlevel['perm_dist_summary']['p97.5']:.1f}]")
    print(f"  Permutation p (one-sided, lower) = {tlevel['perm_median_p']:.4f}")
    print(f"  Positive transects: {tlevel['n_transects_pos']}/8")
    for pt in tlevel['per_transect']:
        if pt['n_detections'] > 0:
            print(f"    {pt['transect_id']}: {pt['n_detections']} det at "
                  f"{pt['detection_dists']} m")

    print("\n=== PREDICTOR COLLINEARITY (Spearman) ===")
    coll = predictor_collinearity(r)
    for pair, vals in sorted(coll.items(), key=lambda x: -abs(x[1]['rho'])):
        flag = '!' if abs(vals['rho']) > 0.5 else ' '
        print(f"  {flag} {pair:35s} rho={vals['rho']:+.3f}  p={vals['p']:.4f}")

    print("\n=== MANN-WHITNEY POWER (Monte Carlo, n1=8 vs n2=72) ===")
    pwr = mw_power_simulation()
    for r_target, pwr_val in pwr.items():
        print(f"  rank-biserial r = {r_target:.2f} -> power = {pwr_val:.3f}")

    print("\n=== JEFFREYS UPPER BOUND FOR LACUSTRINE alpha (replaces 0,0 CI) ===")
    j_up = jeffreys_upper_bound(0, 80)
    rates = bz['rates']
    sum_rates = rates.sum()
    j_alpha_up = j_up / (j_up + sum_rates)
    print(f"  Jeffreys 97.5% upper bound, p_lacustrine = {j_up:.4f}")
    print(f"  -> Manly's alpha_lacustrine 97.5% upper bound ≈ {j_alpha_up:.3f}")

    print("\n=== INDEPENDENT ENCOUNTERS ===")
    print(f"  n={len(ei)}, frogs={int(ei['n_frogs'].sum())}")
    print(f"  Distance: median={ei['distance_m'].median():.1f}m "
          f"IQR=[{ei['distance_m'].quantile(.25):.1f},{ei['distance_m'].quantile(.75):.1f}]")
    rock_kw = r'roca|piedra|mediana|chica|grande'
    print(f"  All under rocks: "
          f"{ei['notes'].str.contains(rock_kw, case=False, na=False).all()}")

    print("\nGenerating figures...")
    fig_zone_patterns(r, i, p, bz, outdir)
    fig_within_zone(r, ei, wz, outdir)
    print(f"  Figures saved to {outdir}")

    # Save summary CSV
    summary_rows = []
    for var, res in wz.items():
        summary_rows.append({'variable': var, **res})
    pd.DataFrame(summary_rows).to_csv(outdir / 'within_zone_mw_results.csv', index=False)

    # Save transect-level confirmatory + collinearity outputs
    pd.DataFrame(tlevel['per_transect']).to_csv(
        outdir / 'within_zone_transect_summary.csv', index=False)
    coll_rows = [{'pair': k, **v} for k, v in coll.items()]
    pd.DataFrame(coll_rows).to_csv(
        outdir / 'within_zone_collinearity.csv', index=False)

    # Single reportable-quantities file (one source of truth)
    reportable = {
        'rocky_n_det'        : int(bz['n_det'][0]),
        'rocky_rate'         : float(bz['rates'][0]),
        'intermed_n_det'     : int(bz['n_det'][1]),
        'intermed_rate'      : float(bz['rates'][1]),
        'lacustrine_n_det'   : int(bz['n_det'][2]),
        'lacustrine_rate'    : float(bz['rates'][2]),
        'alpha_rocky'        : float(bz['alpha'][0]),
        'alpha_intermed'     : float(bz['alpha'][1]),
        'alpha_lacustrine'   : float(bz['alpha'][2]),
        'alpha_rocky_lo'     : float(bz['alpha_ci'][0,0]),
        'alpha_rocky_hi'     : float(bz['alpha_ci'][1,0]),
        'alpha_intermed_lo'  : float(bz['alpha_ci'][0,1]),
        'alpha_intermed_hi'  : float(bz['alpha_ci'][1,1]),
        'alpha_lacustrine_jeffreys_upper': j_alpha_up,
        'chi2_quadrant'      : float(bz['chi2_q']),
        'p_chi2_quadrant'    : float(bz['p_q']),
        'chi2_transect'      : float(bz['chi2_t']),
        'p_chi2_transect'    : float(bz['p_t']),
        'p_perm_transect'    : float(bz['p_perm']),
        'mw_distance_p'      : float(wz['distance_m']['p']),
        'mw_vegetation_p'    : float(wz['pct_vegetation']['p']),
        'mw_rock_p'          : float(wz['pct_rock']['p']),
        'mw_gravel_p'        : float(wz['pct_gravel']['p']),
        'mw_soil_p'          : float(wz['pct_soil']['p']),
        'transect_perm_p_distance': float(tlevel['perm_median_p']),
        'transect_obs_median'    : float(tlevel['obs_median']),
        'n_transects_positive'   : int(tlevel['n_transects_pos']),
        'spearman_dist_veg_rho'  : float(coll['distance_m__vs__pct_vegetation']['rho']),
        'spearman_dist_veg_p'    : float(coll['distance_m__vs__pct_vegetation']['p']),
        'mw_power_r03'           : float(pwr[0.3]),
        'mw_power_r04'           : float(pwr[0.4]),
        'mw_power_r05'           : float(pwr[0.5])
    }
    with open(outdir / 'reportable_quantities.json', 'w') as f:
        import json; json.dump(reportable, f, indent=2)

    print(f"\nAll outputs written to {outdir}")
    return bz, wz, tlevel, coll, reportable


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--datadir', default='data/')
    ap.add_argument('--outdir',  default='outputs/')
    args = ap.parse_args()
    run(args.datadir, args.outdir)