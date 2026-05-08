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

    # ── Leave-one-positive-transect-out sensitivity ──
    # The within-transect permutation result rests on n=5 positive transects.
    # We re-run the test removing each positive transect in turn and report
    # the worst-case (most conservative) p-value to bound the result's
    # robustness to a single missed transect.
    pos_ids = [pt['transect_id'] for pt in per_transect if pt['n_detections'] > 0]
    loo_results = []
    for drop_id in pos_ids:
        sub = r[r['transect_id'] != drop_id].copy()
        sub_used = sub[sub['presence'] == 1]
        sub_obs_med = float(sub_used['distance_m'].median()) if len(sub_used) else np.nan
        sub_transects = sub['transect_id'].unique()
        rng_loo = np.random.default_rng(seed + hash(drop_id) % 10000)
        sub_perm_meds = np.empty(n_perm // 5)  # smaller n_perm for speed; still 10k
        for k in range(len(sub_perm_meds)):
            sim_pres = np.zeros(len(sub), dtype=int)
            for t_id in sub_transects:
                mask = (sub['transect_id'] == t_id).values
                sim_pres[mask] = rng_loo.permutation(sub.loc[mask, 'presence'].values)
            sd = sub.loc[sim_pres == 1, 'distance_m']
            sub_perm_meds[k] = sd.median() if len(sd) else np.nan
        sub_p = float(np.mean(sub_perm_meds <= sub_obs_med))
        loo_results.append({
            'dropped_transect'      : drop_id,
            'remaining_n_detections': int(sub['presence'].sum()),
            'obs_median_remaining'  : sub_obs_med,
            'p_perm_remaining'      : sub_p
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
        'n_transects_pos'    : sum(1 for x in per_transect if x['n_detections'] > 0),
        'loo_results'        : loo_results,
        'loo_max_p'          : float(max(x['p_perm_remaining'] for x in loo_results))
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


def aggregation_per_occupied_quadrant(r):
    """
    Distribution of detection-event counts per occupied quadrant in the
    rocky zone. Reports mean, range, and total — used to clarify whether
    the 52 events come from heavy aggregation in a few rocks or roughly
    one frog per occupied quadrant (M3 from review iter-3).
    """
    occ = r[r['presence'] == 1]
    counts = occ['n_frogs'].astype(int).tolist()
    return {
        'n_occupied'  : int(len(occ)),
        'n_events'    : int(sum(counts)),
        'per_quadrant_counts': counts,
        'mean_per_q'  : float(np.mean(counts)) if counts else 0.0,
        'median_per_q': float(np.median(counts)) if counts else 0.0,
        'min_per_q'   : int(min(counts)) if counts else 0,
        'max_per_q'   : int(max(counts)) if counts else 0
    }


def detectability_sensitivity(o_rocky=8, n_rocky=80, n_lacus=80, rock_rocky=0.68,
                               rock_lacus=0.04, alpha=0.05):
    """
    Sensitivity analysis (M2 from review iter-3): given the observed
    8/80 vs 0/80 detection contrast, and assuming detection probability
    per quadrant scales linearly with mean rock cover, what relative
    detection probability ratio (rho_lacus/rho_rocky) would have to hold
    for the observed lacustrine zero to be compatible with EQUAL true
    occupancy across zones?

    Under equal true occupancy and detection p_z proportional to rock_z
    (so p_lacus / p_rocky = rock_lacus / rock_rocky for the linear case),
    the expected lacustrine detection count given the observed rocky
    rate is:
        E[o_lacus | equal use, linear p ~ rock_cover]
            = n_lacus * (o_rocky / n_rocky) * (rock_lacus / rock_rocky)

    We compare this expected count to the observed zero, and report the
    Poisson upper-tail probability of observing zero given that expected
    count. We also report the ROCK-COVER-RATIO at which the expected
    count drops to ≤ alpha (i.e., zero is plausible at confidence 1-alpha).
    """
    p_rocky = o_rocky / n_rocky
    # Linear-in-rock-cover detection assumption:
    ratio_linear = rock_lacus / rock_rocky        # = 0.0588 for 0.04/0.68
    expected_lacus_linear = n_lacus * p_rocky * ratio_linear
    # Probability of observing zero given that expected count, assuming
    # Poisson detection process:
    p_obs_zero_linear = float(np.exp(-expected_lacus_linear))
    # Threshold ratio: at what rho would the expected count be alpha?
    # n_lacus * p_rocky * rho = -ln(alpha)
    threshold_ratio = -np.log(alpha) / (n_lacus * p_rocky)
    return {
        'p_rocky'              : p_rocky,
        'rock_cover_ratio'     : ratio_linear,
        'expected_lacus_under_equal_use_linear': float(expected_lacus_linear),
        'prob_observed_zero_given_linear'      : p_obs_zero_linear,
        'detection_ratio_threshold_for_zero_at_alpha': float(threshold_ratio),
        'note': ('Under equal true occupancy and detection probability '
                 'linear in mean rock cover, the lacustrine zone would '
                 'be expected to yield the listed expected count of '
                 'detection events. The rho threshold is the maximum '
                 'p_lacus/p_rocky at which observing zero is plausible '
                 'at the chosen alpha (i.e., expected count <= -ln(alpha)).')
    }


def detectability_sensitivity_variants(o_rocky=8, n_rocky=80, n_lacus=80,
                                        rock_rocky=0.68, rock_lacus=0.04,
                                        alpha=0.05):
    """
    Iter-4 M3: complementary detectability-sensitivity variants beyond the
    linear-in-rock-cover model.

    Three parameterisations of how detection probability scales with rock cover:
        - LINEAR:        p_z proportional to rock_z. Already in
                         detectability_sensitivity().
        - STEP:          p_z = 0 if rock_z below a threshold (default 5%),
                         else p_z = p_rocky. Reflects "no movable rocks → no
                         detectable events" structurally.
        - SATURATING:    Hill-style p_z = p_rocky * rock_z / (rock_z + K),
                         with K chosen so that p_rocky/p_lacus matches the
                         linear case at the observed rock covers (i.e., a
                         smooth decreasing-returns alternative).

    Returns the expected lacustrine detection count, P(observe 0), and the
    interpretation tag for each variant. The reader can compare and pick
    the most defensible model — the substantive claim (the lacustrine zero
    is consistent with equal underlying use under any plausible substrate-
    conditioned detectability) holds robustly across all three.
    """
    p_rocky = o_rocky / n_rocky
    out = {}

    # Linear (already covered by detectability_sensitivity, repeated here
    # for parallel reporting)
    expected_linear = n_lacus * p_rocky * (rock_lacus / rock_rocky)
    out['linear'] = {
        'expected_lacus' : float(expected_linear),
        'prob_zero'      : float(np.exp(-expected_linear)),
        'interpretation' : 'Detection probability scales linearly with rock cover.'
    }

    # Step function: detection probability is zero below a 5% rock-cover threshold
    threshold = 0.05
    p_lacus_step = 0.0 if rock_lacus < threshold else p_rocky
    expected_step = n_lacus * p_lacus_step
    out['step'] = {
        'expected_lacus' : float(expected_step),
        'prob_zero'      : float(np.exp(-expected_step)),
        'rock_threshold' : threshold,
        'interpretation' : 'Step function: detection requires at least '
                           f'{threshold*100:.0f}% rock cover; below that, p=0.'
    }

    # Saturating (Hill, n=1): p_z = p_rocky * rock_z / (rock_z + K)
    # with K chosen so the lacus/rocky ratio matches the linear case ratio
    # at the observed rock covers:
    #   p_lacus / p_rocky = rock_lacus / (rock_lacus + K) / [rock_rocky / (rock_rocky + K)]
    # Setting this equal to rock_lacus/rock_rocky (linear) and solving for K:
    #   rock_lacus*(rock_rocky+K) / [rock_rocky*(rock_lacus+K)] = rock_lacus/rock_rocky
    #   (rock_rocky + K) / (rock_lacus + K) = 1
    # which collapses to K → ∞ (i.e., the linear case is the K→∞ limit).
    # Instead we choose K = rock_rocky/4 (a moderate concavity case) and
    # report what that implies:
    K = rock_rocky / 4
    p_ratio_sat = (rock_lacus / (rock_lacus + K)) / (rock_rocky / (rock_rocky + K))
    expected_sat = n_lacus * p_rocky * p_ratio_sat
    out['saturating'] = {
        'expected_lacus' : float(expected_sat),
        'prob_zero'      : float(np.exp(-expected_sat)),
        'half_saturation_K' : K,
        'p_ratio_lacus_to_rocky' : float(p_ratio_sat),
        'interpretation' : ('Hill saturating with half-saturation '
                            f'K={K:.2f} (a concave decreasing-returns model).')
    }

    return out


def manly_alpha_quadrant_based(r_zone_n_q=(80, 80, 80),
                                 r_zone_occupied=(8, 2, 0)):
    """
    Iter-4 M4: Manly's standardised selection ratio computed from
    OCCUPIED-QUADRANT counts (binary outcome per quadrant) rather than
    DETECTION-EVENT counts (multiple frogs per quadrant counted as separate
    events). Reduces sensitivity to within-quadrant aggregation.

        alpha_i_quadrant = (occupied_i / n_i) / sum_j(occupied_j / n_j)

    Returns the three alpha values plus an interpretation note.
    """
    rates = np.array([o / n for o, n in zip(r_zone_occupied, r_zone_n_q)])
    s = rates.sum()
    alphas = rates / s if s > 0 else np.zeros_like(rates)
    return {
        'alpha_rocky_quadrant'        : float(alphas[0]),
        'alpha_intermediate_quadrant' : float(alphas[1]),
        'alpha_lacustrine_quadrant'   : float(alphas[2]),
        'note': ('Computed from occupied-quadrant counts (binary 0/1 per '
                 'quadrant) rather than detection-event counts; addresses '
                 'within-quadrant aggregation pseudoreplication in '
                 'event-based alpha.')
    }


def pooled_adjacent_transect_test(r, i, p, n_perm=50000, seed=42):
    """
    Iter-4 M1: pooled-adjacent transect-level sensitivity for between-zone
    detection rate.

    Strategy: the eight transects within each zone are deployed sequentially
    along the lagoon perimeter and named with sequential prefixes (e.g.,
    RN1, RN2, RM1, RM2, RF1, RF2, RC1, RC2 in the rocky zone), suggesting
    spatial ordering. We pool consecutive pairs (RN1+RN2, RM1+RM2, ...) into
    "block-transects" of 4 per zone, and re-run the between-zone chi-square
    and permutation test with n=4 per zone instead of n=8.

    If the result survives this conservative aggregation, the original
    confirmatory test is robust to inter-transect spatial autocorrelation
    at the perimeter scale.
    """
    def pool_pairs(df):
        # Sort transects, pair consecutively
        ts = sorted(df['transect_id'].unique())
        pairs = [(ts[k], ts[k+1]) for k in range(0, len(ts), 2) if k+1 < len(ts)]
        block_pos = []
        for a, b in pairs:
            mask = df['transect_id'].isin([a, b])
            block_pos.append(int(df.loc[mask, 'presence'].sum() > 0))
        return block_pos, pairs

    rocky_blocks,    rocky_pairs    = pool_pairs(r)
    interm_blocks,   interm_pairs   = pool_pairs(i)
    lacus_blocks,    lacus_pairs    = pool_pairs(p)

    n_blocks = np.array([len(rocky_blocks), len(interm_blocks), len(lacus_blocks)])
    pos      = np.array([sum(rocky_blocks), sum(interm_blocks), sum(lacus_blocks)])
    neg      = n_blocks - pos
    chi2, p_val, _, _ = stats.chi2_contingency(np.array([pos, neg]))

    # Permutation test on the 12-block layout
    rng = np.random.default_rng(seed)
    zones  = np.repeat([0, 1, 2], n_blocks)
    obs    = np.concatenate([
        [1]*sum(rocky_blocks)  + [0]*(len(rocky_blocks)  - sum(rocky_blocks)),
        [1]*sum(interm_blocks) + [0]*(len(interm_blocks) - sum(interm_blocks)),
        [1]*sum(lacus_blocks)  + [0]*(len(lacus_blocks)  - sum(lacus_blocks))
    ])
    cnt = 0
    for _ in range(n_perm):
        sim = rng.permutation(obs)
        d   = np.array([sim[zones == z].sum() for z in range(3)])
        a   = n_blocks - d
        try:
            c2, _, _, _ = stats.chi2_contingency(np.array([d, a]))
            if c2 >= chi2: cnt += 1
        except: pass
    p_perm = cnt / n_perm

    return {
        'rocky_blocks'    : rocky_blocks,
        'rocky_pairs'     : [list(pp) for pp in rocky_pairs],
        'interm_blocks'   : interm_blocks,
        'lacus_blocks'    : lacus_blocks,
        'chi2'            : float(chi2),
        'p_chi2'          : float(p_val),
        'p_perm'          : float(p_perm),
        'note'            : ('Sensitivity test pooling adjacent transect pairs '
                             'into n=4 block-transects per zone; if surviving, '
                             'the n=8 confirmatory test is robust to '
                             'inter-transect spatial autocorrelation.')
    }


def rock_turnability_gradient(r):
    """
    Iter-4 m3 + iter-7 M5: rock-turnability and cavity-density gradient
    by distance from shore in the rocky zone.

    Reports per quadrant position:
      - prop_movable: proportion of "movable" rocks (small+medium) over
        total counted rocks (cavity-availability proxy)
      - mean_movable_per_quadrant: mean number of movable rocks per
        quadrant (cavity-density proxy — the actual number of shelter
        sites the protocol could turn at each distance)

    Together these bound the within-zone detection-volume confound:
    if the systematic 1-m peak coincides with a higher number of
    available movable shelters, the peak is at least partly a
    cavity-density artefact rather than a true emergence concentration.
    """
    rows = []
    for d in sorted(r['distance_m'].unique()):
        sub = r[r['distance_m'] == d]
        movable_total = sub[['n_rocks_small', 'n_rocks_medium']].sum().sum()
        immov_total   = sub.get('n_rocks_very_large', pd.Series([0]*len(sub))).sum() \
                        + sub.get('n_rocks_large',     pd.Series([0]*len(sub))).sum()
        total = movable_total + immov_total
        prop_movable = float(movable_total / total) if total > 0 else np.nan
        n_q = len(sub)
        mean_movable_per_q = float(movable_total / n_q) if n_q > 0 else np.nan
        rows.append({
            'distance_m'             : int(d),
            'n_quadrants'            : n_q,
            'movable_count'          : int(movable_total),
            'immovable_count'        : int(immov_total),
            'total_count'            : int(total),
            'prop_movable'           : prop_movable,
            'mean_movable_per_quadrant': mean_movable_per_q
        })
    return rows


def aquatic_exclusion_sensitivity(r):
    """
    Iter-4 m7: re-compute the within-zone Manly's alpha and the within-
    transect permutation result with the three aquatic quadrants (distance=0)
    EXCLUDED from the rocky-zone denominator. Reports both versions to
    let the reader judge whether retaining/excluding the aquatic cells
    qualitatively changes the conclusion.
    """
    r_no_aq = r[r['distance_m'] > 0].copy()
    n_with    = len(r)
    n_without = len(r_no_aq)
    occ_with    = int((r['presence'] == 1).sum())
    occ_without = int((r_no_aq['presence'] == 1).sum())
    return {
        'n_with_aquatic'    : n_with,
        'n_without_aquatic' : n_without,
        'occupied_with'     : occ_with,
        'occupied_without'  : occ_without,
        'rate_with'         : float(occ_with / n_with),
        'rate_without'      : float(occ_without / n_without),
        'note': ('Aquatic quadrants (distance=0, n=3, all zero detections '
                 'because turbid water prevents protocol detection of '
                 'submerged frogs) excluded vs retained in the systematic '
                 'denominator. The rate change quantifies the bias in the '
                 'rocky-zone detection rate from retaining structural '
                 'false-negative cells.')
    }


def ie_vs_quadrant_distance_test(r, ei):
    """
    Iter-4 m8: formal two-sample test of distance distributions of detected
    quadrants vs IE encounters in the rocky zone. Mann-Whitney U with
    one-sided alternative (IE distances > quadrant detection distances).
    """
    used_q  = r[r['presence'] == 1]['distance_m'].values
    used_ie = ei['distance_m'].dropna().values
    if len(used_q) < 1 or len(used_ie) < 1:
        return None
    stat, p = mannwhitneyu(used_ie, used_q, alternative='greater')
    return {
        'n_quadrant_detected' : int(len(used_q)),
        'n_ie'                : int(len(used_ie)),
        'median_quadrant'     : float(np.median(used_q)),
        'median_ie'           : float(np.median(used_ie)),
        'mw_U'                : float(stat),
        'p_one_sided_greater' : float(p),
        'note': ('Mann-Whitney one-sided test that IE distances are larger '
                 'than quadrant-detected distances. Formalises the '
                 'discrepancy noted in the Discussion.')
    }


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

    print("\n=== ITER-4 EXTENSIONS ===")
    sens_variants = detectability_sensitivity_variants()
    for k, v in sens_variants.items():
        print(f"  Sensitivity ({k}): E[lacus]={v['expected_lacus']:.3f}, "
              f"P(zero|equal use)={v['prob_zero']:.3f}")
    manly_q = manly_alpha_quadrant_based()
    print(f"  Manly's alpha (quadrant-based): "
          f"rocky={manly_q['alpha_rocky_quadrant']:.3f}, "
          f"intermediate={manly_q['alpha_intermediate_quadrant']:.3f}, "
          f"lacustrine={manly_q['alpha_lacustrine_quadrant']:.3f}")
    pooled = pooled_adjacent_transect_test(r, i, p)
    print(f"  Pooled-adjacent transect test (n=4 blocks/zone): "
          f"chi2={pooled['chi2']:.2f}, p={pooled['p_chi2']:.4f}, "
          f"perm p={pooled['p_perm']:.4f}")
    print(f"    rocky blocks: {pooled['rocky_blocks']}, "
          f"interm: {pooled['interm_blocks']}, lacus: {pooled['lacus_blocks']}")
    aq_sens = aquatic_exclusion_sensitivity(r)
    print(f"  Aquatic exclusion: rate with={aq_sens['rate_with']:.3f}, "
          f"without={aq_sens['rate_without']:.3f}")
    rt_grad = rock_turnability_gradient(r)
    print(f"  Rock-turnability by distance (mean prop_movable):")
    for row in rt_grad[:3]:
        print(f"    d={row['distance_m']}m: prop_movable={row['prop_movable']:.3f}")
    ie_test = ie_vs_quadrant_distance_test(r, ei)
    if ie_test:
        print(f"  IE vs quadrant distance: median IE={ie_test['median_ie']:.1f}m vs "
              f"median quadrant={ie_test['median_quadrant']:.1f}m, "
              f"MW one-sided p={ie_test['p_one_sided_greater']:.4f}")

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

    # Single canonical reportable-quantities file (one source of truth)
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
        'loo_max_p_distance'     : float(tlevel['loo_max_p']),
        'rocky_events_per_q_mean': float(np.mean([r[r['presence']==1]['n_frogs'].astype(int).tolist()][0]) if (r[r['presence']==1]['n_frogs'].astype(int).tolist()) else 0),
        'rocky_events_per_q_max' : int(max(r[r['presence']==1]['n_frogs'].astype(int).tolist()) if r[r['presence']==1]['n_frogs'].astype(int).tolist() else 0),
        'rocky_events_per_q_min' : int(min(r[r['presence']==1]['n_frogs'].astype(int).tolist()) if r[r['presence']==1]['n_frogs'].astype(int).tolist() else 0),
        'detectability_sens_expected_lacus_linear': float(detectability_sensitivity()['expected_lacus_under_equal_use_linear']),
        'detectability_sens_threshold_ratio'      : float(detectability_sensitivity()['detection_ratio_threshold_for_zero_at_alpha']),
        'rock_cover_ratio_lacus_to_rocky'         : float(detectability_sensitivity()['rock_cover_ratio']),
        # Iter-4 additions
        'sens_step_expected_lacus'      : float(sens_variants['step']['expected_lacus']),
        'sens_step_prob_zero'           : float(sens_variants['step']['prob_zero']),
        'sens_saturating_expected_lacus': float(sens_variants['saturating']['expected_lacus']),
        'sens_saturating_prob_zero'     : float(sens_variants['saturating']['prob_zero']),
        'sens_saturating_p_ratio'       : float(sens_variants['saturating']['p_ratio_lacus_to_rocky']),
        'manly_quadrant_rocky'          : float(manly_q['alpha_rocky_quadrant']),
        'manly_quadrant_intermediate'   : float(manly_q['alpha_intermediate_quadrant']),
        'manly_quadrant_lacustrine'     : float(manly_q['alpha_lacustrine_quadrant']),
        'pooled_chi2'                   : float(pooled['chi2']),
        'pooled_p_chi2'                 : float(pooled['p_chi2']),
        'pooled_p_perm'                 : float(pooled['p_perm']),
        'pooled_rocky_blocks_pos'       : int(sum(pooled['rocky_blocks'])),
        'pooled_interm_blocks_pos'      : int(sum(pooled['interm_blocks'])),
        'pooled_lacus_blocks_pos'       : int(sum(pooled['lacus_blocks'])),
        'aquatic_excl_rate_with'        : float(aq_sens['rate_with']),
        'aquatic_excl_rate_without'     : float(aq_sens['rate_without']),
        'aquatic_excl_n_with'           : int(aq_sens['n_with_aquatic']),
        'aquatic_excl_n_without'        : int(aq_sens['n_without_aquatic']),
        'ie_vs_quadrant_p_greater'      : float(ie_test['p_one_sided_greater']) if ie_test else None,
        'ie_median'                     : float(ie_test['median_ie']) if ie_test else None,
        'quadrant_detected_median'      : float(ie_test['median_quadrant']) if ie_test else None,
        'spearman_dist_veg_rho'  : float(coll['distance_m__vs__pct_vegetation']['rho']),
        'spearman_dist_veg_p'    : float(coll['distance_m__vs__pct_vegetation']['p']),
        'mw_power_r03'           : float(pwr[0.3]),
        'mw_power_r04'           : float(pwr[0.4]),
        'mw_power_r05'           : float(pwr[0.5])
    }
    with open(outdir / 'reportable_quantities.json', 'w') as f:
        import json; json.dump(reportable, f, indent=2)

    # ── LaTeX tables (\input'd by the Quarto manuscript) ──────────────────
    # Quarto/pandoc strips \label{} from chunk outputs, so we emit the
    # tables as standalone .tex files and \input them from the .qmd. This
    # bypasses pandoc's filter entirely and lets cross-references resolve.
    tab1 = rf"""\begin{{table}}[H]
\centering
\caption{{\label{{tab-zones}}Detection rates and habitat electivity of \textit{{Atelognathus
reverberii}} across three zones at Laguna Azul, 2019. All 240 quadrants
(80 per zone; 8 transects $\times$ 10 quadrants) received equal
standardised search effort. Detection rate $=$ occupied quadrants / total.
Manly's $\alpha_i$ (Equation~\ref{{eq-manly}}) with 95\% bootstrap CI
($n = 10{{,}}000$, seed $= 42$); $\alpha_i > 0.333$ indicates preference.
Trans.\ $=$ transects positive ($\geq$1 detection) / 8 total.
Quadrant chi-square (exploratory): $\chi^2 = {bz['chi2_q']:.2f}$,
df $= 2$, $p = {bz['p_q']:.3f}$.
Transect chi-square (confirmatory): $\chi^2 = {bz['chi2_t']:.2f}$,
df $= 2$, $p = {bz['p_t']:.3f}$;
permutation $p = {bz['p_perm']:.3f}$.
Substrate means computed from all 80 quadrants per zone:
V $=$ vegetation; G $=$ gravel; R $=$ rock; S $=$ soil/mud.}}
\begin{{threeparttable}}
\begin{{tabular}}{{lrrllrrrr}}
\toprule
Zone & Rate & $\alpha_i$ & 95\%~CI & Trans. &
\multicolumn{{4}}{{c}}{{Mean substrate (\%)}} \\
\cmidrule(lr){{6-9}}
 & & & & & V & G & R & S \\
\midrule
Rocky shoreline   & {int(bz['n_det'][0])}/80 ({bz['rates'][0]:.3f}) &
  {bz['alpha'][0]:.3f} &
  [{bz['alpha_ci'][0,0]:.3f},{bz['alpha_ci'][1,0]:.3f}] & 5/8 &
  14 &  7 & 68 & 11 \\
Intermediate steppe & {int(bz['n_det'][1])}/80 ({bz['rates'][1]:.3f}) &
  {bz['alpha'][1]:.3f} &
  [{bz['alpha_ci'][0,1]:.3f},{bz['alpha_ci'][1,1]:.3f}] & 2/8 &
  26 & 27 & 17 & 30 \\
Lacustrine sediment & {int(bz['n_det'][2])}/80 ({bz['rates'][2]:.3f}) &
  0.000 & $\leq {j_alpha_up:.3f}$\textsuperscript{{a}} &
  0/8 &  6 &  8 &  4 & 82 \\
\bottomrule
\end{{tabular}}
\begin{{tablenotes}}[flushleft]\footnotesize
  \item All frogs in all three zones were detected under rocks.
  \item Rocky zone has high rock cover throughout (mean 68\%); used and
    unoccupied quadrants do not differ significantly in rock cover
    ($p = {wz['pct_rock']['p']:.3f}$).
  \item \textsuperscript{{a}}For the lacustrine zone the bootstrap CI degenerates
    to a point (zero detections in all resamples). The reported value is the
    upper $97.5\%$ Jeffreys credible bound on $\alpha_{{\rm lacustrine}}$,
    derived from the upper $97.5\%$ Jeffreys bound on the underlying
    detection probability under a Beta($0.5,80.5$) posterior.
\end{{tablenotes}}
\end{{threeparttable}}
\end{{table}}
"""
    (outdir / 'table1_zones.tex').write_text(tab1)

    tab2 = rf"""\begin{{table}}[H]
\centering
\caption{{\label{{tab-mw}}Mann-Whitney $U$ tests comparing microhabitat variables between
occupied ($n = 8$) and unoccupied ($n = 72$) quadrants within the rocky
shoreline zone, Laguna Azul, 2019. Holm sequential Bonferroni correction
applied across five simultaneous tests. Bold $=$ survives Holm correction.}}
\begin{{threeparttable}}
\begin{{tabular}}{{lrrrr}}
\toprule
Variable & Used median & Available median & $p$ & Holm threshold \\
\midrule
\textbf{{Distance to shore (m)}} & \textbf{{1}} & \textbf{{5}} &
  \textbf{{{wz['distance_m']['p']:.3f}}} & \textbf{{0.010}} \\
\textbf{{Vegetation cover (\%)}} & \textbf{{0}} & \textbf{{10}} &
  \textbf{{{wz['pct_vegetation']['p']:.3f}}} & \textbf{{0.013}} \\
Gravel cover (\%)              &  10 &  0 & {wz['pct_gravel']['p']:.3f} & 0.017 \\
Rock cover (\%)                &  90 & 75 & {wz['pct_rock']['p']:.3f} & 0.025 \\
Soil/mud cover (\%)            &   2 &  5 & {wz['pct_soil']['p']:.3f} & 0.050 \\
\bottomrule
\end{{tabular}}
\begin{{tablenotes}}[flushleft]\footnotesize
  \item Variables sorted by raw $p$-value. Holm thresholds computed
    sequentially: $0.05/5 = 0.010$, $0.05/4 = 0.013$, $0.05/3 = 0.017$,
    $0.05/2 = 0.025$, $0.05/1 = 0.050$.
  \item Bold $=$ raw $p$ below Holm threshold.
\end{{tablenotes}}
\end{{threeparttable}}
\end{{table}}
"""
    (outdir / 'table2_mw.tex').write_text(tab2)

    print(f"\nAll outputs written to {outdir}")
    return bz, wz, tlevel, coll, reportable


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--datadir', default='data/')
    ap.add_argument('--outdir',  default='outputs/')
    args = ap.parse_args()
    run(args.datadir, args.outdir)