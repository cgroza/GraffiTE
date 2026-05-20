"""
HERV-K SV polymorphism classifier — v3

Five hypotheses:
    H_C : null     <-> solo-LTR        (canonical |SVLEN| ~  968 bp)
    H_T : truncated proviral           (1500 <= |SVLEN| <= 8000, one LTR + partial INT)
    H_B : solo-LTR <-> proviral        (canonical |SVLEN| ~ 8504 bp)
    H_A : null     <-> proviral        (canonical |SVLEN| ~ 9472 bp)
    H_X : non-transposition / other    (flat background)

Family whitelists:
    LTR bucket: LTR5_Hs, LTR5A, LTR5B   (HML-2 LTRs only)
    INT bucket: HERVK-int               (HML-2 internal only)

Intended for use inside GraffiTE when the --human flag is set; produces
posteriors that downstream code can write into the trusted/human-filtered
VCF/TSV outputs.
"""

import numpy as np
import pandas as pd

# ---- Reference architecture ------------------------------------------------
LTR_LEN, INT_LEN = 968, 7536
SOLO_PROV  = LTR_LEN + INT_LEN          # 8504
NULL_PROV  = 2 * LTR_LEN + INT_LEN      # 9472

# Truncated proviral plausibility window
T_MIN, T_MAX = 1500, 8000

# Expected (s*, lam*, nu*) for the structured Gaussian hypotheses ------------
EXPECTED = {
    'C': {'s': LTR_LEN,    'lam': LTR_LEN,   'nu': 0       },
    'B': {'s': SOLO_PROV,  'lam': LTR_LEN,   'nu': INT_LEN },
    'A': {'s': NULL_PROV,  'lam': 2*LTR_LEN, 'nu': INT_LEN },
}

# ---- Family whitelists -----------------------------------------------------
LTR_FAMILY = {'LTR5_Hs', 'LTR5A', 'LTR5B'}
INT_FAMILY = {'HERVK-int'}

# ---- Defaults (see plan §4 for justification) ------------------------------
DEFAULT_SIGMAS = {
    's_C':  30.0,    # solo-LTRs are very length-uniform
    's_B':  800.0,   # allow ~10% INT truncation in proviral
    's_A':  800.0,
    'lam':  300.0,   # absorbs (x)-merged LTR/INT calls
    'nu':   1000.0,
    't':    200.0,   # SV must be mostly HERV-K
}

# Literature-informed defaults (plan §5):
#   solo-LTR dominates polymorphism; full-length proviral insertions are
#   rare; truncated proviral is intermediate; H_X reserves enough mass to
#   reject SVs that merely contain HERV-K fragments.
DEFAULT_PRIORS = {'C': 0.55, 'T': 0.08, 'B': 0.05, 'A': 0.02, 'X': 0.30}

# Background flat density for H_X
S_RANGE, LAM_RANGE, NU_RANGE = 30000.0, 5000.0, 30000.0
LOG_BACKGROUND = -np.log(S_RANGE * LAM_RANGE * NU_RANGE)

# H_T flat density on s within window (lam, nu still Gaussian-constrained)
LOG_T_S_DENSITY = -np.log(T_MAX - T_MIN)


# ---- Annotation parsing ----------------------------------------------------
def parse_hits(match_lengths, repeat_ids):
    """Return (lambda, nu) — bp matching LTR family, bp matching INT family.

    Handles the GraffiTE '(x)' suffix (RepeatMasker processRepeat merge tag)
    and comma-separated multi-hit annotations.
    """
    if pd.isna(match_lengths) or pd.isna(repeat_ids):
        return 0.0, 0.0
    lens = [float(x) for x in str(match_lengths).split(',')]
    ids  = [r.strip().replace('(x)', '') for r in str(repeat_ids).split(',')]
    lam = sum(L for L, r in zip(lens, ids) if r in LTR_FAMILY)
    nu  = sum(L for L, r in zip(lens, ids) if r in INT_FAMILY)
    return lam, nu


# ---- Likelihoods -----------------------------------------------------------
def log_likelihood_gaussian(s, lam, nu, hyp, sigmas):
    """Gaussian log-likelihood for H_C, H_B, H_A (up to additive constant)."""
    e, sig_s, t = EXPECTED[hyp], sigmas[f's_{hyp}'], lam + nu
    return (
        -0.5 * ((s   - e['s'])   / sig_s        ) ** 2
        -0.5 * ((lam - e['lam']) / sigmas['lam']) ** 2
        -0.5 * ((nu  - e['nu'])  / sigmas['nu'] ) ** 2
        -0.5 * ((s   - t)        / sigmas['t']  ) ** 2
    )


def log_likelihood_truncated(s, lam, nu, sigmas):
    """H_T: flat in s over (T_MIN, T_MAX); Gaussian on lam ≈ LTR_LEN; coverage."""
    if s < T_MIN or s > T_MAX:
        return -1e6  # effectively excluded
    t = lam + nu
    return (
        LOG_T_S_DENSITY
        - 0.5 * ((lam - LTR_LEN) / sigmas['lam']) ** 2
        - 0.5 * ((s   - t)       / sigmas['t']  ) ** 2
    )


# ---- Classification --------------------------------------------------------
def classify_row(s, lam, nu, priors=DEFAULT_PRIORS, sigmas=DEFAULT_SIGMAS):
    """Return posterior dict over {'C','T','B','A','X'} for one SV."""
    lp = {k: log_likelihood_gaussian(s, lam, nu, k, sigmas) + np.log(priors[k])
          for k in 'ABC'}
    lp['T'] = log_likelihood_truncated(s, lam, nu, sigmas) + np.log(priors['T'])
    lp['X'] = LOG_BACKGROUND                               + np.log(priors['X'])
    m = max(lp.values())
    e = {k: np.exp(v - m) for k, v in lp.items()}
    Z = sum(e.values())
    return {k: v / Z for k, v in e.items()}


def classify_table(df, priors=DEFAULT_PRIORS, sigmas=DEFAULT_SIGMAS):
    """Add lambda_LTR, nu_INT, P_H_*, MAP_class, MAP_posterior to df."""
    out = df.copy()
    rows = []
    for _, row in df.iterrows():
        lam, nu = parse_hits(row['match_lengths'], row['repeat_ids'])
        s = abs(row['SVLEN'])
        post = classify_row(s, lam, nu, priors, sigmas)
        rows.append((lam, nu,
                     post['C'], post['T'], post['B'], post['A'], post['X']))
    cols = ['lambda_LTR','nu_INT','P_H_C','P_H_T','P_H_B','P_H_A','P_H_X']
    out[cols] = pd.DataFrame(rows, index=df.index)
    pcols = ['P_H_C','P_H_T','P_H_B','P_H_A','P_H_X']
    out['MAP_class']     = out[pcols].idxmax(axis=1).str.replace('P_H_', '')
    out['MAP_posterior'] = out[pcols].max(axis=1)
    return out


# ---- Allelic-state interpretation -----------------------------------------
ALLELE_INTERPRETATION = {
    'C': {'ref': 'null',  'alt': 'solo' },
    'T': {'ref': 'null',  'alt': 'truncated_prov'},
    'B': {'ref': 'solo',  'alt': 'prov' },
    'A': {'ref': 'null',  'alt': 'prov' },
    'X': {'ref': '?',     'alt': '?'    },
}

def interpret_genotype(map_class, gt_value, svtype):
    """Map a 0/1 genotype call to a HERV-K allelic state.

    Convention: SVTYPE=INS means alt allele is the longer (inserted) state;
    SVTYPE=DEL means alt allele is the shorter (deleted) state. The
    ALLELE_INTERPRETATION table is written from the perspective of "alt =
    inserted/encoded state" — for SVTYPE=DEL we swap ref<->alt so that
    gt=1 still corresponds to the shorter-allele state on the haplotype.
    """
    if pd.isna(gt_value) or gt_value in ('.', './.'):
        return '?'
    interp = ALLELE_INTERPRETATION[map_class]
    if svtype == 'DEL':
        interp = {'ref': interp['alt'], 'alt': interp['ref']}
    return interp['alt' if str(gt_value).strip() in ('1','1|1','1|0','0|1') else 'ref']


# ---- Script entry point ----------------------------------------------------
if __name__ == '__main__':
    import sys
    src = sys.argv[1] if len(sys.argv) > 1 \
          else '/mnt/user-data/uploads/HERVK_annot_test.tsv'
    df = pd.read_csv(src, sep='\t')
    out = classify_table(df)
    keep = ['CHROM','POS','SVTYPE','SVLEN','repeat_ids',
            'lambda_LTR','nu_INT',
            'P_H_C','P_H_T','P_H_B','P_H_A','P_H_X',
            'MAP_class','MAP_posterior']
    pd.set_option('display.width', 220)
    print(out[keep].to_string(index=False, float_format=lambda x: f'{x:.3f}'))
