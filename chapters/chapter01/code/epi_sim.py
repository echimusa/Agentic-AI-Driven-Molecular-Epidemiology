"""
Chapter 1 modelling foundations: a stochastic model of the three layers behind a genomic epidemiology data set.
Layer 1, transmission: a branching process. Each case produces a negative-binomial number of secondary cases (mean R, dispersion k; k = infinity is Poisson) at times drawn from a gamma generation-time distribution.
Layer 2, evolution: mutations arise along each host's lineage as a Poisson process with rate lam = mu * L per day (infinite-sites model: every mutation is new). A child inherits the parent's genome as it stood at the moment of infection.
Layer 3, observation: each case is sequenced with probability rho, some days after infection.
Default parameters are anchored on Chapter 5: Ebola virus substitution rate 1.25e-3 per site per year (about 23.7 per genome-year, 18,961 sites), and an assumed mean generation time of 15 days.
"""
import numpy as np
MU_SITE_YEAR, L_GENOME = 1.25e-3, 18961; LAM_DEFAULT = MU_SITE_YEAR * L_GENOME / 365.0          # mutations per genome per day (0.065)
H = 60.0                                                                                       # horizon (days) over which a host's mutations are generated

def simulate_outbreak(rng, R=1.5, k=1.0, Tg=15.0, shape=4.0, lam=LAM_DEFAULT, n_cap=300, min_cases=60, root_genome=frozenset(), mut_counter=None, max_tries=200):
    mut_counter = mut_counter if mut_counter is not None else [0]
    for _ in range(max_tries):
        parent = [-1]; t_inf = [0.0]; queue = [0]; qi = 0
        while qi < len(queue) and len(parent) < n_cap:
            i = queue[qi]; qi += 1
            z = rng.poisson(R) if k is None else rng.negative_binomial(k, k / (k + R))
            for _ in range(z):
                if len(parent) >= n_cap: break
                g = min(rng.gamma(shape, Tg / shape), H); parent.append(i); t_inf.append(t_inf[i] + g); queue.append(len(parent) - 1)
        if len(parent) >= min_cases: break
    n = len(parent); events = []; inherited = []
    for i in range(n):
        m = rng.poisson(lam * H); tt = np.sort(rng.uniform(t_inf[i], t_inf[i] + H, m)); ids = []
        for _ in range(m): mut_counter[0] += 1; ids.append(mut_counter[0])
        events.append((tt, ids))
        if i == 0: inherited.append(set(root_genome))
        else:
            p = parent[i]; tp, ip = events[p]; inherited.append(set(inherited[p]) | {ip[j] for j in range(len(tp)) if tp[j] <= t_inf[i]})
    return dict(parent=np.array(parent), t_inf=np.array(t_inf), events=events, inherited=inherited, n=n, processed=qi)      # processed = cases whose offspring were all drawn

def sample_genome(o, i, delay):
    tt, ids = o["events"][i]; return frozenset(o["inherited"][i] | {ids[j] for j in range(len(tt)) if tt[j] <= o["t_inf"][i] + delay})

def depths(parent):
    d = np.zeros(len(parent), int)
    for i in range(1, len(parent)): d[i] = d[parent[i]] + 1
    return d

def steps_between(parent, dep, u, v):
    a, b = u, v; s = 0
    while a != b:
        if dep[a] >= dep[b]: a = parent[a]
        else: b = parent[b]
        s += 1
    return s
