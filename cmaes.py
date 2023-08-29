import scipy.special
import numpy as np
import math
import random

try:
    import scipy.special
except ImportError:
    scipy = None
try:
    import numpy as np
except ImportError:
    np = None

def cmaes(fun, x0, sigma, g_max, trace=False):
    """CMA-ES optimization

        Parameters
        ----------
        fun : callable
              a target function
        x0 : tuple
              the initial point
        sigma : double
              initial variance
        g_max : int
              maximum generation
        trace : bool
              return a trace of the algorithm (default: False)

        Return
        ----------
        xmin : tuple"""

    def cumulation(c, A, B):
        alpha = 1 - c
        beta = math.sqrt(c * (2 - c) * mueff)
        return [alpha * a + beta * b for a, b in zip(A, B)]

    def wsum(A):
        return [
            math.fsum(w * a[i] for w, a in zip(weights, A)) for i in range(N)
        ]

    if scipy == None:
        raise ModuleNotFoundError("cmaes needs scipy")
    if np == None:
        raise ModuleNotFoundError("cmaes needs nump")
    xmean, N = x0[:], len(x0)
    lambd = 4 + int(3 * math.log(N))
    mu = lambd // 2
    weights = [math.log((lambd + 1) / 2) - math.log(i + 1) for i in range(mu)]
    weights = [e / math.fsum(weights) for e in weights]
    mueff = 1 / math.fsum(e**2 for e in weights)
    cc = (4 + mueff / N) / (N + 4 + 2 * mueff / N)
    cs = (mueff + 2) / (N + mueff + 5)
    c1 = 2 / ((N + 1.3)**2 + mueff)
    cmu = min(1 - c1, 2 * (mueff - 2 + 1 / mueff) / ((N + 2)**2 + mueff))
    damps = 1 + 2 * max(0, math.sqrt((mueff - 1) / (N + 1)) - 1) + cs
    chiN = math.sqrt(2) * math.gamma((N + 1) / 2) / math.gamma(N / 2)
    ps, pc, C = [0] * N, [0] * N, np.identity(N)
    Trace = []
    for gen in range(1, g_max + 1):
        sqrtC = np.real(scipy.linalg.sqrtm(C))
        x0 = [[random.gauss(0, 1) for d in range(N)] for i in range(lambd)]
        x1 = [sqrtC @ e for e in x0]
        xs = [xmean + sigma * e for e in x1]
        ys = [fun(e) for e in xs]
        ys, x0, x1, xs = zip(*sorted(zip(ys, x0, x1, xs)))
        xmean = wsum(xs)
        ps = cumulation(cs, ps, wsum(x0))
        pssq = math.fsum(e**2 for e in ps)
        sigma *= math.exp(cs / damps * (math.sqrt(pssq) / chiN - 1))
        Cmu = sum(w * np.outer(d, d) for w, d in zip(weights, x1))
        if (N + 1) * pssq < 2 * N * (N + 3) * (1 - (1 - cs)**(2 * gen)):
            pc = cumulation(cc, pc, wsum(x1))
            C1 = np.outer(pc, pc)
            C = (1 - c1 - cmu) * C + c1 * C1 + cmu * Cmu
        else:
            pc = [(1 - cc) * e for e in pc]
            C1 = np.outer(pc, pc)
            C = (1 - c1 - cmu) * C + c1 * (C1 + cc * (2 - cc) * C) + cmu * Cmu
        if trace:
            Trace.append(
                (gen * lambd, ys[0], xs[0], sigma, C, ps, pc, Cmu, C1, xmean))
    return Trace if trace else xmean
