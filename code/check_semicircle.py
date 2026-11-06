"""I0082 independent check: mu = semicircle of variance v (R-transform v*g), kappa = p delta_0 + q U[0,1].
G_sigma(z) = G_kappa(z - v G_sigma(z)) for sigma = mu boxplus kappa. Density f = -Im G / pi on the real line.
A unimodal law cannot have f' < 0 at x1 < x2 with f' > 0 at x2. Also scan f on a grid for an interior local minimum.
Also mu = symmetric Bernoulli(+-a): solve the two subordination equations by iteration on the upper half-plane."""
import numpy as np
p, q = 0.01, 0.99
OUT = __file__.replace(".py", ".out")


def Gk(z):
    return p / z + q * (np.log(z) - np.log(z - 1))


def G_semi(x, v, eta=1e-13):
    z = x + 1j * eta
    G = Gk(z + 1j * 1e-3)
    for _ in range(400):
        Gn = Gk(z - v * G)
        if abs(Gn - G) < 1e-17:
            break
        G = Gn
    return Gn


def G_bern(x, a, eta=1e-13):
    # mu = (delta_a + delta_-a)/2, F_mu(w) = (w^2 - a^2)/w. Subordination: F_mu(w1) = F_k(w2) = F, w1 + w2 = z + F.
    # iterate w2 -> z + h_mu(z + h_k(w2)), h = F - id.
    z = x + 1j * eta
    Fk = lambda w: 1 / Gk(w)
    hmu = lambda w: -a * a / w
    hk = lambda w: Fk(w) - w
    w2 = z + 1j * 1e-3
    for _ in range(2000):
        w2n = z + hmu(z + hk(w2))
        if abs(w2n - w2) < 1e-16:
            break
        w2 = w2n
    return Gk(w2n)


def dens(Gfun, x, par):
    return -Gfun(x, par).imag / np.pi


lines = []
for name, Gfun, par in (("semicircle var 1e-6", G_semi, 1e-6), ("semicircle var 1e-4", G_semi, 1e-4),
                        ("Bernoulli +-1e-3", G_bern, 1e-3)):
    h = 1e-5
    d = [(dens(Gfun, x + h, par) - dens(Gfun, x - h, par)) / (2 * h) for x in (0.01, 0.25)]
    xs = np.linspace(0.004, 0.6, 3000)
    f = np.array([dens(Gfun, x, par) for x in xs])
    i = int(np.argmin(f[: len(f) // 2])) if True else 0
    imin = np.argmin(f)
    loc = [k for k in range(1, len(f) - 1) if f[k] < f[k - 1] and f[k] < f[k + 1]]
    lines.append("%s: f'(0.01)=%.6e f'(0.25)=%.6e  interior local minima at x=%s (f=%s)  mass check sum f dx over (0.004,0.6)=%.4f" % (
        name, d[0], d[1], [round(xs[k], 4) for k in loc], [round(f[k], 6) for k in loc], np.trapezoid(f, xs)))
for l in lines:
    print(l)
open(OUT, "w").write("\n".join(lines) + "\n")
