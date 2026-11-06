"""I0082: exact / high-precision check of the constants in Section 3 of the paper.
p = 1/100, q = 99/100, g(x) = p/x + q log(x/(1-x)), x1 = p, x2 = 1/4, c = pi q,
D = disc centre b = g(x2) - 10 i, radius r = 99/10, w_j = g(x_j) - i c."""
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 50
p, q = Fr(1, 100), Fr(99, 100)
x1, x2 = p, Fr(1, 4)
gp = lambda x: (x - p) / (x * x * (1 - x))                          # g'
gpp = lambda x: 2 * p / x ** 3 + q * (-1 / x ** 2 + 1 / (1 - x) ** 2)  # g''
out = []
out.append("g'(x1) = %s, g''(x1) = %s, g'(x2) = %s, g''(x2) = %s" % (gp(x1), gpp(x1), gp(x2), gpp(x2)))
# check g' formula against symbolic derivative numerically
g = lambda x: mp.mpf(p.numerator) / p.denominator / x + mp.mpf(99) / 100 * mp.log(x / (1 - x))
for x in (x1, x2):
    xm = mp.mpf(x.numerator) / x.denominator
    out.append("  numeric g'(%s)=%s g''=%s" % (x, mp.nstr(mp.diff(g, xm), 12), mp.nstr(mp.diff(g, xm, 2), 12)))
c = mp.pi * mp.mpf(99) / 100
d = g(mp.mpf(1) / 100) - g(mp.mpf(1) / 4)
out.append("d = g(x1)-g(x2) = %s (closed form 24/25 - (99/100) log 33 = %s)" % (mp.nstr(d, 15), mp.nstr(mp.mpf(24) / 25 - mp.mpf(99) / 100 * mp.log(33), 15)))
r = mp.mpf(99) / 10
b = g(mp.mpf(1) / 4) - 10j
w1, w2 = g(mp.mpf(1) / 100) - 1j * c, g(mp.mpf(1) / 4) - 1j * c
out.append("|w1-b| = %s, |w2-b| = %s, r = 9.9; top of closed disc Im = %s" % (mp.nstr(abs(w1 - b), 10), mp.nstr(abs(w2 - b), 10), mp.nstr(-10 + r, 5)))
Delta = r ** 2 - (10 - c) ** 2
rho = r * abs(d) / abs(r ** 2 - (w1 - b) * mp.conj(w2 - b))
out.append("Delta = %s (> 48.5891: %s); 2r/Delta = %s (< 5/12=0.41667: %s)" % (mp.nstr(Delta, 12), Delta > 48.5891, mp.nstr(2 * r / Delta, 10), 2 * r / Delta < mp.mpf(5) / 12))
out.append("rho = %s (< 2/3: %s); bound r|d|/Delta = %s; (1-rho)/(1+rho) = %s (> 1/5)" % (mp.nstr(rho, 10), rho < mp.mpf(2) / 3, mp.nstr(r * abs(d) / Delta, 10), mp.nstr((1 - rho) / (1 + rho), 10)))
# pseudohyperbolic distance in a disc: |(w1-w2) r| / |r^2 - (w1-b) conj(w2-b)|  -- check against Mobius formula
z1, z2 = (w1 - b) / r, (w2 - b) / r
out.append("Mobius rho = %s" % mp.nstr(abs((z1 - z2) / (1 - z1 * mp.conj(z2))), 10))
# step 7: derivative bounds, with K(w1) >= t/5 and |K'(w2)| <= (5/12) t, K(w2) = t
c1 = -gpp(x1) * Fr(1, 5)                 # f'(x1) <= -(g''(x1) K(w1)) / pi, g'(x1)=0
c2 = -gpp(x2) - gp(x2) ** 2 * Fr(5, 12)  # f'(x2) >= [-g''(x2) t - g'(x2)^2 |K'|] / pi
out.append("f'(x1) pi / t <= %s (= -200000/99: %s); f'(x2) pi / t >= %s (= 704/375: %s)" % (c1, c1 == Fr(-200000, 99), c2, c2 == Fr(704, 375)))
for l in out:
    print(l)
open(__file__.replace(".py", ".out"), "w").write("\n".join(out) + "\n")
