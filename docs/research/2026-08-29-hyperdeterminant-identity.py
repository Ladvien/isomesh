#!/usr/bin/env python3
"""The symbolic check behind P-123/R-123.

Claim: `b*b - 4*a*c`, built from the coefficients in
`crates/isomesh/src/marching_cubes/trilinear.rs:200-214` and consumed at :246, is
IDENTICALLY Cayley's 2x2x2 hyperdeterminant of the eight corner values, under the
crate's own corner indexing `f[u + 2v + 4w]`.

Also checked here, because each one is a claim in the ledger:
  * it equals disc(det(A0 + lam*A1)) for the two opposite-face 2x2 corner
    matrices, for ALL THREE axis pairings   -> the mechanism behind M-206;
  * Delta(g.A) = (det g1 * det g2 * det g3)^2 * Delta(A), a GL(2)^3 relative
    invariant of SQUARE weight, so sign(Delta) -- and hence the body-saddle
    count -- is invariant under every per-axis affine reparametrisation, cell
    aspect ratio included;
  * Delta is invariant under all 48 octahedral relabellings of the cube's
    corners, and under negating the field.

Requires sympy.  Runs in a few seconds.  Exits non-zero if any check fails.
"""

import itertools
import random
import sys

import sympy as sp

F = sp.symbols("f0:8")  # corner[0..7], index = u + 2*v + 4*w
LAM = sp.symbols("lam")


def repo_discriminant(f):
    """b*b - 4*a*c, transcribed from BodySaddles::coefficients."""
    twist_lo = (f[0] + f[3]) - (f[1] + f[2])
    twist_hi = (f[4] + f[7]) - (f[5] + f[6])
    du_lo, du_hi = f[1] - f[0], f[5] - f[4]
    dv_lo, dv_hi = f[2] - f[0], f[6] - f[4]

    a = du_hi * twist_lo - du_lo * twist_hi
    b = (f[4] * twist_lo - f[0] * twist_hi) + (du_hi * dv_lo - du_lo * dv_hi)
    c = f[2] * f[4] - f[0] * f[6]
    return sp.expand(b * b - 4 * a * c)


def cayley(f):
    """Cayley's hyperdeterminant of the 2x2x2 tensor A[i][j][k] = f[i + 2j + 4k]."""
    a = lambda i, j, k: f[i + 2 * j + 4 * k]  # noqa: E731
    return sp.expand(
        a(0, 0, 0) ** 2 * a(1, 1, 1) ** 2
        + a(0, 0, 1) ** 2 * a(1, 1, 0) ** 2
        + a(0, 1, 0) ** 2 * a(1, 0, 1) ** 2
        + a(1, 0, 0) ** 2 * a(0, 1, 1) ** 2
        - 2
        * (
            a(0, 0, 0) * a(0, 0, 1) * a(1, 1, 0) * a(1, 1, 1)
            + a(0, 0, 0) * a(0, 1, 0) * a(1, 0, 1) * a(1, 1, 1)
            + a(0, 0, 0) * a(1, 0, 0) * a(0, 1, 1) * a(1, 1, 1)
            + a(0, 0, 1) * a(0, 1, 0) * a(1, 0, 1) * a(1, 1, 0)
            + a(0, 0, 1) * a(1, 0, 0) * a(0, 1, 1) * a(1, 1, 0)
            + a(0, 1, 0) * a(1, 0, 0) * a(0, 1, 1) * a(1, 0, 1)
        )
        + 4
        * (
            a(0, 0, 0) * a(0, 1, 1) * a(1, 0, 1) * a(1, 1, 0)
            + a(0, 0, 1) * a(0, 1, 0) * a(1, 0, 0) * a(1, 1, 1)
        )
    )


def act(f, g1, g2, g3):
    """(g1, g2, g3) acting on the 2x2x2 tensor."""
    out = [0] * 8
    for i, j, k in itertools.product(range(2), repeat=3):
        out[i + 2 * j + 4 * k] = sum(
            g1[i][p] * g2[j][q] * g3[k][r] * f[p + 2 * q + 4 * r]
            for p, q, r in itertools.product(range(2), repeat=3)
        )
    return out


def octahedral_permutations():
    """The 48 relabellings of the cube's corners: 3! axis permutations x 2^3 flips."""
    perms = []
    for pi in itertools.permutations(range(3)):
        for flips in itertools.product((0, 1), repeat=3):
            p = [0] * 8
            for idx in range(8):
                bits = [(idx >> t) & 1 for t in range(3)]
                nb = [bits[pi[t]] ^ flips[t] for t in range(3)]
                p[nb[0] + 2 * nb[1] + 4 * nb[2]] = idx
            perms.append(tuple(p))
    return perms


def main() -> int:
    failures = []
    disc = repo_discriminant(F)
    det = cayley(F)

    # 1. the identity itself
    if sp.simplify(disc - det) != 0:
        failures.append("repo discriminant != Cayley hyperdeterminant")
    if sp.Poly(disc, F).total_degree() != 4 or len(disc.args) != 12:
        failures.append("expected a 12-term degree-4 form")

    # 2. the pencil, on all three axis pairings
    pairings = {
        "w-slices 0123|4567": ((0, 1, 2, 3), (4, 5, 6, 7)),
        "v-slices 0145|2367": ((0, 1, 4, 5), (2, 3, 6, 7)),
        "u-slices 0246|1357": ((0, 2, 4, 6), (1, 3, 5, 7)),
    }
    for name, (s0, s1) in pairings.items():
        A = sp.Matrix(2, 2, [F[i] for i in s0])
        B = sp.Matrix(2, 2, [F[i] for i in s1])
        r, q, p = sp.Poly(sp.expand((A + LAM * B).det()), LAM).all_coeffs()
        if sp.simplify((q * q - 4 * r * p) - disc) != 0:
            failures.append(f"pencil discriminant mismatch on {name}")

    # 3. GL(2)^3 relative invariance, weight (det g1 det g2 det g3)^2
    random.seed(20260829)
    rnd = lambda: sp.Rational(random.randint(-9, 9), random.randint(1, 5))  # noqa: E731
    trials = 0
    while trials < 25:
        v = [rnd() for _ in range(8)]
        gs = [[[rnd(), rnd()], [rnd(), rnd()]] for _ in range(3)]
        dets = [sp.Matrix(g).det() for g in gs]
        if any(d == 0 for d in dets):
            continue
        trials += 1
        lhs = cayley(act(v, *gs))
        rhs = (dets[0] * dets[1] * dets[2]) ** 2 * cayley(v)
        if sp.simplify(lhs - rhs) != 0:
            failures.append("GL(2)^3 weight is not (det g1 det g2 det g3)^2")
            break

    # 4. octahedral and negation invariance
    perms = octahedral_permutations()
    assert len(set(perms)) == 48
    base = cayley(F)
    for p in perms:
        if sp.simplify(cayley([F[p[i]] for i in range(8)]) - base) != 0:
            failures.append("not invariant under an octahedral relabelling")
            break
    if sp.simplify(cayley([-x for x in F]) - base) != 0:
        failures.append("not invariant under negating the field")

    for line in (
        "identity          b*b - 4*a*c == Det_2,2,2",
        "pencil            == disc(det(A0 + lam A1)), all three axis pairings",
        "GL(2)^3           relative invariant, weight a perfect square",
        "octahedral        invariant under all 48 relabellings",
        "negation          invariant under f -> -f",
    ):
        print(f"  {line}")
    if failures:
        print("\nFAILED:")
        for f in failures:
            print(f"  - {f}")
        return 1
    print("\nall checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
