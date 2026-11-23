#!/usr/bin/env python3
"""Independent re-check of the terrain boundary-normal analysis.

Re-derives the stencils and the analytic surface from scratch (no shared code
with src/terrain.c), verifies:
  1. truncation-error bound constants on pure cubics (M3 h^2/6, M3 h^2/3,
     and the old M2 h/2 leading term);
  2. normal angles at the four edges against the propagated angle bound
     alpha = asin(eps/(1 - s eps));
  3. convergence orders O(h^2) (new) vs O(h) (old).
Run: python3 -I verify_boundary_normals.py
"""
import math

W = H = 16
M3X, M3Y = 0.06, 0.12  # fxxx, fyyy of the analytic surface


def f(x, y):
    return 0.01 * x**3 + 0.02 * y**3 + 0.3 * x**2 + 0.4 * y**2 + 0.1 * x * y + 1.0


def grad_f(x, y):
    return (0.03 * x**2 + 0.6 * x + 0.1 * y,
            0.06 * y**2 + 0.8 * y + 0.1 * x)


def build(h):
    return [[f(x * h, y * h) for y in range(H)] for x in range(W)]


def grad_new(z, h, x, y):
    """centered interior, 3-point one-sided 2nd-order at edges."""
    if 0 < x < W - 1:
        gx = (z[x + 1][y] - z[x - 1][y]) / (2 * h)
    elif x == 0:
        gx = (-3 * z[0][y] + 4 * z[1][y] - z[2][y]) / (2 * h)
    else:
        gx = (3 * z[W - 1][y] - 4 * z[W - 2][y] + z[W - 3][y]) / (2 * h)
    if 0 < y < H - 1:
        gy = (z[x][y + 1] - z[x][y - 1]) / (2 * h)
    elif y == 0:
        gy = (-3 * z[x][0] + 4 * z[x][1] - z[x][2]) / (2 * h)
    else:
        gy = (3 * z[x][H - 1] - 4 * z[x][H - 2] + z[x][H - 3]) / (2 * h)
    return gx, gy


def grad_old(z, h, x, y):
    """clamped stencil: degrades to 1st-order one-sided difference."""
    xm, xp = (x - 1 if x > 0 else x), (x + 1 if x < W - 1 else x)
    ym, yp = (y - 1 if y > 0 else y), (y + 1 if y < H - 1 else y)
    return ((z[xp][y] - z[xm][y]) / ((xp - xm) * h),
            (z[x][yp] - z[x][ym]) / ((yp - ym) * h))


def normal(gx, gy):
    inv = 1.0 / math.sqrt(gx * gx + gy * gy + 1.0)
    return (-gx * inv, -gy * inv, inv)


def angle(a, b):
    d = max(-1.0, min(1.0, sum(u * v for u, v in zip(a, b))))
    return math.degrees(math.acos(d))


def alpha_deg(eps, s):
    v = eps / (1.0 - s * eps)
    return math.degrees(math.asin(v)) if v < 1.0 else 180.0


def main():
    ok = True

    # --- 1. bound constants on exact cubics (float64: rounding negligible) ---
    h = 0.5
    # g(t) = t^3: g'=3t^2, g''=6t (M2 over [0,2h]=6h), g'''=6 (M3)
    g = [((i * h) ** 3) for i in range(4)]
    centered = (g[2] - g[0]) / (2 * h)            # at t=h
    onesided = (-3 * g[0] + 4 * g[1] - g[2]) / (2 * h)  # at t=0
    forward = (g[1] - g[0]) / h                   # old clamp stencil at t=0
    e_centered = abs(centered - 3 * h**2)
    e_onesided = abs(onesided - 0.0)
    e_forward = abs(forward - 0.0)
    b_centered = 6 * h**2 / 6
    b_onesided = 6 * h**2 / 3
    b_forward = 6 * h * h / 2 + 6 * h**2 / 6      # M2 h/2 + M3 h^2/6
    print("== truncation bound constants (g(t)=t^3, h=%.2f) ==" % h)
    print(f"  centered  err={e_centered:.3e} bound M3 h^2/6={b_centered:.3e} "
          f"hold={e_centered <= b_centered + 1e-15}")
    print(f"  one-sided err={e_onesided:.3e} bound M3 h^2/3={b_onesided:.3e} "
          f"hold={e_onesided <= b_onesided + 1e-15}")
    print(f"  forward   err={e_forward:.3e} bound M2h/2+M3h^2/6={b_forward:.3e} "
          f"hold={e_forward <= b_forward + 1e-15}")
    ok &= e_centered <= b_centered + 1e-15
    ok &= e_onesided <= b_onesided + 1e-15
    ok &= e_forward <= b_forward + 1e-15

    # quadratic q(t)=3t^2+2t: f'''=0, so the M2 h/2 leading term is exactly tight
    q = [3 * (i * h) ** 2 + 2 * (i * h) for i in range(2)]
    e_fwd_q = abs((q[1] - q[0]) / h - 2.0)       # q'(0)=2, M2=6
    b_fwd_q = 6 * h / 2
    print(f"  forward(quad) err={e_fwd_q:.3e} M2 h/2={b_fwd_q:.3e} "
          f"tight={abs(e_fwd_q - b_fwd_q) < 1e-15}")
    ok &= abs(e_fwd_q - b_fwd_q) < 1e-15

    # --- 2. four-edge normal angle check -------------------------------------
    z = build(h)
    ebx, eby = h**2 / 3 * M3X, h**2 / 3 * M3Y
    eps_bnd = math.hypot(ebx, eby)
    edges = [("west", 0, 7), ("east", W - 1, 7),
             ("south", 7, 0), ("north", 7, H - 1)]
    print("\n== four-edge normal angles (h=%.2f, eps_bnd=%.6f) ==" %
          (h, eps_bnd))
    for name, x, y in edges:
        ax, ay = grad_f(x * h, y * h)
        s = math.hypot(ax, ay)
        bound = alpha_deg(eps_bnd, s)
        n_ref = normal(ax, ay)
        a_old = angle(n_ref, normal(*grad_old(z, h, x, y)))
        gxn, gyn = grad_new(z, h, x, y)
        a_new = angle(n_ref, normal(gxn, gyn))
        gerr = math.hypot(gxn - ax, gyn - ay)
        new_ok = a_new <= bound + 1e-12
        old_bad = a_old > bound
        ok &= new_ok and old_bad
        print(f"  {name:<5} bound={bound:7.4f}  old={a_old:7.4f}  "
              f"new={a_new:7.4f}  |gerr|={gerr:.2e}  new_inside={new_ok} "
              f"old_exceeds={old_bad}")

    # --- 3. interior centered bound ------------------------------------------
    x = y = 7
    ax, ay = grad_f(x * h, y * h)
    gxn, gyn = grad_new(z, h, x, y)
    eps_in = math.hypot(h**2 / 6 * M3X, h**2 / 6 * M3Y)
    ierr = math.hypot(gxn - ax, gyn - ay)
    in_ok = ierr <= eps_in + 1e-15
    ok &= in_ok
    print(f"\n== interior (7,7): |gerr|={ierr:.3e} eps_in={eps_in:.3e} "
          f"hold={in_ok} ==")

    # --- 4. convergence orders -----------------------------------------------
    print("\n== convergence at west edge (0,7) ==")
    en, eo = [], []
    for hh in (1.0, 0.5, 0.25):
        zz = build(hh)
        ax, ay = grad_f(0.0, 7 * hh)
        en.append(math.hypot(*[a - b for a, b in
                               zip(grad_new(zz, hh, 0, 7), (ax, ay))]))
        eo.append(math.hypot(*[a - b for a, b in
                               zip(grad_old(zz, hh, 0, 7), (ax, ay))]))
    rn = [en[i] / en[i + 1] for i in range(2)]
    ro = [eo[i] / eo[i + 1] for i in range(2)]
    print(f"  new ratios={[round(v, 3) for v in rn]} (expect ~4, O(h^2))")
    print(f"  old ratios={[round(v, 3) for v in ro]} (expect ~2, O(h))")
    ok &= all(3.9 < v < 4.1 for v in rn) and all(1.9 < v < 2.1 for v in ro)

    print("\npython_check_all_ok=%s" % ok)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
