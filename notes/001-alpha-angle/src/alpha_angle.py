# ITU-R S.1503-4 — Cálculo do ângulo α (alpha) mínimo ao arco GSO visível
#
# Implementação para Terra esférica (||rE|| = Re), consistente com a Recomendação.
#
# Saídas principais:
#   - |alpha|      : mínimo no arco GSO visível do ângulo entre ES→NGSO e ES→ponto no arco
#   - sinal(alpha) : conforme §D6.4.4.3 (critério geométrico via interseção com o plano XY)
#   - DeltaLong    : LongAlpha − LongNGSO, com desempate:
#         * se houver dois mínimos, escolher min |DeltaLong|
#         * se empatar com sinais opostos, escolher DeltaLong positivo
#
# Interface recomendada:
#   alpha_to_gso_arc(rE, rN, signed=True, output="deg", return_details=True)

import math
from typing import Dict, Literal, NamedTuple, Optional, Tuple, Union

import numpy as np


# ---------------------------------------------------------------------
# Constantes (S.1503, valores típicos usados na Tabela 2)
# ---------------------------------------------------------------------
RE_KM_DEFAULT: float = 6378.145
RGEO_KM_DEFAULT: float = 42164.2

AngleUnit = Literal["rad", "deg"]


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------
def _wrap_pi(x: float) -> float:
    """Wrap (rad) to (-pi, pi]."""
    return math.atan2(math.sin(x), math.cos(x))


def _lon_from_ecef(r: np.ndarray) -> float:
    r = np.asarray(r, float).reshape(3)
    return math.atan2(float(r[1]), float(r[0]))


def _to_unit(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, float).reshape(3)
    n = float(np.linalg.norm(v))
    if n == 0.0:
        raise ValueError("Cannot normalize zero vector.")
    return v / n


# ---------------------------------------------------------------------
# Cache esférico por estação (usa lat/lon geocêntricos)
# ---------------------------------------------------------------------
class ESCacheSph(NamedTuple):
    lon: float        # rad
    clon: float
    slon: float
    px: float         # Re*cos(lat)
    pz: float         # Re*sin(lat)
    Rgeo: float
    E0: float         # Rgeo^2 + Re^2
    F0: float         # -2*Rgeo*px
    thmax: float      # theta_max >= 0
    has_arc: bool


def es_cache_spherical(
    lat_es: float,
    lon_es: float,
    Re: float,
    Rgeo: float,
    eps_vis: float = 1e-15
) -> ESCacheSph:
    """Cache para Terra esférica. theta_max: cos(theta_max)=Re/(Rgeo*cos(lat_es))."""
    clat = math.cos(lat_es)
    slat = math.sin(lat_es)
    clon = math.cos(lon_es)
    slon = math.sin(lon_es)

    has_arc = (clat > 0.0) and (Re < Rgeo * clat - eps_vis)
    if has_arc:
        ratio = Re / (Rgeo * clat)
        ratio = max(-1.0, min(1.0, ratio))
        thmax = math.acos(ratio)
    else:
        thmax = 0.0

    px = Re * clat
    pz = Re * slat
    E0 = Rgeo * Rgeo + Re * Re
    F0 = -2.0 * Rgeo * px

    return ESCacheSph(lon_es, clon, slon, px, pz, Rgeo, E0, F0, thmax, has_arc)


# ---------------------------------------------------------------------
# Núcleo analítico (quartica + Newton) para obter |alpha| e lambda*
# ---------------------------------------------------------------------
def _dot_ratio(A: float, B: float, C: float, E0: float, F0: float, th: float) -> float:
    """cos(alpha) = (A + B cos(th) + C sin(th)) / sqrt(E0 + F0 cos(th))."""
    ct = math.cos(th)
    st = math.sin(th)
    num = A + B * ct + C * st
    den2 = E0 + F0 * ct
    if den2 <= 0.0:
        return -1.0
    return num / math.sqrt(den2)


def _newton_quartic_roots(q4, q3, q2, q1, q0, tol=1e-6, max_iter=50):
    """Newton-Raphson em f(x)=q4 x^4+...+q0, com sementes x=±1 (e x=0 como robustez)."""
    roots = []
    for x0 in (1.0, -1.0, 0.0):
        x = float(x0)
        for _ in range(max_iter):
            fx  = (((q4*x + q3)*x + q2)*x + q1)*x + q0
            fpx = ((4.0*q4*x + 3.0*q3)*x + 2.0*q2)*x + q1
            if fpx == 0.0:
                break
            x_next = x - fx / fpx
            if abs(x_next - x) < tol:
                x = x_next
                break
            x = x_next

        if -1.0 <= x <= 1.0 and math.isfinite(x):
            if not any(abs(x - r) < 1e-7 for r in roots):
                roots.append(x)
    return roots


def alpha_s1503_spherical_analytic_with_tiebreak(
    uN_x: float, uN_y: float, uN_z: float,
    cache: ESCacheSph,
    lon_ngso: float,
    *,
    uN_already_unit: bool = True,
    fallback_golden: bool = True,
    eps_dot: float = 1e-12,
) -> Tuple[float, float, float, float, float]:
    """
    Retorna:
      alpha_abs [rad], lambda_star [rad], theta_max [rad], dot_max, delta_long [rad]
    Aplicando o desempate por DeltaLong quando houver múltiplos mínimos.
    """
    if not cache.has_arc:
        return math.pi, cache.lon, 0.0, -1.0, 0.0

    # normalização (se necessário)
    if not uN_already_unit:
        n = math.sqrt(uN_x*uN_x + uN_y*uN_y + uN_z*uN_z)
        if n == 0.0:
            raise ValueError("uN must be non-zero.")
        uN_x, uN_y, uN_z = uN_x/n, uN_y/n, uN_z/n

    # rotaciona para frame lon_ES = 0
    uxp = uN_x * cache.clon + uN_y * cache.slon
    uyp = -uN_x * cache.slon + uN_y * cache.clon
    uzp = uN_z

    # f = A + B cos(th) + C sin(th)
    A = - (uxp * cache.px + uzp * cache.pz)
    B = cache.Rgeo * uxp
    C = cache.Rgeo * uyp

    E0, F0, thmax = cache.E0, cache.F0, cache.thmax

    # parâmetros {a,b,c,d,e} com G=0
    a = -2.0 * C * E0
    b = B * F0
    c = 2.0 * C * F0
    d = A * F0 - 2.0 * B * E0
    e = -C * F0

    # quartica em x=sin(th)
    q4 = e*e + b*b
    q3 = 2.0*d*e + 2.0*a*b
    q2 = d*d + 2.0*c*e + a*a - b*b
    q1 = 2.0*c*d - 2.0*a*b
    q0 = c*c - a*a

    # candidatos: bordas ±thmax + 0 + raízes Newton
    cand_th = [-thmax, thmax, 0.0]
    for x in _newton_quartic_roots(q4, q3, q2, q1, q0):
        th = math.asin(max(-1.0, min(1.0, x)))
        if -thmax - 1e-15 <= th <= thmax + 1e-15:
            cand_th.append(th)

    # avalia cos(alpha) e encontra máximo
    scored = []
    best_dot = -1e300
    for th in cand_th:
        dotv = _dot_ratio(A, B, C, E0, F0, th)
        if math.isfinite(dotv):
            scored.append((th, dotv))
            if dotv > best_dot:
                best_dot = dotv

    # fallback iterativo (raro; mas permitido)
    if fallback_golden and (not math.isfinite(best_dot)):
        phi = (1.0 + math.sqrt(5.0)) / 2.0
        inv_phi = 1.0 / phi
        a0, b0 = -thmax, thmax
        c0 = b0 - inv_phi * (b0 - a0)
        d0 = a0 + inv_phi * (b0 - a0)
        fc = _dot_ratio(A, B, C, E0, F0, c0)
        fd = _dot_ratio(A, B, C, E0, F0, d0)
        for _ in range(35):
            if (b0 - a0) < 1e-6:
                break
            if fc < fd:
                a0, c0, fc = c0, d0, fd
                d0 = a0 + inv_phi * (b0 - a0)
                fd = _dot_ratio(A, B, C, E0, F0, d0)
            else:
                b0, d0, fd = d0, c0, fc
                c0 = b0 - inv_phi * (b0 - a0)
                fc = _dot_ratio(A, B, C, E0, F0, c0)

        th_mid = 0.5 * (a0 + b0)
        cand2 = [(a0, _dot_ratio(A,B,C,E0,F0,a0)),
                 (b0, _dot_ratio(A,B,C,E0,F0,b0)),
                 (c0, fc),
                 (d0, fd),
                 (th_mid, _dot_ratio(A,B,C,E0,F0,th_mid))]
        th_best, best_dot = max(cand2, key=lambda t: t[1])
        scored = [(th_best, best_dot)]

    # desempate por DeltaLong: pegue candidatos com dot dentro de eps_dot do máximo
    tied = []
    for th, dotv in scored:
        if dotv >= best_dot - eps_dot:
            lam = _wrap_pi(cache.lon + th)      # LongAlpha
            dlong = _wrap_pi(lam - lon_ngso)    # DeltaLong = LongAlpha - LongNGSO
            tied.append((th, dotv, lam, dlong))

    if not tied:
        lam = cache.lon
        dlong = _wrap_pi(lam - lon_ngso)
        return math.pi, lam, thmax, -1.0, dlong

    # escolha: min |DeltaLong|; empate -> DeltaLong positivo
    def key_tiebreak(item):
        _, dotv, _, dlong = item
        return (abs(dlong), 0 if dlong > 0.0 else 1, -dotv)

    th_sel, dot_sel, lam_sel, dlong_sel = min(tied, key=key_tiebreak)

    dot_sel = max(-1.0, min(1.0, dot_sel))
    alpha_abs = math.acos(dot_sel)
    return alpha_abs, lam_sel, thmax, dot_sel, dlong_sel


# ---------------------------------------------------------------------
# Sinal de alpha (S.1503 §D6.4.4.3)
# ---------------------------------------------------------------------
def alpha_sign_s1503(uN_x: float, uN_y: float, uN_z: float, cache: ESCacheSph) -> int:
    """Retorna +1, 0, -1 conforme o critério geométrico do §D6.4.4.3."""
    # RES em ECEF (terra esférica)
    RES_x = cache.px * cache.clon
    RES_y = cache.px * cache.slon
    RES_z = cache.pz

    # Caso equador: alpha = -sign(uN_z)
    if abs(RES_z) < 1e-12:
        if abs(uN_z) < 1e-18:
            return 0
        return -1 if (uN_z > 0.0) else +1

    # λz=0 = -RES(z)/REN(z) (REN ~ uN)
    if abs(uN_z) < 1e-18:
        return -1  # quase paralelo ao plano XY

    lam_z0 = -RES_z / uN_z

    # λz=0 < 0 => interseção "atrás" da ES -> rho = infinity
    if lam_z0 < 0.0:
        rho = float("inf")
    else:
        Rx = RES_x + lam_z0 * uN_x
        Ry = RES_y + lam_z0 * uN_y
        rho = math.hypot(Rx, Ry)

    if RES_z > 0.0:  # norte
        if rho < cache.Rgeo - 1e-12:
            return +1
        if abs(rho - cache.Rgeo) <= 1e-12:
            return 0
        return -1
    else:            # sul
        if rho > cache.Rgeo + 1e-12:
            return +1
        if abs(rho - cache.Rgeo) <= 1e-12:
            return 0
        return -1


# ---------------------------------------------------------------------
# Função pública: recebe rE e rN, retorna alpha (e detalhes)
# ---------------------------------------------------------------------
def alpha_to_gso_arc(
    rE_km: np.ndarray,
    rN_km: np.ndarray,
    *,
    Re_km: float = RE_KM_DEFAULT,
    Rgeo_km: float = RGEO_KM_DEFAULT,
    signed: bool = True,
    output: AngleUnit = "rad",
    return_details: bool = True,
    eps_dot: float = 1e-12,
) -> Union[float, Dict[str, float]]:
    """
    Calcula o ângulo alpha mínimo ao arco GSO visível.

    rE_km, rN_km: vetores ECEF [km] da estação e do satélite NGSO.
    """
    rE = np.asarray(rE_km, float).reshape(3)
    rN = np.asarray(rN_km, float).reshape(3)

    nE = float(np.linalg.norm(rE))
    if nE == 0.0:
        raise ValueError("Earth station position vector cannot be zero.")
    if np.allclose(rE, rN):
        raise ValueError("Earth station and NGSO satellite cannot be coincident.")

    lon_es = _lon_from_ecef(rE)
    lat_es = math.asin(float(rE[2]) / nE)

    cache = es_cache_spherical(lat_es, lon_es, Re_km, Rgeo_km)

    uN = _to_unit(rN - rE)  # unit LOS ES->NGSO
    lon_ngso = _lon_from_ecef(rN)

    alpha_abs, lam_star, thmax, dotmax, dlong = alpha_s1503_spherical_analytic_with_tiebreak(
        float(uN[0]), float(uN[1]), float(uN[2]),
        cache, lon_ngso,
        uN_already_unit=True,
        fallback_golden=True,
        eps_dot=eps_dot,
    )

    if signed:
        sgn = alpha_sign_s1503(float(uN[0]), float(uN[1]), float(uN[2]), cache)
        alpha = sgn * alpha_abs
    else:
        sgn = +1
        alpha = alpha_abs

    # unit conversion
    if output == "deg":
        k = 180.0 / math.pi
        alpha_out = alpha * k
        alpha_abs_out = alpha_abs * k
        lam_star_out = lam_star * k
        thmax_out = thmax * k
        dlong_out = dlong * k
        lon_ngso_out = lon_ngso * k
    elif output == "rad":
        alpha_out = alpha
        alpha_abs_out = alpha_abs
        lam_star_out = lam_star
        thmax_out = thmax
        dlong_out = dlong
        lon_ngso_out = lon_ngso
    else:
        raise ValueError("output must be 'rad' or 'deg'")

    if not return_details:
        return float(alpha_out)

    return {
        "alpha": float(alpha_out),
        "alpha_abs": float(alpha_abs_out),
        "sign_alpha": float(sgn),
        "lambda_star": float(lam_star_out),
        "theta_max": float(thmax_out),
        "delta_long": float(dlong_out),
        "long_ngso": float(lon_ngso_out),
        "dot_max": float(dotmax),
        "has_arc": float(cache.has_arc),
    }
