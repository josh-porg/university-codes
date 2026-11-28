"""Symbolic thermodynamic property relations from an equation of state.

Ported from AE 573 (``GasLaw`` and ``Relation``). Requires sympy:
``pip install "unicodes[symbolic]"``.

    gl = GasLaw("p*v = R*T", cp="cp0")
    gl.partials["PrTcV"]        # (dp/dT) at constant v  ->  R/v
    gl.cv                       # cp - T (dv/dT)_p^2 (dp/dv)_T  ->  cp0 - R

Partial-derivative keys read ``<A>r<B>c<C>`` = (dA/dB) holding C constant,
as in the MATLAB struct ``formulas.Partial``.
"""

from __future__ import annotations

from collections.abc import Iterable


def _sympy():
    try:
        import sympy
    except ImportError as exc:
        raise ImportError('Symbolic relations need sympy: pip install "unicodes[symbolic]"') from exc
    return sympy


_FUNCTIONS = {"sqrt", "exp", "log", "sin", "cos", "tan", "asin", "acos", "atan", "sinh", "cosh", "tanh", "pi", "E"}


def _parse_equation(text: str, local=None):
    """Expression equal to zero for ``lhs = rhs``. MATLAB ``^`` is accepted, and every identifier
    that is not a standard function is a plain symbol (so ``gamma`` or ``beta`` are not sympy functions)."""
    import re

    sp = _sympy()
    text = text.replace("^", "**")
    names = set(re.findall(r"[A-Za-z_]\w*", text)) - _FUNCTIONS
    local = {**{n: sp.Symbol(n) for n in names}, **(local or {})}
    lhs, _, rhs = text.replace("==", "=").partition("=")
    expr = sp.sympify(lhs, locals=local)
    if rhs:
        expr = expr - sp.sympify(rhs, locals=local)
    return expr  # expression equal to zero


class GasLaw:
    """Partial derivatives, Maxwell relations and property changes for an EoS.

    Parameters
    ----------
    equation_of_state : relation between ``p``, ``v`` and ``T``, e.g.
        ``"p*v = R*T"`` or ``"(p + a/v**2)*(v - b) = R*T"``.
    cp : expression for the constant-pressure specific heat.
    """

    def __init__(self, equation_of_state: str, cp: str):
        sp = _sympy()
        self.p, self.v, self.T = p, v, T = sp.symbols("p v T")
        local = {"p": p, "v": v, "T": T}
        F = _parse_equation(equation_of_state, local)
        self.equation_of_state = sp.Eq(F, 0)

        def d(a, b):
            # implicit derivative (da/db) at constant third variable, from F(p, v, T) = 0
            return sp.simplify(-sp.diff(F, b) / sp.diff(F, a))

        P = {
            "PrTcV": d(p, T),
            "PrVcT": d(p, v),
            "TrPcV": d(T, p),
            "TrVcP": d(T, v),
            "VrTcP": d(v, T),
            "VrPcT": d(v, p),
        }
        self.cp = sp.sympify(cp, locals=local)

        # Maxwell relations and energy/enthalpy derivatives
        P["SrVcT"] = P["PrTcV"]
        P["SrPcT"] = -P["VrTcP"]
        P["UrVcT"] = sp.simplify(T * P["PrTcV"] - p)  # internal pressure
        P["UrPcT"] = sp.simplify(P["UrVcT"] * P["VrPcT"])
        P["HrPcT"] = sp.simplify(v - T * P["VrTcP"])
        P["HrVcT"] = sp.simplify(T * P["PrTcV"] + v * P["PrVcT"])

        self.cv = sp.simplify(self.cp + T * P["VrTcP"] ** 2 * P["PrVcT"])
        P["HrTcP"] = self.cp
        P["UrTcV"] = self.cv
        P["SrTcP"] = self.cp / T
        P["SrTcV"] = self.cv / T
        P["TrVcU"] = sp.simplify((p - T * P["PrTcV"]) / self.cv)  # Joule coefficient
        P["TrPcH"] = sp.simplify((T * P["VrTcP"] - v) / self.cp)  # Joule-Thomson coefficient
        P["PrVcS"] = sp.simplify(self.cp / self.cv * P["PrVcT"])
        self.partials = P

        # Reciprocity checks (the MATLAB printed a message when these failed)
        self.reciprocity_ok = all(
            sp.simplify(P[a] * P[b] - 1) == 0
            for a, b in (("TrVcP", "VrTcP"), ("PrVcT", "VrPcT"), ("TrPcV", "PrTcV"))
        )

        # Property changes as (indefinite) integrals along each variable,
        # as in the MATLAB ``formulas.Delta``
        self.delta = {
            "SrTV": sp.integrate(P["SrTcV"], T) + sp.integrate(P["SrVcT"], v),
            "SrTP": sp.integrate(P["SrTcP"], T) + sp.integrate(P["SrPcT"], p),
            "HrTP": sp.integrate(P["HrTcP"], T) + sp.integrate(P["HrPcT"], p),
            "UrTV": sp.integrate(P["UrTcV"], T) + sp.integrate(P["UrVcT"], v),
        }


    def in_terms_of(self, expr, variable: str = "p"):
        """Eliminate ``variable`` from ``expr`` using the equation of state.

        Implicit derivatives contain all of ``p, v, T``; e.g. the ideal-gas
        internal pressure ``R*T/v - p`` only becomes 0 once ``p = R*T/v``
        is substituted.
        """
        sp = _sympy()
        sym = {"p": self.p, "v": self.v, "T": self.T}[variable]
        solutions = sp.solve(self.equation_of_state, sym)
        if len(solutions) != 1:
            raise ValueError(f"Equation of state is not uniquely solvable for {variable}")
        return sp.simplify(expr.subs(sym, solutions[0]))


class Relation:
    """An equation that can be solved for whichever single variable is unknown.

        r = Relation("p*v = R*T", positive=["p", "v", "T", "R"])
        r.solve(p=101325, R=287, T=288)   # -> [v]

    Parameters
    ----------
    equation : the relation, written with ``=``.
    name, assumptions, tags : free-form metadata (equation name, conditions
        under which it holds, search tags).
    positive : variables assumed positive (and real) when solving.
    """

    def __init__(
        self,
        equation: str,
        name: str | None = None,
        assumptions: Iterable[str] = (),
        tags: Iterable[str] = (),
        positive: Iterable[str] = (),
    ):
        sp = _sympy()
        self.name = name or equation
        self.equation_text = equation
        self.assumptions = list(assumptions)
        self.tags = list(tags)
        plain = _parse_equation(equation)
        positive = set(positive)
        self.symbols = {
            s.name: sp.Symbol(s.name, positive=True) if s.name in positive else sp.Symbol(s.name, real=True)
            for s in plain.free_symbols
        }
        self.expression = plain.subs({sp.Symbol(n): s for n, s in self.symbols.items()})

    @property
    def variables(self) -> list[str]:
        return sorted(self.symbols)

    def solve(self, **known) -> list:
        """Substitute the known values and solve for the one remaining variable."""
        sp = _sympy()
        unknown = [n for n in self.symbols if n not in known]
        if len(unknown) != 1:
            raise ValueError(f"Exactly one unknown required, got {unknown}")
        expr = self.expression.subs({self.symbols[n]: val for n, val in known.items()})
        return sp.solve(expr, self.symbols[unknown[0]])

    def solve_numeric(self, **known) -> list[float]:
        """Real roots for the one unknown, found numerically.

        Polynomials of degree <= 2 are solved exactly; anything else is
        scanned on a logarithmic grid (both signs unless the variable is
        assumed positive) and every sign change refined with ``brentq``.
        Much faster than :meth:`solve` for relations with non-integer powers.
        """
        import numpy as np
        from scipy.optimize import brentq

        sp = _sympy()
        unknown = [n for n in self.symbols if n not in known]
        if len(unknown) != 1:
            raise ValueError(f"Exactly one unknown required, got {unknown}")
        x = self.symbols[unknown[0]]
        expr = self.expression.subs({self.symbols[n]: val for n, val in known.items()})
        try:
            poly = sp.Poly(expr, x)
            if poly.degree() <= 2:
                out = []
                for r in np.roots([complex(c) for c in poly.all_coeffs()]):
                    if abs(r.imag) < 1e-9 * max(1, abs(r)):
                        out.append(float(r.real))
                return out
        except (sp.PolynomialError, sp.GeneratorsNeeded, TypeError, ValueError):
            pass
        f = sp.lambdify(x, expr, "numpy")

        def fr(z):
            with np.errstate(all="ignore"):
                v = complex(f(z))
            return v.real if abs(v.imag) < 1e-12 * max(1.0, abs(v.real)) else np.nan

        mags = np.logspace(-10, 10, 2001)
        grid = mags if x.is_positive else np.concatenate([-mags[::-1], [0.0], mags])
        vals = np.array([fr(g) for g in grid])
        roots = []
        ok = np.isfinite(vals[:-1]) & np.isfinite(vals[1:]) & (np.sign(vals[:-1]) != np.sign(vals[1:]))
        for i in np.flatnonzero(ok):
            if vals[i] == 0:
                roots.append(float(grid[i]))
                continue
            r = brentq(fr, grid[i], grid[i + 1], xtol=1e-14, rtol=1e-13)
            if abs(fr(r)) < 1e-6 * max(1.0, abs(vals[i]), abs(vals[i + 1])):  # reject poles
                roots.append(float(r))
        return roots


POSITIVE_PREFIXES = ("M_", "p_", "pt_", "T_", "Tt_", "a_", "rho_")


def expand_stations(templates: Iterable[str], stations: Iterable, token: str = "_st") -> list[str]:
    """Repeat each template equation for every station, replacing ``token`` (``equationRepeater``)."""
    return [t.replace(token, f"_{s}") for t in templates for s in stations]


def _is_positive(name, prefixes):
    return name == "alpha" or any(name.startswith(p) for p in prefixes)


def build_relations(equations: Iterable[str], positive_prefixes=POSITIVE_PREFIXES) -> list[Relation]:
    """Relations with Mach numbers, pressures, temperatures, speeds of sound and densities assumed positive,
    as the AE 573 scripts set up through their assumption dictionaries."""
    out = []
    for eq in equations:
        names = {s.name for s in _parse_equation(eq).free_symbols}
        out.append(Relation(eq, positive=[n for n in names if _is_positive(n, positive_prefixes)]))
    return out


def solve_relations(relations: Iterable[Relation], knowns: dict, verbose: bool = False) -> dict:
    """Repeatedly solve every relation that has exactly one unknown until nothing changes (``solveRelations``).

    ``knowns`` maps names to values (``None`` or NaN = unknown). Real roots
    are kept; with several, the first positive one is taken. Relations that
    fail to solve are skipped, as the MATLAB caught and reported the error.
    """
    import math

    values = {k: v for k, v in knowns.items() if v is not None and not (isinstance(v, float) and math.isnan(v))}
    relations = list(relations)
    changed = True
    while changed:
        changed = False
        for rel in relations:
            unknown = [n for n in rel.variables if n not in values]
            if len(unknown) != 1:
                continue
            try:
                real = rel.solve_numeric(**{n: values[n] for n in rel.variables if n in values})
            except Exception as exc:  # noqa: BLE001
                if verbose:
                    print(f"could not solve {rel.equation_text} for {unknown[0]}: {exc}")
                continue
            if not real:
                continue
            pos = [r for r in real if r > 0]
            values[unknown[0]] = (pos or real)[0]
            changed = True
            if verbose:
                print(f"solved {unknown[0]} = {values[unknown[0]]:g}")
    return values
