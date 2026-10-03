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


def _parse_equation(text: str, local=None):
    sp = _sympy()
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
