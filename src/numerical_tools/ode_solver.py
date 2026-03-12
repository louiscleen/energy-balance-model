"""
ode_solver.py

Petit module pour résoudre des EDO du premier ordre avec la méthode d'Euler explicite.
"""

from __future__ import annotations
from typing import Callable, Sequence
import numpy as np


ArrayLike = Sequence[float] | np.ndarray
RHSFunction = Callable[[float, np.ndarray], np.ndarray]


def euler_explicit(
    f: RHSFunction,
    t_span: tuple[float, float],
    y0: ArrayLike,
    h: float,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Résout y'(t) = f(t, y) sur [t0, tf] avec la méthode d'Euler explicite.

    Paramètres
    ----------
    f : callable
        Fonction du membre de droite. Signature : f(t, y).
    t_span : tuple(float, float)
        Intervalle (t0, tf).
    y0 : array-like
        Condition initiale (scalaire ou vecteur).
    h : float
        Pas de temps.

    Retour
    ------
    t : ndarray
        Vecteur des temps.
    y : ndarray
        Solution approchée.
    """

    t0, tf = t_span

    if h <= 0:
        raise ValueError("Le pas h doit être > 0")
    if tf < t0:
        raise ValueError("t_span doit être (t0, tf) avec tf >= t0")

    y0 = np.atleast_1d(np.asarray(y0, dtype=float))
    is_scalar = y0.size == 1

    n_steps = int(np.ceil((tf - t0) / h))

    t = np.zeros(n_steps + 1)
    y = np.zeros((n_steps + 1, y0.size))

    t[0] = t0
    y[0] = y0

    for n in range(n_steps):

        dt = min(h, tf - t[n])

        f_val = np.asarray(f(t[n], y[n]), dtype=float)

        if f_val.shape != y[n].shape:
            raise ValueError("f(t, y) doit retourner un tableau de même taille que y")

        y[n + 1] = y[n] + dt * f_val
        t[n + 1] = t[n] + dt

    if is_scalar:
        return t, y[:, 0]

    return t, y


if __name__ == "__main__":

    # Exemple : y' = -y
    def f(t, y):
        return -y

    t, y = euler_explicit(
        f=f,
        t_span=(0.0, 2.0),
        y0=1.0,
        h=0.1
    )

    print("t =", t[:5], "...")
    print("y =", y[:5], "...")