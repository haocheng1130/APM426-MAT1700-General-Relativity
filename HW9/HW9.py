"""
radial_geodesic.py
------------------
This script computes the value of E and the corresponding r_max for a radial
geodesic in Schwarzschild spacetime under the condition that the round-trip
coordinate time is T = 160*pi*M. It then calculates the proper time for Peter.
All calculations are done in geometric units (G = c = 1) with M = 1.
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

# Define the Schwarzschild metric function f(r) = 1 - 2M/r
def f(r, M=1.0):
    return 1.0 - 2.0*M/r

# Compute the round-trip coordinate time T(E)
def T_of_E(E, M=1.0):
    """
    Compute T(E) = 2 * ∫[r=4M to r_max] { E / ( f(r)*sqrt(E^2 - f(r)) ) } dr,
    where the turning point r_max = 2M/(1-E^2) (from E^2 = f(r_max)).
    """
    if E >= 1.0:
        return np.inf

    r_i = 4.0 * M
    r_max = 2.0 * M / (1.0 - E**2)
    if r_max <= r_i:
        return np.inf

    def integrand(r):
        return E / ( f(r, M) * np.sqrt(E**2 - f(r, M)) )

    val, err = quad(integrand, r_i, r_max, limit=100)
    return 2.0 * val

# Compute the proper time (tau) for Peter along the radial geodesic
def tau_of_E(E, M=1.0):
    """
    Compute tau(E) = 2 * ∫[r=4M to r_max] { 1 / sqrt(E^2 - f(r)) } dr,
    where r_max = 2M/(1-E^2) is determined by the turning point condition.
    """
    if E >= 1.0:
        return np.inf

    r_i = 4.0 * M
    r_max = 2.0 * M / (1.0 - E**2)
    if r_max <= r_i:
        return np.inf

    def integrand(r):
        return 1.0 / np.sqrt(E**2 - f(r, M))
    
    val, err = quad(integrand, r_i, r_max, limit=100)
    return 2.0 * val

# Auto-search for a bracket [E_left, E_right] where T_of_E(E) - T_target changes sign.
def find_bracket_for_root(func, start=0.9, step=1e-4, max_e=0.9999999999):
    """
    Starting from 'start', increment E by 'step' until func(E) changes sign.
    Returns a tuple (E_left, E_right) if found, otherwise None.
    """
    E_left = start
    f_left = func(E_left)
    
    E = E_left
    while E < max_e:
        E += step
        f_right = func(E)
        if f_left * f_right < 0:
            return (E - step, E)
        f_left = f_right
    return None

# Wrapper function for root-finding: f(E) = T_of_E(E) - T_target
def func_to_solve(E, T_target, M=1.0):
    return T_of_E(E, M) - T_target

def find_E_for_T(T_target, M=1.0, initial_bracket=(0.9, 0.999999)):
    """
    Finds E in (0,1) such that T_of_E(E) = T_target using Brentq.
    It first searches for a bracket where the function changes sign.
    """
    E_left, E_right = initial_bracket
    f_left = func_to_solve(E_left, T_target, M)
    f_right = func_to_solve(E_right, T_target, M)
    
    if f_left * f_right >= 0:
        print("Initial bracket does not bracket a root. Auto-searching for bracket...")
        bracket = find_bracket_for_root(lambda E: func_to_solve(E, T_target, M),
                                        start=E_left, step=1e-4, max_e=0.9999999999)
        if bracket is None:
            raise ValueError("Unable to find a bracket where the function changes sign.")
        E_left, E_right = bracket
        print(f"Found bracket: E_left = {E_left}, E_right = {E_right}")
    else:
        print(f"Using initial bracket: E_left = {E_left}, E_right = {E_right}")
    
    E_solution = brentq(lambda E: func_to_solve(E, T_target, M), E_left, E_right)
    return E_solution

if __name__ == "__main__":
    # Set parameters: M = 1 (geometric units) and target coordinate time T_target = 160*pi*M
    M_val = 1.0
    T_target = 160.0 * np.pi * M_val

    # Find E such that T_of_E(E) = T_target.
    try:
        E_sol = find_E_for_T(T_target, M_val, initial_bracket=(0.9, 0.999999))
    except ValueError as e:
        print("Error:", e)
        exit(1)
    
    # Compute r_max from the turning point condition: r_max = 2M/(1-E^2)
    r_max_sol = 2.0 * M_val / (1.0 - E_sol**2)

    # Compute the proper time for Peter
    tau_sol = tau_of_E(E_sol, M_val)

    # Print the results
    print("=== Computation Results ===")
    print(f"Target coordinate time T_target = {T_target:.6f}")
    print(f"Found E = {E_sol:.8f}")
    print(f"Corresponding r_max = {r_max_sol:.6f}")
    print(f"Coordinate time T(E_sol) = {T_of_E(E_sol, M_val):.6f}")
    print(f"Proper time for Peter, tau = {tau_sol:.6f}")

    # For comparison, note that Paul's proper time (for 10 orbits at r = 4M) is 80*pi*M
    paul_tau = 80.0 * np.pi * M_val
    print(f"Paul's proper time = {paul_tau:.6f}")
    print(f"Proper time difference (Peter - Paul) = {tau_sol - paul_tau:.6f}")
