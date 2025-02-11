# APM426-MAT1700-General-Relativity

# General Relativity Tensors in Python

This repository provides a collection of Python functions to compute fundamental tensors in General Relativity. Using symbolic computation with [Sympy](https://www.sympy.org/), the code calculates the Christoffel symbols, Riemann tensor, Ricci tensor, and Ricci scalar for a given metric tensor. The implementation is generalized to work in any \( n+1 \) dimensional space.

## Features

- **Generalized to \(n+1\) Dimensions:**  
  The code dynamically adapts to the number of coordinate variables provided.
- **Orientation Flexibility:**  
  Specify whether your metric is provided in the "down" (covariant) or "up" (contravariant) orientation.
- **Symbolic Computation:**  
  Utilizes Sympy to perform symbolic differentiation and simplification.
- **Modular Functions:**  
  Includes functions for:
  - Calculating all Christoffel symbols: `ChristoffelSymbols`
  - Retrieving a specific Christoffel symbol: `GiveChristoffel`
  - Computing the Riemann tensor: `RiemannTensor`
  - Retrieving a specific Riemann tensor component: `GiveRiemann`
  - Computing the Ricci tensor: `RicciTensor`
  - Retrieving a specific Ricci tensor component: `GiveRicci`
  - Computing the Ricci scalar: `RicciScalar`

## Requirements

- **Python 3.x**
- **Sympy:** Install via pip:
  ```bash
  pip install sympy
