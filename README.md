# Quantum Mechanics of Gravitational Potential Wells

This project explores the quantum behavior of a particle in a one-dimensional gravitational potential well, where the potential is defined as:

> **V(x) = α / x**,  
> with α = -G * m₁ * m₂ (Jupiter and Ganymede system)

We solve the **time-independent Schrödinger equation** analytically and simulate the bound states using MATLAB.

---

## 📌 Objectives
- Solve the Schrödinger equation for a gravitational potential
- Derive dimensionless wavefunctions using Laguerre polynomials
- Visualize normalized wavefunctions and probability densities
- Apply **SVD** and **optimization** to refine superpositions

---

## 🧮 Technologies Used
- MATLAB (Symbolic Math Toolbox, Optimization Toolbox)
- Linear Algebra (SVD, matrix formulations)
- Special Functions (Laguerre polynomials)
- Numerical integration & plotting

---

## 📂 File Structure

| File / Folder         | Description                                             |
|-----------------------|---------------------------------------------------------|
| `project.m`           | Main MATLAB script with symbolic + numeric simulations  |
| `report.pdf`          | Full technical report                                   |
| `ppt.pdf`             | Summary presentation                                    |

---

## 🔬 Mathematical Highlights
- Use of **dimensionless transformation** to handle scale mismatches
- Wavefunctions expressed as:

\[
\psi_n(x) = \text{(Normalization)} \cdot x \cdot e^{-x/(2na)} \cdot L_{n-1}^{(1)}(x/na)
\]

- Energy quantization:

\[
E_n = \frac{\mu \alpha^2}{2 n^2 \hbar^2}
\]
