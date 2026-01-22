# QuantumWalks

**QuantumWalks** is a Wolfram Language (Mathematica) package designed for the efficient simulation of Discrete Time Quantum Walks (DTQW).

The package offers optimized tools for simulations on 1D infinite lines and 2D confined geometries (quantum billiards), leveraging `SparseArray` and vectorized operations to maximize computational performance and memory efficiency.

## 🚀 Key Features

* **1D DTQW:** Implementation of walks on infinite lines with standard Hadamard coins and arbitrary SU(2) coins. Includes analysis tools for Parrondo's paradox (`L`, `W` inequalities).
* **2D Quantum Billiards:** Modular framework for defining confined geometries.
    * Built-in support for the **Bunimovich Stadium**.
    * Automatic generation of shift operators ($W_m, W_n$) based on coordinate mapping.
* **High Performance:** Extensive use of sparse linear algebra (`SparseArray`) and `KroneckerProduct` to handle large Hilbert spaces efficiently.
* **Analysis:** Integrated tools to calculate position probability distributions and expected values.

## 📦 Installation

1.  Download the `QuantumWalks` folder.
2.  Move the folder to your Mathematica user applications directory. You can find this path by running the following command in a notebook:
    ```wolfram
    FileNameJoin[{$UserBaseDirectory, "Applications"}]
    ```
3.  Restart Mathematica or the Kernel so the package is recognized.

## 🛠️ Dependencies

This package requires the **ForScience** library for usage message formatting.
* The package attempts to install it automatically if not found (`PacletInstall["ForScience"]`).
* Repository: [MMA-ForScience on GitHub](https://github.com/MMA-ForScience/ForScience).

## 💻 Quick Start Guide

To load the full package:

```wolfram
<< QuantumWalks`
```

## 📂 Package Structure

```text
QuantumWalks/
├── Kernel/
│   └── init.m              # Initialization and dependency loader
├── Billiards/
│   ├── Common.wl           # Base 2D logic (BuildGridShiftOperators)
│   └── Bunimovich.wl       # Specific Bunimovich Stadium geometry
└── DQWL.wl                 # 1D Logic (Shift, Coin, DTQW, Parrondo)

## ⚠️ Technical Notes

* [cite_start]**2D Computational Basis:** The Hilbert space for billiards is organized as $|m, n, s\rangle$, where $(m,n)$ are coordinates and $s \in \{\uparrow, \downarrow\}$[cite: 245, 246, 247].
* [cite_start]**Boundaries:** The `BuildGridShiftOperators` function automatically handles boundary conditions (perfect reflection) when a node has no neighbor in the direction of movement[cite: 248, 250].
* [cite_start]**Efficiency:** For large systems, avoid displaying the generated operators using `MatrixForm`, as they are large sparse arrays[cite: 271].

## 📄 License

[Insert your preferred license here, e.g., MIT, GPL]
