# QuantumWalks: a Wolfram Language package

The package offers optimized tools for simulations on 1D infinite lines and 2D confined geometries (like quantum billiards).

## 🚀 Key Features

* **1D DTQW:** Optimized implementation for the state evolution of a discrete-time quantum walk on an infinite line.
* **2D quantum billiards:** Modular framework for defining confined geometries.
    * Built-in support for the **Bunimovich Stadium**.
    * Automatic generation of shift operators ($W_m, W_n$) based on coordinate mapping.

## 📦 How to use the package?

1. Clone the `libs/` repository in your machine. 

2.  Load the package in your notebook:
    ```mathematica
    Get["route/to/libs/QuantumWalks/Kernel/init.m"]
    ```
    
    Make sure to change `route/to/libs/` by the route in your machine to `libs/`.

## Functions details and usage examples

See the notebooks in `UsageExamples/`.

## 🛠️ Dependencies

This package requires the **ForScience** library for usage message formatting.
* The package attempts to install it automatically if not found (`PacletInstall["ForScience"]`).
* Repository: [MMA-ForScience on GitHub](https://github.com/MMA-ForScience/ForScience).

## ⚠️ Technical Notes

* **2D Computational Basis:** The Hilbert space for billiards is organized as $|m, n, s\rangle$, where $(m,n)$ are coordinates and $s \in \{\uparrow, \downarrow\}$.
* **Boundaries:** The `BuildShiftOperators` function automatically handles boundary conditions (perfect reflection) when a node has no neighbor in the direction of movement.
* **Efficiency:** For large systems, avoid displaying the generated operators using `MatrixForm`, as they are large sparse arrays.

## 📜 License & Citation

This project is licensed under the **CC BY-NC-SA 4.0** (Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International).

### Use and Redistribution
* **You can** modify and redistribute this code for non-commercial purposes.
* **You must** give appropriate credit to the original author.
* **You cannot** use this package for commercial purposes (profit-making activities).

### 🎓 How to Cite
If you use **QuantumWalks** for research or published results, please cite this repository (or the associated paper) as follows:

> **de Leon, J. A. (2026).** *QuantumWalks: a WL package for Discrete-Time quantum walks*. GitHub Repository. https://github.com/deleonja/libs/tree/main/QuantumWalks

BibTeX entry:
```bibtex
@misc{QuantumWalksJAdeLeon,
  author = {de Leon, Jose Alfredo},
  title = {QuantumWalks: a WL package for Discrete-Time quantum walks},
  year = {2026},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{[https://github.com/deleonja/libs/tree/main/QuantumWalks](https://github.com/deleonja/libs/tree/main/QuantumWalks)}}
}
