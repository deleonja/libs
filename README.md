# QuantumWalks and QMB Mathematica Packages

This repository contains two Mathematica packages for quantum physics simulations:

1. **QuantumWalks.wl**: Simulates discrete-time quantum walks (DTQW) on 1D infinite lines.
2. **QMB.wl**: Provides tools for quantum many-body systems, quantum information, and quantum chaos.

---

## Table of Contents

- [Installation](#installation)
- [QuantumWalks Package](#quantumwalks-package)
- [QMB Package](#qmb-package)
  - [General Quantum Mechanics](#general-quantum-mechanics)
  - [Quantum Chaos & RMT](#quantum-chaos--rmt)
  - [Bose-Hubbard Model](#bose-hubbard-model)
  - [Spin Chains & Symmetries](#spin-chains--symmetries)
  - [Quantum Channels](#quantum-channels)
- [Requirements](#requirements)
- [Usage Examples](#usage-examples)
- [Contributing](#contributing)

---

## Installation

1. Place both `.wl` files and the `ForScience-0.88.45.paclet` in your working directory
2. Load packages in Mathematica:
   ```mathematica
   << "QuantumWalks`"
   << "QMB`"
   ```

## QuantumWalks Package

Functions for 1D discrete-time quantum walks:

| Function                               | Description                                               |
| -------------------------------------- | --------------------------------------------------------- |
| `Shift[t]`                             | Shift operator for DTQW at time `t`                       |
| `Coin[t]`                              | Hadamard coin operator at time `t` (accepts custom coins) |
| `DTQWStep[t]`                          | Unitary matrix for DTQW step                              |
| `DTQW[ψ0,t]`                           | Quantum state evolution after `t` steps                   |
| `PositionProbabilityDistribution[ψ,t]` | Position probabilities of state `ψ`                       |
| `ExpValPosition[ψ,t]`                  | Expected position value                                   |
| `L[θ,θa,θb]`                           | Losing strategy inequality (Parrondo's paradox)           |
| `W[θ,θa,θb]`                           | Winning strategy inequality (Parrondo's paradox)          |
| `CriticalAngle[avgPos]`                | Finds `θ` where average position is closest to zero       |

## QMB Package

### General Quantum Mechanics

| Function                           | Description                                                |
| ---------------------------------- | ---------------------------------------------------------- |
| `DensityMatrix[ψ]`                 | Density matrix of state vector `ψ`                         |
| `Pauli`                            | Pauli matrices (supports multi-qubit strings)              |
| `MatrixPartialTrace[mat,n,d]`      | Partial trace over subsystem `n`                           |
| `RandomQubitState[]`               | Haar-random qubit state                                    |
| `RandomChainProductState[L]`       | Haar-random product state of `L` qubits                    |
| `Commutator[A,B]`                  | Computes `AB - BA`                                         |
| `CommutationQ[A,B]`                | Checks if matrices `A` and `B` commute                     |
| `MutuallyCommutingSetQ[{A,B,...}]` | Checks if set of matrices `{A,B,..}` is mutually commuting |
| `Braket[ψ,ϕ]`                      | Inner product between state vectors `ψ` and `ϕ`            |
| `BlochVector[ρ]`                   | Bloch vector of single-qubit density matrix                |
| `Purity[ρ]`                        | Purity of density matrix `ρ`                               |
| `Concurrence[ρ]`                   | Two-qubit concurrence                                      |
| `Qubit[θ,φ]`                       | Qubit state with parameters `(θ,φ)` in the Bloch sphere    |
| `SU2Rotation`                      | SU(2) rotation matrix                                      |

---

### Quantum Chaos & RMT

| Function                             | Description                          |
| ------------------------------------ | ------------------------------------ |
| `MeanLevelSpacingRatio[eigenvalues]` | Average level spacing ratio 〈rₙ〉   |
| `IPR[ψ]`                             | Inverse participation ratio          |
| `kthOrderSpacings[spectrum,k]`       | k-th order level spacings            |
| `SpacingRatios[spectrum,k]`          | Level spacing ratios of order `k`    |
| `RatiosDistribution[r,β]`            | RMT level spacing ratio distribution |

---

### Bose-Hubbard Model

| Function                          | Description                                          |
| --------------------------------- | ---------------------------------------------------- |
| `BoseHubbardHamiltonian[N,L,J,U]` | Bose-Hubbard Hamiltonian for `N` bosons on `L` sites |
| `FockBasis[N,L]`                  | Fock basis states                                    |
| `BosonicPartialTrace[ρ]`          | Partial trace for bosonic systems                    |
| `RenyiEntropy[α,ρ]`               | α-th order Rényi entropy                             |

---

### Spin Chains & Symmetries

| Function                                   | Description                                             |
| ------------------------------------------ | ------------------------------------------------------- |
| `TranslationEigenvectorRepresentatives[L]` | Translation eigenvectors representatives for `L` qubits |
| `BlockDiagonalize[matrix]`                 | Block-diagonalizes matrices by symmetry                 |
| `IsingHamiltonian[hx,hz,J,L]`              | Ising model Hamiltonian                                 |
| `HeisenbergXXXwNoise[hz,L]`                | Noisy Heisenberg XXX chain Hamiltonian                  |

---

### Quantum Channels

| Function                           | Description                                |
| ---------------------------------- | ------------------------------------------ |
| `Reshuffle[m]`                     | Performs reshuffle operation over matrices |
| `FuzzyMeasurementChannel[ρ,p,...]` | Implements fuzzy measurement channels      |

---

## Requirements

- **Mathematica 12.0+**
- **ForScience paclet** (included, auto-installs on first load)

## Usage Examples

See demonstration notebooks in the repository for detailed usage examples of key functions.

## Contributing

Contributions are welcome! Please fork the repository and submit pull requests.
