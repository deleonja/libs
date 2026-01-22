# QuantumWalks

**QuantumWalks** es un paquete de Wolfram Language (Mathematica) diseñado para la simulación eficiente de Caminatas Cuánticas en Tiempo Discreto (DTQW).

El paquete ofrece herramientas optimizadas para simulaciones en 1D (líneas infinitas) y geometrías confinadas en 2D (billares cuánticos), aprovechando `SparseArray` y operaciones vectorizadas para maximizar el rendimiento computacional.

## 🚀 Características Principales

* **DTQW en 1D:** Implementación clásica con moneda de Hadamard y monedas arbitrarias. Incluye análisis de paradojas de Parrondo.
* **Billares Cuánticos 2D:** Framework modular para definir geometrías confinadas.
    * Soporte para **Estadio de Bunimovich**.
    * Generación automática de operadores de desplazamiento ($W_m, W_n$) basada en mapas de coordenadas.
* **Alto Rendimiento:** Uso extensivo de álgebra lineal dispersa (`SparseArray`) y `KroneckerProduct` para minimizar el uso de memoria en grandes espacios de Hilbert.
* **Análisis:** Herramientas integradas para calcular distribuciones de probabilidad y valores esperados.

## 📦 Instalación

1. Descarga la carpeta `QuantumWalks`.
2. Mueve la carpeta al directorio de aplicaciones de usuario de Mathematica. Puedes encontrar esta ruta ejecutando el siguiente comando en un notebook:
   ```wolfram
   FileNameJoin[{$UserBaseDirectory, "Applications"}]
