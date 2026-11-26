# SEMFC: Structural Equation Modeling with Factors and Composites

##### Version: 1.0

##### Authors:
Alain LUONG, Arthur TENENHAUS

##### Key-words:
Structural Equation Modeling with Factors and Composites, (Restricted) Maximum Likelihood, Singular Value Decomposition

##### Contact:
arthur.tenenhaus@centralesupelec.fr

## Description
The SEMFC package implements Structural Equation Modeling (SEM) with factors and composites within the framework of the basic design (Tenenhaus et al., 2026). SEMFC supports both SVD-based and (restricted) maximum likelihood estimation methods. 

The SVD-based approach relies on a non-iterative algorithm for parameter estimation and yields consistent and asymptotically normal estimators, providing a statistically sound and computationally efficient solution.

The package also implements the Restricted Maximum Likelihood SEM (RML-SEM) approach for the basic design with factors and composites. SVD-based estimates are used as starting values for RML-SEM estimation. The RML-SEM estimator is obtained by solving a constrained nonlinear optimization problem (Kraft, 1988) using the optimization routines available in the nloptr package (Johnson, 2008), enabling efficient estimation under parameter constraints.

In addition, SEMFC provides comprehensive tools for statistical inference, including bootstrap procedures for SVD-based estimation and asymptotic inference for maximum likelihood estimation. The package also offers a wide range of goodness-of-fit measures for model evaluation, including the chi-square statistic, CFI, TLI, RMSEA, SRMR, AIC, and BIC.

## Installation
Required:

- Software: R (≥ 3.5.0)

- R libraries: see the [DESCRIPTION](https://github.com/Tenenhaus/SEMFC/blob/main/DESCRIPTION) file.

```
install.packages("SEMFC")
```

## Installation of a development branch from the git repository
Required:

- Software: R (≥ 3.5.0)

- R libraries: see the [DESCRIPTION](https://github.com/Tenenhaus/SEMFC/blob/main/DESCRIPTION) file.

- The R library `devtools`.

```
remove.packages("SEMFC")
devtools::install_github(repo="https://github.com/Tenenhaus/SEMFC.git", ref = "main")
```

## References
1. Tenenhaus, A., Tenenhaus, M., & Dijkstra, T. K. (2026). Structural equation modeling with factors and composites within the framework of the basic design. Advances in Data Analysis and Classification, 20(1), 227-254.
2. Johnson SG (2008). The NLopt nonlinear-optimization package. https://github.com/stevengj/nlopt.
3. Kraft D (1988). A software package for sequential quadratic programming. Technical Report DFVLR-FB 88-28, Institut für Dynamik der Flugsysteme, Oberpfaffenhofen.
