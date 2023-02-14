# basis
Implementation of basis sets for galactic dynamics.

Contents of this repository:
- `Isochrone.jl`
  - Implemention in Julia of the isochrone-adapted basis set using the
    'discretized Stieltjes procedure'. The comments describe where to
    make change so as to adapt to a different zeroth-order.
  - Helper functions are provided that hook into JBF's mode-recovery code
  - A function is provided (`checkOrthogIso`) that checks
    (bi)orthogonality for the constructed basis set (integrating up an
    N*N table of potential/density functions)
- `isochrone_expressions.m`, 
  - Mathematica code that generates symbolic expressions for the
    isochrone-adapted basis set using the 'modified Chebyshev
    algorithm'. These agree with the Julia code implementation (up to
    numerical error which is inherent in the Chebyshev approach).
- `spherical_expressions.m`, `disc_expressions.m`
  - Mathematica code demonstrating symbolic expressions for the
    'known' basis sets mentioned in the paper, implemented explicitly
    using hypergeometric polynomials in a certain differential
    operator.






			       
