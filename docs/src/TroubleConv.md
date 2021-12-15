# Troubleshooting and Conventions 

Geometrical optics usually relies on some choice of conventions.  This package is no exception. Some conventions describe choice of geometry and some describe potential pitfalls in analysis.  Some helpful tips are below.  Most of these tips are illustrated in the [Quickstart tutorial](@ref) and in [Examples (2D and 3D)](@ref).

- Physical coordinates are arbitrary; units need only be kept self-consistent.
- Rays may fail to properly trace with respect to surfaces when sufficiently steep.
- Index of refraction is defined __before__ its corresponding surface, where the light originates.
- Ray bundles and surfaces are all defined in global coordinates to minimize ambiguity.
- Reflection mode is defined by a sign change on the index of refraction after a given surface.
- Currently, rays must be known to intersect the surfaces in the prescribed sequence *a priori*. 
    - This is especially important to bear in mind for gradient index media where curved rays may significantly diverge.
- Surfaces, index functions, and phase functions defined by the user are not limited by form, but are limited by automatic-differentiability.
    - This will not likely be an issue with compositions of continuous math functions.
    - Failures are likely to arise if the chosen function cannot be understood by ForwardDiff.jl.
- 3D plotting does not always maintain proper aspect ratio, but 2D plotting does.  This is a consequence of dependency.
    - Where 1:1 aspect ratio is important, 2D plots should be used.
