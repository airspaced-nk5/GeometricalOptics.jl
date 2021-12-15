# Bundle, surface, and trace options

Each system for analysis takes some configuration of bundle and surfaces as input, and produces a trace.  The built-in options are summarized as follows (more details in [Reference](@ref)):

## Bundles
Physically, `x,y,z` are position and `angx,angy` are direction tangent angles (angles in radians). `zpos` is a fixed z position.  `Dx,Dy,Dz` are direction cosine (unit vector) components.
- [`bundle(x,y,angx,angy,zpos)`](@ref) is a bundle which takes either vectors in `x` and `y` or vectors in `angx` and `angy`.  All other arguments are scalar.  Vectors may be different length to form rectanglular array sampling. 
- [`bundle_as_array(x,y,angx,angy,zpos)`](@ref) takes matrices `Matrix{T} where T<:Real` for the first four arguments, and a scalar for `zpos::T where T<:Real`.
- [`bundle_as_array_big(x,y,z,Dx,Dy,Dz)`](@ref) takes all arguments as type `Matrix{T} where T<:Real`.


## Surfaces 
- [`zplane`](@ref) is a z-plane, with one coefficient used to describe the fixed z-location of the plane.
- [`spherical`](@ref) is a sphere, with one coefficient used to describe the z-location of the apex, and the other to describe the signed radius of curvature relative to that apex.
- [`aspherical`](@ref) is an aspherical surface, with coefficients describing the z-location, signed radius, conic constant, and arbitrary number of even order radial polynomial contributions.
- [`zernike`](@ref) is a zernike polynomial series representation. The normalization radius and arbitrary number of OSA/ANSI coefficients can be provided in `coeffs`. __Rays will terminate with `NaN` if they strike outside the normalization radius for this built-in function.__

## Traces and trace methods
Traces are generated when a type of bundle is sent through an [`opticalstack`](@ref) type.  Methods on traces either return organized information on the rays traced from the bundle, or directly return some basic diagnostic.
- [`trace`](@ref) is a Type which contains matrices of ray coordinates at all points in propagation.
- [`bigtrace`](@ref) is a Type which contains matrices of ray coordinates and direction cosines at all points in propagation.
- [`bigtrace_to_bundle`](@ref) returns a bundle object from the ray information at a particular position along propagation in a [`bigtrace`](@ref). 
Built-in methods to extract data from a trace:
- [`trace_extract_terminus`](@ref) returns requested ray information for all rays in [`trace`](@ref) or [`bigtrace`](@ref) as `Matrix` at the requested propagation position.
- [`trace_extract_ray`](@ref) returns requested ray information for the `i,j` ray in the [`trace`](@ref) or [`bigtrace`](@ref) at all positions in propagation as a vector.
Built-in evaluation methods on [`trace`](@ref) or [`bigtrace`](@ref) types return numbers or plots:
- [`rms_spot`](@ref) returns the rms-sense spot size of all rays in the trace at a particular position in propagation. 
- [`splot`](@ref) is a spot plot, and returns a plot of lateral ray coordinates of the bundle at a particular position in propagation. 
- [`rac`](@ref) returns a line plot of ray coordinates at one position in propagation as a function of the same coordinates at a different point in propagation.