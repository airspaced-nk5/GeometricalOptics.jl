module GeometricalOptics

export zplane, spherical, aspherical, zernike, bundle, bundle_fast, bundle_as_array, bundle_as_array_fast, opticalstack, 
rms_spot, splot, rac, trace, bigtrace, bigtrace_to_bundle, bundle_as_array_big, bundle_as_array_big_fast,
trace_extract_ray, trace_extract_terminus

using Plots
using PlotlyBase
using PlotlyKaleido
# Plots.plotly()
# Plots.PlotlyBackend()
using ForwardDiff
using LinearAlgebra
using ZernikePolynomials

"""
    zplane(x, y, coeffs)

Composite type representing a plane z(x, y) = zpos. Argument coeffs of type `Vector{T} where T<:Real` has
entries `coeffs = [zpos]`. 

# Examples
```julia-repl
julia> zplane(0., 2., [3.])
3.0
```
"""
function zplane(x::T where T<:Real,y::T where T<:Real,coeffs::Vector{T} where T<:Real)
    return coeffs[1]
end


"""
    zernike(x, y, coeffs)

Composite type representing OSA/ANSI type Zernike polynomial surface as a function of lateral coordinate `x, y`. Argument `coeffs` of type `Vector{T} where T<:Real` has
entries `coeffs = [rad, C0, C1, C2, ...]`. The normalization radius is `rad` and the corresponding Zernike surface is given
by the contributions of the Zⱼ with arguments over the cartesian coordinates normalized to the unit disk, i.e.:

``z(x, y) = ∑ⱼCⱼZⱼ(r, θ)`` where ``r ≡ √( x² + y² ) / rad`` and ``θ ≡ atan(y / x)``

__NOTE this function evaluates to `NaN` outside of the normalization radius.__

# Examples
```julia-repl
julia> zernike(0., 2., [2., 0., 1.])
2.0
```
"""
function zernike(x::T where T<:Real,y::T where T<:Real,coeffs::Vector{T} where T<:Real)
    normRad=coeffs[1]
    if (x/normRad)^2 + (y/normRad)^2 >1.0
        return NaN
    else
        z=0.0
        for i=1:(length(coeffs)-1)
            Z=Zernike(i-1,index=:OSA,coord=:cartesian)
            z=z+coeffs[i+1]*Z(x/normRad,y/normRad)
        end
        return z
    end
end


"""
    spherical(x, y, coeffs)

Composite type representing spherical surface with apex at defined point. Argument `coeffs` of type `Vector{T} where T<:Real` has
entries `coeffs = [zpos, signed_radius]`. Here `zpos` is the apex and `signed_radius` has magnitude equal to the radius of the sphere 
and sign which sets the direction of curvature of the sphere.

# Examples
```julia-repl
julia> spherical(0., 2., [1., -10.])
0.7979589711327115
```
"""
function spherical(x::T where T<:Real,y::T where T<:Real,coeffs::Vector{T} where T<:Real)
    z=coeffs[1]+coeffs[2]-sign(coeffs[2])*(sqrt(coeffs[2]^2-x^2-y^2))
    return z
end


"""
    aspherical(x, y, coeffs)

Composite type representing aspherical surface as a function of lateral coordinate `x,y`. Argument `coeffs` of type `Vector{T} where T<:Real` has
entries `coeffs = [zpos, R, ``κ``, A4, A6, ...]` where the surface is 

``z(x,y) = zpos + r² / R / (1 + √(1 - (1 + κ)r² / R²)) + A₄r⁴ + A₆r⁶ + ...`` where ``r² ≡ x² + y²``

# Examples
```julia-repl
julia> aspherical(0., 2., [1., -10., -4., 0.001])
1.0616649185805458
```
"""
function aspherical(x::T where T<:Real,y::T where T<:Real,coeffs::Vector{T} where T<:Real)
    rsq=x^2+y^2
    Rsq=coeffs[2]^2
    zb=coeffs[1]+rsq/coeffs[2]/(1+sqrt(1-(1+coeffs[3])*rsq/coeffs[2]^2))
    for j=1:length(coeffs)-3 
        zb=zb+coeffs[j+3]*rsq^(2*(1+j))
    end 
    return zb
end


struct impz{coeffs,surfDef}
    coeffs::coeffs 
    surfDef::surfDef
end 
function (iz::impz)(r)
    zret=iz.surfDef(r[1],r[2],iz.coeffs)
    return (r[3]-zret)^2
end

struct normalVec{impSurfIn}
    impSurfIn::impSurfIn
end
function (nv::normalVec)(R)
    gs=R -> ForwardDiff.gradient(nv.impSurfIn,R)
    gs(R) ./ sqrt(dot(gs(R),gs(R)))
end

struct refrSf{impSurf}
    impSurf::impSurf 
end
function (rf::refrSf)(R0,D0,n0,n1;isReorient=true)
    nvf=normalVec(rf.impSurf)
    N=nvf(R0)
    if n1!=-n0
        if isReorient
            sg=sign(dot(N,D0))
            N=sg*N
        end
        if (n1^2-n0^2+n0^2*dot(N,D0)^2)<0; @error "TIR violates requested refraction for current ray"; end;
        k = sqrt(n1^2-n0^2+n0^2*dot(N,D0)^2)-n0*dot(N,D0)
    else
        k=sqrt(n1^2-n0^2+n0^2*dot(N,D0)^2)-n0*dot(N,D0)
    end
    (k*N+n0*D0)/n1
end

struct getIntersect_preInit{impSurfIn,sInit}
    impSurfIn::impSurfIn
    sInit::sInit
end
function (gi::getIntersect_preInit)(R,D)
    # for start position R and direction cosine vector D, implicit surface impSurfIn, and coefficient coeff,
    f(s)=gi.impSurfIn(R.+D.*s)
    gf=s->ForwardDiff.derivative(f,s)
    sg=gi.sInit
    tol=0.000001
    adj=tol+1.0
    while abs(adj)>tol
        adj=f(sg)/gf(sg)
        sg=sg-adj
    end
    sg
end

VecOrNum=Union{AbstractVector{T},T} where T

"""
    bundle(x, y, angx, angy, zpos)

Composite type representing the origination of a bundle of rays. Input zpos is type `T<:Real` 
and is z-plane start position for the rays.
Other four arguments are typed in one of the following two patterns:
1. x and y are Type `AbstractArray{T,1} where T<:Real`, angx and angy are Type `T where T<:Real`
2. x and y are Type `T where T<:Real`, angx and angy are Type `AbstractArray{T,1} where T<:Real`}

When a variable is declared this type with specific parameters and called 
as a function with no arguments, it returns a tuple (position, direction) each 
entry of which is an output array of vectors corresponding to the input information.

# Examples
```julia-repl
julia> eval_bundle = bundle([0., 1.], [0.,1.], 0., 0., 0.)
bundle([0.0, 1.0], [0.0, 1.0], 0.0, 0.0, 0.0)

julia> eval_bundle()
([[0.0, 0.0, 0.0] [0.0, 1.0, 0.0]; [1.0, 0.0, 0.0] [1.0, 1.0, 0.0]], [[0.0, 0.0, 1.0] [0.0, 0.0, 1.0]; [0.0, 0.0, 1.0] [0.0, 0.0, 1.0]])
```
"""
struct bundle{x<:VecOrNum,
                y<:VecOrNum,
                angx<:VecOrNum,
                angy<:VecOrNum,
                zpos<:Real}
    x::x
    y::y
    angx::angx# these are direction tangents.
    angy::angy
    zpos::zpos
end
struct bundle_fast{T}
    x::VecOrNum{T}
    y::VecOrNum{T}
    angx::VecOrNum{T}# these are direction tangents.
    angy::VecOrNum{T}
    zpos::T
end
function (cb::Union{bundle,bundle_fast})()
    if typeof(cb.x)<:AbstractArray && typeof(cb.y)<:AbstractArray && typeof(cb.angx)<:Real && typeof(cb.angy)<:Real
        R0=[[cb.x[i].+eps(),cb.y[j].+eps(),cb.zpos.+eps()]  for i in 1:length(cb.x),j in 1:length(cb.y)]
        D0=[[sin(cb.angx),sin(cb.angy),sqrt(1-sin(cb.angx)^2-sin(cb.angy)^2)] for i in 1:length(cb.x),j in 1:length(cb.y)]
    elseif typeof(cb.angx)<:AbstractArray && typeof(cb.angy)<:AbstractArray && typeof(cb.x)<:Real && typeof(cb.y)<:Real
        R0=[[cb.x.+eps(),cb.y.+eps(),cb.zpos.+eps()]  for i in 1:length(cb.angx),j in 1:length(cb.angy)]
        D0=[[sin(cb.angx[i]),sin(cb.angy[j]),sqrt(1-sin(cb.angx[i])^2-sin(cb.angy[j])^2)] for i in 1:length(cb.angx),j in 1:length(cb.angy)]
    else 
        error("First four attributes only allowed in one of two configurations for bundle: \n 
            1. x and y are Type{AbstractArray{T,1} where T<:Real}, angx and angy are Type{T} where T<:Real \n
            2. x and y are Type{T} where T<:Real, angx and angy are Type{AbstractArray{T,1} where T<:Real} ")
    end
    return R0,D0 
end

"""
    bundle_as_array(x, y, angx, angy, zpos)

Composite type representing the origination of a set of rays. First four arguments are 
Type `AbstractMatrix{T} where T<:Real`; last argument is a scalar for the z position.

When a variable is declared this type with specific parameters and called 
as a function with no arguments, it returns a tuple (position, direction) each 
entry of which is an output array of vectors corresponding to the input information.

__NOTE__ if types of entries in inputs are known to all be the same, use `bundle_as_array_fast` with same syntax for increased speed.

# Examples
```julia-repl
julia> eval_bundle = bundle_as_array([0.5 1.], [0.5 1.], [0. 0.], [0. 0.], 0.)
bundle_as_array{Float64}([0.5 1.0], [0.5 1.0], [0.0 0.0], [0.0 0.0], 0.0)

julia> eval_bundle()
([[0.5000000000000002, 0.5000000000000002, 2.220446049250313e-16] [1.0000000000000002, 1.0000000000000002, 2.220446049250313e-16]], [[0.0, 0.0, 1.0] [0.0, 0.0, 1.0]])
```
"""
struct bundle_as_array{ x<:AbstractMatrix{T} where T<:Real,
                        y<:AbstractMatrix{T} where T<:Real,
                        angx<:AbstractMatrix{T} where T<:Real,
                        angy<:AbstractMatrix{T} where T<:Real,
                        zpos<:Real}
    x::x
    y::y
    angx::angx
    angy::angy
    zpos::zpos
end
struct bundle_as_array_fast{T,V}
    x::T
    y::T
    angx::T# these are direction tangents.
    angy::T
    zpos::V
end
function (cb::Union{bundle_as_array,bundle_as_array_fast})()

    if length(size(cb.zpos))==2
        R0=[[cb.x[i,j].+eps(),cb.y[i,j].+eps(),cb.zpos[i,j].+eps()] for i in 1:size(cb.x)[1],j in 1:size(cb.y)[2]]
    else
        R0=[[cb.x[i,j].+eps(),cb.y[i,j].+eps(),cb.zpos.+eps()] for i in 1:size(cb.x)[1],j in 1:size(cb.y)[2]]
    end
    D0=[[sin(cb.angx[i,j]),sin(cb.angy[i,j]),sqrt(1-sin(cb.angx[i,j])^2-sin(cb.angy[i,j])^2)] for i in 1:size(cb.x)[1],j in 1:size(cb.y)[2]]
    return R0,D0 
end

"""
    bundle_as_array_big(x, y, z, Dx, Dy, Dz)

Composite type representing the origination of a set of rays. All arguments are Type `AbstractArray{T,2} where T<:Real`.

__NOTE__ if all input arrays are guaranteed to have the same type, then use `bundle_as_array_big_fast` with same syntax for increased speed.

"""
struct bundle_as_array_big{ x<:AbstractMatrix{T} where T<:Real,
                            y<:AbstractMatrix{T} where T<:Real,
                            z<:AbstractMatrix{T} where T<:Real,
                            Dx<:AbstractMatrix{T} where T<:Real,
                            Dy<:AbstractMatrix{T} where T<:Real,
                            Dz<:AbstractMatrix{T} where T<:Real}
    x::x
    y::y
    z::z
    Dx::Dx# these are direction cosines.
    Dy::Dy
    Dz::Dz
end
struct bundle_as_array_big_fast{T}
    x::T
    y::T
    z::T
    Dx::T# these are direction cosines.
    Dy::T
    Dz::T
end
function (cb::Union{bundle_as_array_big,bundle_as_array_big_fast})()
    R0=[[cb.x[i,j].+eps(),cb.y[i,j].+eps(),cb.z[i,j].+eps()] for i in 1:size(cb.x)[1],j in 1:size(cb.y)[2]]
    D0=[[cb.Dx[i,j],cb.Dy[i,j],cb.Dz[i,j]] for i in 1:size(cb.x)[1],j in 1:size(cb.y)[2]]
    return R0,D0 
end

""" 
    trace(xr, yr, zr, nsr)

Composite type holding position and accumulated optical path for parsing. Each is Type `Matrix{Vector{T}} where T<:Real`, 
vector holding information along propagation and the matrix mapping to the chosen input bundle subscripting.
"""
struct trace{T}
    xr::T
    yr::T
    zr::T
    nsr::T
end

""" 
    bigtrace(xr, yr, zr, Dxr, Dyr, Dzr, nsr)

Composite type holding position, direction, and accumulated optical path for parsing. Each is Type `Matrix{Vector{T}} where T<:Real`, 
vector holding information along propagation and the matrix mapping to the chosen input bundle subscripting.

"""
struct bigtrace{T}
    xr::T
    yr::T
    zr::T
    Dxr::T
    Dyr::T
    Dzr::T
    nsr::T
end

"""
    bigtrace_to_bundle(tracething::bigtrace, pos_prop::Int)

Takes a `bigtrace` object and integer position `pos_prop` and returns a bundle object.

"""
function bigtrace_to_bundle(bT::bigtrace, pos_prop::Int)
    x=[(bT.xr[i,j][pos_prop] +eps()) for i in 1:size(bT.xr)[1], j in 1:size(bT.xr)[2]]
    y=[(bT.yr[i,j][pos_prop] +eps()) for i in 1:size(bT.xr)[1], j in 1:size(bT.xr)[2]]
    z=[(bT.zr[i,j][pos_prop] +eps()) for i in 1:size(bT.xr)[1], j in 1:size(bT.xr)[2]]
    Dx=[bT.Dxr[i,j][pos_prop] for i in 1:size(bT.xr)[1], j in 1:size(bT.xr)[2]]
    Dy=[bT.Dyr[i,j][pos_prop] for i in 1:size(bT.xr)[1], j in 1:size(bT.xr)[2]]
    Dz=[bT.Dzr[i,j][pos_prop] for i in 1:size(bT.xr)[1], j in 1:size(bT.xr)[2]]
    return bundle_as_array_big(x,y,z,Dx,Dy,Dz)
end

"""
    trace_extract_terminus(tracething::Union{trace, bigtrace}, pos_prop::Int; coord="all")

Takes a `trace` or `bigtrace` object and integer position pos_prop and returns matrices of coordinates according to the bundle.

Coordinates are `x,y,z`, direction cosines are `Dx,Dy,Dz` and accumulated optical path is `ns`.

Optional argument coord can be set "all" to return all information as Type `Matrix`.  This will be passed out as
- x, y, z, ns for `trace` object
- x, y, z, Dx, Dy, Dz, ns for `bigtrace` object

Or optional argument can be set to a different string to only pass out one `Matrix` of information.
- "x", "y", "z", "ns" for `trace` object
- "x", "y", "z", "Dx", "Dy", "Dz", "ns" for `bigtrace` object

"""
function trace_extract_terminus(bT::bigtrace, pos_prop;coord="all")
    x=[bT.xr[i,j][pos_prop] for i in 1:size(bT.xr)[1], j in 1:size(bT.xr)[2]]
    y=[bT.yr[i,j][pos_prop] for i in 1:size(bT.xr)[1], j in 1:size(bT.xr)[2]]
    z=[bT.zr[i,j][pos_prop] for i in 1:size(bT.xr)[1], j in 1:size(bT.xr)[2]]
    Dx=[bT.Dxr[i,j][pos_prop] for i in 1:size(bT.xr)[1], j in 1:size(bT.xr)[2]]
    Dy=[bT.Dyr[i,j][pos_prop] for i in 1:size(bT.xr)[1], j in 1:size(bT.xr)[2]]
    Dz=[bT.Dzr[i,j][pos_prop] for i in 1:size(bT.xr)[1], j in 1:size(bT.xr)[2]]
    ns=[bT.nsr[i,j][pos_prop] for i in 1:size(bT.xr)[1], j in 1:size(bT.xr)[2]]
    if coord=="all"
        return x,y,z,Dx,Dy,Dz
    elseif coord=="x"
        return x
    elseif coord=="y"
        return y
    elseif coord=="z"
        return z
    elseif coord=="Dx"
        return Dx
    elseif coord=="Dy"
        return Dy
    elseif coord=="Dz"
        return Dz
    elseif coord=="ns"
        return ns
    end
end
function trace_extract_terminus(bT::trace, pos_prop;coord="all")
    x=[bT.xr[i,j][pos_prop] for i in 1:size(bT.xr)[1], j in 1:size(bT.xr)[2]]
    y=[bT.yr[i,j][pos_prop] for i in 1:size(bT.xr)[1], j in 1:size(bT.xr)[2]]
    z=[bT.zr[i,j][pos_prop] for i in 1:size(bT.xr)[1], j in 1:size(bT.xr)[2]]
    ns=[bT.nsr[i,j][pos_prop] for i in 1:size(bT.xr)[1], j in 1:size(bT.xr)[2]]
    if coord=="all"
        return x,y,z
    elseif coord=="x"
        return x
    elseif coord=="y"
        return y
    elseif coord=="z"
        return z
    elseif coord=="ns"
        return ns
    end
end

"""
    trace_extract_ray(tracething::Union{trace, bigtrace}, i::Int, j::Int; coord = "all")

Takes a `trace` or `bigtrace` object and integers i,j and returns the ray at position in the traced bundle array i,j.

Coordinates are `x, y, z`, direction cosines are `Dx, Dy, Dz` and accumulated optical path is `ns`.

Optional argument coord can be set "all" to return all information as Matrices.  This will be passed out as
- x, y, z, ns for `trace` object
- x, y, z, Dx, Dy, Dz, ns for `bigtrace` object

Or optional argument can be set to a different string to only pass out one matrix of information.
- "x", "y", "z", "ns" for `trace` object
- "x", "y", "z", "Dx", "Dy", "Dz", "ns" for `bigtrace` object


"""
function trace_extract_ray(bT::bigtrace,i::Int,j::Int;coord="all")
    x=bT.xr[i,j]
    y=bT.yr[i,j]
    z=bT.zr[i,j]
    Dx=bT.Dxr[i,j]
    Dy=bT.Dyr[i,j]
    Dz=bT.Dzr[i,j]
    ns=bT.nsr[i,j]
    if coord=="all"
        return x,y,z,Dx,Dy,Dz
    elseif coord=="x"
        return x
    elseif coord=="y"
        return y
    elseif coord=="z"
        return z
    elseif coord=="Dx"
        return Dx
    elseif coord=="Dy"
        return Dy
    elseif coord=="Dz"
        return Dz
    elseif coord=="ns"
        return ns
    end
end
function trace_extract_ray(bT::trace,i,j;coord="all")
    x=bT.xr[i,j]
    y=bT.yr[i,j]
    z=bT.zr[i,j]
    ns=bT.nsr[i,j]
    if coord=="all"
        return x,y,z
    elseif coord=="x"
        return x
    elseif coord=="y"
        return y
    elseif coord=="z"
        return z
    elseif coord=="ns"
        return ns
    end
end

struct diffSf{impSurf,volPhase}
    impSurf::impSurf 
    volPhase::volPhase # volPhase(x,y,z) is the definition here. grating period is defined as the deriv phase fn divided by 2pi
end
function (rf::diffSf)(R0,D0,m,lambda;isReorient=true)
    nvf=normalVec(rf.impSurf)
    N=nvf(R0)
    vpc(R)=rf.volPhase(R[1],R[2],R[3])
    gradPhase(R)=ForwardDiff.gradient(vpc,R)
    kvec=cross(cross(N,gradPhase(R0)),N)
    D1=D0+kvec*lambda/2/pi*m
    D1./sqrt(dot(D1,D1))
end

function propSingleGRIN2(R,D,indexRef,indexRefGrad;dt=0.1)
    # Appl.Opt. Sharma et al 1982
    n = indexRef(R)
    ∇n = indexRefGrad(R)
    TT = n.*D
    A = dt.*∇n .*n
    RB = R .+ dt./2. .*TT.+dt./8 .*A
    nB = indexRef(RB)
    ∇nB = indexRefGrad(RB)
    B = dt.*∇nB .*nB
    RC = R .+ dt .* TT .+ dt./ 2. .* B
    nC = indexRef(RC)
    ∇nC = indexRefGrad(RC)
    C = dt.*∇nC .*nC
    Tnew = TT .+ 1. ./ 6. .* (A .+ 4. .* B .+ C)
    Ro=R .+ dt*(TT .+ 1. ./ 6. .* (A.+ 2. .* B))
    Dro=Tnew./n
    return Ro,Dro
end


CoeffLikeOption=Union{Vector{Vector{T}},Nothing} where T
VectorLikeOption=Union{Vector,Nothing}
NumLikeOption=Union{Real,Nothing}

"""
    opticalstack_instance = opticalstack( ... )

    ... = opticalstack_instance( ... )

For surface reflection and refraction use this method to define an instance of the `opticalstack` type

    opticalstack(coeffslist, surfslist, nlist)

For gradient index (and surface reflection and refraction) use this method to define an instance of the `opticalstack` type

    opticalstack(coeffslist, surfslist, nlist, ncoeffslist, grindt)

For grating/hologram/metasurface (and surface reflection and refraction) use this method to define an instance of the `opticalstack` type

    opticalstack(coeffslist, surfslist, nlist, diffractiveslist, diffcoeffslist, diffmlist, wavelength)

To use all features together use this method to define an instance of the `opticalstack` type

    opticalstack(coeffslist, surfslist, nlist, diffractiveslist, diffcoeffslist, diffmlist, wavelength, ncoeffslist, grindt)

Arguments to create an instance of `opticalstack` are 

- `coeffslist` is a `Vector` set of coefficients for a corresponding surface; each set of coefficients is itself also a `Vector`.
- `surfslist` is a `Vector` of functions callable on `(x, y, coeffs)` to return a z-zurface position.
- `nlist` is a `Vector` of index values or index functions corresponding to the medium before the corresponding surface.
    - Vector values T<:Real work universally
    - Vector input of mixed type, `Real` and functions `n` where `n(x, y, z, coeffs) = ...`, only works for declarations of `opticalstack` which define `ncoeffslist` and `grindt` as well.

- `ncoeffslist::Vector{Vector{T}} where T<:Real` is a `Vector` set of coefficients for a corresponding gradient index; each set of coefficients is itself also a `Vector`.
- `grindt::Real` is the Sharma ray increment in the gradient index medium. 
    - This step parameter scales the physical quantity ``OP/n²`` (optical path length and index of refraction) at every step, and raytracing is by the algorithm of Sharma et al. 1982.

-  `diffractiveslist` is a `Vector` of zeros and functions. Where a vector entry is not zero, the entry must be a valid phase function, e.g. ``ϕ`` where ``ϕ(x, y, z, coeffs) = ...`` sampled by the surface at a corresponding position to make the diffractive.
- `diffcoeffslist::Vector{Vector{T}} where T<:Real` is a `Vector` set of coefficients for a corresponding diffractive; each set of coefficients is itself also a `Vector`.
- `diffmlist::Vector{Int}` must have the integer corresponding to the diffraction order in the vector position corresponding to the position of the diffractive surface.
- `wavelength::Real` is the wavelength in system units simulated for the definition of the diffractive.

When an instance is declared, `opticalstack_instance = opticalstack( ... )`, the declared instance is called on a bundle instance `bundlething` and optional arguments as

    opticalstack_instance(bundlething; rend = nothing, color = :red, issurfplot = true, numsurfpoints = 20, halfdomain = nothing, xdom = nothing, ydom = nothing, plobj = nothing, isbigtrace = false, plotplane_at = 0., init_prop_guesses=nothing)

- When `rend` is `nothing` and `isbigtrace = false` the output is a `trace` type.
- When `rend` is `nothing` and `isbigtrace = true` the output is a `bigtrace` type.
- When `rend` is not `nothing`, `rend` must be `"YZ"` or `"XZ"` or `"3Dcirc"` or `"3Dsq"`. The output will be a plot object compatible with Plots.jl
- When `rend` is not `nothing` and `issurfplot = false` the surfaces being modeled will not be plotted.

- The domain plot can be scaled globally by entering a vector of desired plot extent into `xdom` or `ydom`. This is preferred over using `numsurfpoints` and `halfdomain.`
- `plotplane_at` is the value at which a 2D plot is sliced in the non-plotted dimension.
- `init_prop_guesses::Vector{T} where T<:Real` is an option to declare a first guess in numerical propagation preceding surface in the corresponding vector entry. 
    - Defaults to zero.
    - Not recommended to use a nonzero value when traversing a gradient index function. 

__Note__ that the data scaling of the output of gradient index analysis is not only with the number of rays but also with the resolution of propagation.  Assuming
that propagation occurs over some position difference ``Δz``, this data scaling is APPROXIMATELY ``Δz``/`dt`. As expected, `dt` trades precision for computation time.

__Note__ for gradient index surfaces the output ray coordinates are stacked according to propagation needed for the "longest" gradient index ray. Therefore a returned `trace` or `bigtrace` will not be easily indexed according to its surface, but all propagation positions will properly align when entering and exiting the gradient index.
Shorter rays will replicate their terminated values until the longest ray terminates. Then all rays will exit the gradient index together at the appropriate surface.

__Note__ for diffractive surfaces, the optical path contribution of the corresponding surface is that due to its preceding propagation PLUS that due to the diffraction phase function.

__Note__ that the index of refraction must be equal on both sides of the diffractive.  To work with cases where the diffractive sits on a substrate of different index of refraction, two surfaces should be collocated in sequence, one for reflection/refraction and one for diffraction.

# Examples
```julia-repl
julia> opticalstack_test = opticalstack([[1., 6.], [2., -6.], [4., 8.]], [spherical, spherical, zernike], [1., 1.5, 1.])
opticalstack{Float64}([[1.0, 6.0], [2.0, -6.0], [4.0, 8.0]], UnionAll[spherical, spherical, zernike], [1.0, 1.5, 1.0])

julia> bundle_test = bundle([0.5, 1.], [0.5, 1.], 0., 0., 0.)
bundle([0.5, 1.0], [0.5, 1.0], 0.0, 0.0, 0.0)

julia> tracething = opticalstack_test(bundle_test)
trace{Float64}([[0.5, 0.5, 0.4742890889441142, -0.024645519233983726] ... ; ...])
```

"""
struct opticalstack
    coeffslist::CoeffLikeOption
    surfslist::VectorLikeOption
    nlist::VectorLikeOption
    diffractiveslist::VectorLikeOption
    diffcoeffslist::CoeffLikeOption
    diffmlist::VectorLikeOption
    wavelength::NumLikeOption
    ncoeffslist::CoeffLikeOption
    grindt::NumLikeOption
end
opticalstack(
    coeffslist::CoeffLikeOption,
    surfslist::VectorLikeOption,
    nlist::VectorLikeOption)=opticalstack(
        coeffslist,
        surfslist,
        nlist,
        nothing,nothing,nothing,nothing,nothing,nothing)
opticalstack(
    coeffslist::CoeffLikeOption,
    surfslist::VectorLikeOption,
    nlist::VectorLikeOption,
    diffractiveslist::VectorLikeOption,
    diffcoeffslist::CoeffLikeOption,
    diffmlist::VectorLikeOption,
    wavelength::NumLikeOption) = opticalstack(
        coeffslist,
        surfslist,
        nlist,
        diffractiveslist,
        diffcoeffslist,
        diffmlist,
        wavelength,
        nothing,nothing)
opticalstack(
    coeffslist::CoeffLikeOption,
    surfslist::VectorLikeOption,
    nlist::VectorLikeOption,
    ncoeffslist::CoeffLikeOption,
    grindt::NumLikeOption) = opticalstack(
        coeffslist,
        surfslist,
        nlist,
        nothing,nothing,nothing,nothing,
        ncoeffslist,
        grindt)
function (optSta::opticalstack)(bundleEv;rend=nothing,color=:red,issurfplot=true,numsurfpoints=20,halfdomain=nothing,xdom=nothing,ydom=nothing,plobj=nothing,isbigtrace=false,plotplane_at=0.,init_prop_guesses=nothing)
    R0list,D0list=bundleEv()
    Rmaster=[R0list]
    Dmaster=[D0list]
    nsmaster=[zeros(size(D0list))]
    nlistApp=reduce(vcat,[optSta.nlist,optSta.nlist[end]]) # add final index to list (equivalent to making sure image surface is immersed.)
    
    if optSta.diffractiveslist===nothing
        diffList=zeros(length(optSta.coeffslist))
    else 
        diffList=optSta.diffractiveslist
    end
    for i in 1:length(optSta.coeffslist)
        impz1=impz(optSta.coeffslist[i],optSta.surfslist[i])
        if init_prop_guesses===nothing
            gi1=getIntersect_preInit(impz1,0.)
        else 
            gi1=getIntersect_preInit(impz1,init_prop_guesses[i])
        end

        
        if !(typeof(nlistApp[i])<:Real)
            nu(R)=nlistApp[i](R[1],R[2],R[3],optSta.ncoeffslist[i])
            ∇nu(R)=ForwardDiff.gradient(nu,R)
            fntraR(R,D)=propSingleGRIN2(R,D,nu,∇nu;dt=optSta.grindt)[1]
            fntraD(R,D)=propSingleGRIN2(R,D,nu,∇nu;dt=optSta.grindt)[2]
            sc=gi1.(Rmaster[end],Dmaster[end])
            while any(sc.>=0)
                Rn=fntraR.(Rmaster[end],Dmaster[end])
                Dn=fntraD.(Rmaster[end],Dmaster[end])
                sc=gi1.(Rn,Dn)
                Rn=Rn.*0.5.*(sign.(sc).+1.) .+ Rmaster[end].*0.5.*(sign.(-sc).+1.)
                Dn=Dn.*0.5.*(sign.(sc).+1.) .+ Dmaster[end].*0.5.*(sign.(-sc).+1.)
                nsmaster=reduce(vcat,[nsmaster,[nsmaster[i].+optSta.grindt.*nu.(Rmaster[end]).^2]])
                Rmaster=reduce(vcat,[Rmaster,[Rn]])
                Dmaster=reduce(vcat,[Dmaster,[Dn]])
            end
        end
        s=gi1.(Rmaster[end],Dmaster[end])
        R1=Rmaster[end].+s.*Dmaster[end]
        
        if diffList[i]!=0.
            if nlistApp[i]!=nlistApp[i+1]; @error "refractive indices must be equal on either side of a diffractive surface function";end
            if nlistApp[i]!=1.0; @info "diffractive is not immersed in n=1, be careful that phase function definition compensates for this"; end
            diffUse(x,y,z)=diffList[i](x,y,z,optSta.diffcoeffslist[i])
            ds1=diffSf(impz1,diffUse)
            D1=ds1.(R1,Dmaster[end],optSta.diffmlist[i],optSta.wavelength) 
        else
            rs1=refrSf(impz1)
            if typeof(nlistApp[i+1])<:Real
                n2=nlistApp[i+1]
            else
                if length(optSta.ncoeffslist)==i; grC=optSta.ncoeffslist[i]; else; grC=optSta.ncoeffslist[i+1]; end
                nR2(R)=nlistApp[i+1](R[1],R[2],R[3],grC)
                n2=nR2.(R1)
            end
            if typeof(nlistApp[i])<:Real
                n1=nlistApp[i]
            else
                nR1(R)=nlistApp[i](R[1],R[2],R[3],optSta.ncoeffslist[i])
                n1=nR1.(R1)
            end
            D1=rs1.(R1,Dmaster[end],n1,n2) 
        end
        if diffList[i]!=0.
            diffUseR(R)=diffUse(R[1],R[2],R[3])
            nsmaster=reduce(vcat,[nsmaster,[nsmaster[end].+nlistApp[i].*(s).+ 1. ./ 2 ./ pi .* diffUseR.(Rmaster[end])]])
        else
            nsmaster=reduce(vcat,[nsmaster,[nsmaster[end].+n1.*(s)]]) # this syntax is more diff-able
        end
        Rmaster=reduce(vcat,[Rmaster,[R1]])
        Dmaster=reduce(vcat,[Dmaster,[D1]])
    end
    Rm=Rmaster
    nsm=nsmaster
    xr=[[Rm[i][j,k][1] for i in 1:length(Rm)] for j in 1:size(Rm[1])[1],k in 1:size(Rm[1])[2]]
    yr=[[Rm[i][j,k][2] for i in 1:length(Rm)] for j in 1:size(Rm[1])[1],k in 1:size(Rm[1])[2]]
    zr=[[Rm[i][j,k][3] for i in 1:length(Rm)] for j in 1:size(Rm[1])[1],k in 1:size(Rm[1])[2]]
    if isbigtrace
        Dm=Dmaster
        Dxr=[[Dm[i][j,k][1] for i in 1:length(Rm)] for j in 1:size(Dm[1])[1],k in 1:size(Dm[1])[2]]
        Dyr=[[Dm[i][j,k][2] for i in 1:length(Rm)] for j in 1:size(Dm[1])[1],k in 1:size(Dm[1])[2]]
        Dzr=[[Dm[i][j,k][3] for i in 1:length(Rm)] for j in 1:size(Dm[1])[1],k in 1:size(Dm[1])[2]]
    end
    nsr=[[nsm[i][j,k] for i in 1:length(nsm)] for j in 1:size(nsm[1])[1],k in 1:size(nsm[1])[2]]
    
    halfdomainl=0
    if halfdomain===nothing
        if rend=="3Dcirc"
            halfdomainv_mx=maximum.(Rm).*1.41
            halfdomainv_mn=minimum.(Rm).*1.41
        else 
            halfdomainv_mx=maximum.(Rm)
            halfdomainv_mn=minimum.(Rm)
        end
        halfdomainl=[max(halfdomainv_mx[j][1],halfdomainv_mx[j][2],halfdomainv_mn[j][1],halfdomainv_mn[j][2]) for j in 1:length(halfdomainv_mx)]
    end 

    if rend!==nothing
        if plobj===nothing 
            plobj=Plots.plot() 
        end
        if rend=="YZ"
            if issurfplot
                for ij in 1:length(optSta.surfslist)
                    zmprof(y)=optSta.surfslist[ij](plotplane_at,y,optSta.coeffslist[ij])
                    if halfdomainl!=0; halfdomain=halfdomainl[ij+1]; end
                    res=2*halfdomain/(numsurfpoints-1)
                    if ydom===nothing; ydom=-halfdomain:res:halfdomain; end
                    Plots.plot!(plobj,zmprof.(ydom),ydom,color=:black,aspect_ratio=:equal)
                end
            end
            Plots.plot!(plobj,zr[:],yr[:],color=color,legend=:none,xlabel="z",ylabel="y")
        elseif rend=="XZ"
            if issurfplot
                for ij in 1:length(optSta.surfslist)
                    zmprof(x)=optSta.surfslist[ij](x,plotplane_at,optSta.coeffslist[ij])
                    if halfdomainl!=0; halfdomain=halfdomainl[ij+1]; end
                    res=2*halfdomain/(numsurfpoints-1)
                    if xdom===nothing; xdom=-halfdomain:res:halfdomain; end
                    Plots.plot!(plobj,zmprof.(xdom),xdom,color=:black,aspect_ratio=:equal)
                end
            end
            Plots.plot!(plobj,zr[:],xr[:],color=color,legend=:none,xlabel="z",ylabel="x")
        elseif rend=="3Dsq" 
            if issurfplot
                for ij in 1:length(optSta.surfslist)
                    if halfdomainl!=0; halfdomain=halfdomainl[ij+1]; end
                    res=2*halfdomain/(numsurfpoints-1)
                    if xdom===nothing; xdom=-halfdomain:res:halfdomain; end
                    if ydom===nothing; ydom=-halfdomain:res:halfdomain; end
                    x=xdom'.*ones(length(ydom))
                    y=ones(length(xdom))'.*ydom
                    zmprof(x,y)=optSta.surfslist[ij](x,y,optSta.coeffslist[ij])
                    Plots.surface!(plobj,xdom,ydom,zmprof.(x,y),color=:blue,alpha=0.7,legend=:none,colorbar=nothing)
                end
            end
            Plots.plot!(plobj,xr[:],yr[:],zr[:],color=color,legend=:none,xlabel="x",ylabel="y",zlabel="z")
        elseif rend=="3Dcirc"     
            
            if issurfplot
                for ij in 1:length(optSta.surfslist)
                    if halfdomainl!=0; halfdomain=halfdomainl[ij+1]; end
                    res=2*halfdomain/(numsurfpoints-1)
                    if xdom===nothing; xdom=-halfdomain:res:halfdomain; end
                    if ydom===nothing; ydom=-halfdomain:res:halfdomain; end
                    x=xdom'.*ones(length(ydom))
                    y=ones(length(xdom))'.*ydom
                    zmprof(x,y)=optSta.surfslist[ij](x,y,optSta.coeffslist[ij])
                    zp=zmprof.(x,y)
                    rad=halfdomain
                    zp[x.^2 .+ y.^2 .>rad.^2].=NaN
                    Plots.surface!(plobj,xdom,ydom,zp,color=:blue,alpha=0.7,legend=:none,colorbar=nothing)
                end
            end
            Plots.plot!(plobj,xr[:],yr[:],zr[:],color=color,legend=:none,xlabel="x",ylabel="y",zlabel="z")
        end
        return plobj
    else
        if !isbigtrace
            return trace(xr,yr,zr,nsr)
        else 
            return bigtrace(xr,yr,zr,Dxr,Dyr,Dzr,nsr)
        end
    end
end


"""
    rms_spot(tracething::Union{trace, bigtrace}; isunsqrt = false, pos = nothing, iszincluded = true)

Evaluates the rms spot size for a trace intersecting a surface. 

Pass in the result of a trace, and calcuate the rms spot based on the ray values for the subscript `pos`. If 
pos===nothing then the final values of the trace are used (i.e. the final surface is taken to be the evaluation surface).
Option `isunsqrt` removes the square root from the calculation (returns the sum of squared error instead). The z coordinate variation is included in the calculation when `iszincluded == true`.

# Examples
```julia-repl
julia> opticalstack_test = opticalstack([[1., 6.], [2., -6.], [4., 8.]], [spherical, spherical, zernike], [1., 1.5, 1.])
opticalstack{Float64}([[1.0, 6.0], [2.0, -6.0], [4.0, 8.0]], UnionAll[spherical, spherical, zernike], [1.0, 1.5, 1.0])

julia> bundle_test = bundle([0.5, 1.], [0.5, 1.], 0., 0., 0.)
bundle([0.5, 1.0], [0.5, 1.0], 0.0, 0.0, 0.0)

julia> tracething = opticalstack_test(bundle_test)
trace{Float64}([[0.5, 0.5, 0.4742890889441142, -0.024645519233983726] [0.5, 0.5, 0.47749522247588994, -0.04245690108918665]; [1.0, 1.0, 0.9549904449517799, -0.0849138021783733] [1.0, 1.0, 0.9616736859140215, -0.1253539680837904]], [[0.5, 0.5, 0.4742890889441142, -0.024645519233983726] [1.0, 1.0, 0.9549904449517799, -0.0849138021783733]; [0.5, 0.5, 0.47749522247588994, -0.04245690108918665] [1.0, 1.0, 0.9616736859140215, -0.1253539680837904]], [[0.0, 1.0418113625438212, 1.9623895629300963, 7.999999280260749] [0.0, 1.105086411777687, 1.904234300100513, 7.999999273328101]; [0.0, 1.105086411777687, 1.904234300100513, 7.999999273328101] [0.0, 1.169047547709123, 1.8438309411610754, 7.999999266127459]], [[0.0, 1.0418113625438212, 2.9446610481344417, 9.023361779810676] [0.0, 1.105086411777687, 2.8587256850369682, 9.064376469832027]; [0.0, 1.105086411777687, 2.8587256850369682, 9.064376469832027] [0.0, 1.169047547709123, 2.7690064459109145, 9.114214613215701]])

julia> rms_spot(tracething)
0.1125476814135891
```
"""
function rms_spot(tracething;isunsqrt=false,pos=nothing,iszincluded=true)
    # extract the last position of ray information x,y, then compute sse for that data
    xru=tracething.xr[:]
    yru=tracething.yr[:]
    zru=tracething.zr[:]
    if pos===nothing
        xe=[xru[j][end] for j =1:length(xru)]
        ye=[yru[j][end] for j =1:length(yru)]
        ze=[zru[j][end] for j =1:length(zru)]
    else
        xe=[xru[j][pos] for j =1:length(xru)]
        ye=[yru[j][pos] for j =1:length(yru)]
        ze=[zru[j][pos] for j =1:length(zru)]
    end
    xe=xe.-sum(xe)/length(xe)
    ye=ye.-sum(ye)/length(ye)
    ze=ze.-sum(ze)/length(ze)
    if iszincluded
        mse=sum(xe.^2 .+ ye.^2 .+ ze.^2)/length(xe)
    else
        mse=sum(xe.^2 .+ ye.^2)/length(xe)
    end
    if isunsqrt
        return mse
    else
        return sqrt(mse)
    end
end

"""
    splot(tracething::Union{trace, bigtrace}; pos=nothing, pl_splot = nothing, color = :red, markersize = 2)

Plots the spot of a particular trace or bigtrace object at position `pos`. If `pos === nothing`, the last position in the trace is 
used. Can optionally pass in a handle to existing plot pl_splot, and change the marker size from default `markersize = 2` and the plot color from default `color = :red` to any 
of the conventional options supported by Plots.jl.

"""
function splot(tracething;pos=nothing,pl_splot=nothing,color=:red,markersize=2)
    
    if pl_splot===nothing 
        pl_splot=Plots.plot() 
    end
    xru=tracething.xr[:]
    yru=tracething.yr[:]
    if pos===nothing
        xe=[xru[j][end] for j =1:length(xru)]
        ye=[yru[j][end] for j =1:length(yru)]
    else
        xe=[xru[j][pos] for j =1:length(xru)]
        ye=[yru[j][pos] for j =1:length(yru)]
    end  
    if size(tracething.xr)[1]==1
        Plots.scatter!(pl_splot,xe,ye,legend=:none,xlims=[-0.1,0.1],xlabel="x",ylabel="y",color=color,markersize=markersize)
    elseif size(tracething.xr)[2]==1
        Plots.scatter!(pl_splot,xe,ye,legend=:none,ylims=[-0.1,0.1],xlabel="x",ylabel="y",color=color,markersize=markersize)
    else
        Plots.scatter!(pl_splot,xe,ye,legend=:none,aspect_ratio=:equal,xlabel="x",ylabel="y",color=color,markersize=markersize)
    end
    return pl_splot
end


"""
    rac(tracething::Union{trace, bigtrace}, pos_in, pos_out; pl_rac = nothing, color = :red, axistype = "y", pmajoraxis = 2)

Plots the chosen physical coordinate `axistype=y` or `x` at an end position `pos_out` of a particular `trace` or `bigtrace` object as a function of the same coordinate at input position `pos_in`. If `pos===nothing`, the last position in the trace is 
used. The coordinate plots are stitched together in lines based on the chosen matrix subscript in the `trace` or `bigtrace`.  This matrix subscript is set by the optional argument `pmajoraxis = 2` or `1`. Can optionally pass in a handle to existing plot pl_rac, and change the plot color from default `color=:red` to any 
of the conventional options supported by Plots.jl.

NOTE that when tracing a ray fan (a bundle with one dimension of the array having length 1), the `pmajoraxis` must be set to the correct number to span the nonunitary length of the array. Else the plot may not appear to show anything.

"""
function rac(tracething,pos_in,pos_out;pl_rac=nothing,color=:red,axistype="y",pmajoraxis=2)
    
    if pl_rac===nothing 
        pl_rac=Plots.plot() 
    end
    if axistype=="x"
        ru=tracething.xr
    elseif axistype=="y"
        ru=tracething.yr
    end
    
    re_i=[ru[i,j][pos_in] for i =1:size(ru)[1],j =1:size(ru)[2]]
    re_o=[ru[i,j][pos_out] for i =1:size(ru)[1],j =1:size(ru)[2]]
 
    if pmajoraxis==1
        Plots.plot!(pl_rac,re_i,re_o,legend=:none,xlabel=axistype*" pos. "*string(pos_in),ylabel=axistype*" pos. "*string(pos_out),color=color)
    elseif pmajoraxis==2
        Plots.plot!(pl_rac,re_i',re_o',legend=:none,xlabel=axistype*" pos. "*string(pos_in),ylabel=axistype*" pos. "*string(pos_out),color=color)
    end
    return pl_rac
end


end # module