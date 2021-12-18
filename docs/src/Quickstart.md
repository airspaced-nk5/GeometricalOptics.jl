# Quickstart tutorial

A lens can easily be defined, plotted, and evaluated in ten lines of code as shown in this tutorial.  The first few lines are for definition of system geometry, then bundle definition defines the rays for analysis, and finally the evaluation of the optical stack on the bundle yields the trace and plots used for analysis.

## System geometry

The system geometry is defined by 

```@example 1; continued=true
using GeometricalOptics

funcsList = [spherical, spherical, zplane]
coeffsList = [[1., 6.], [2., -6.], [8.]]
nList = [1., 1.5, 1.]

optSta = opticalstack(coeffsList, funcsList, nList)
```

To understand how this works, note that each surface in `funcsList` has an associated `Vector` in `coeffsList` and value in `nList`. The positions in all three `Vector`s line up. The coefficient list entry is used for the corresponding surface, and the index of refraction entry is for the space preceding the corresponding surface. The "precedence" in the previous sentence is set by the sequence of the surfaces in `surfList` and therefore the order of traversal of the surfaces by the ray bundle. By contrast, surface precedence is __not__ explicitly implied by the definition of the surface geometry.  This is useful for defining different types of systems.

The `spherical` surface has its own coefficient `Vector` definition. It is [`spherical(x,y,coeffs)`](@ref) where `coeffs=[zpos,signed_radius]` with entries z-position and signed radius, accordingly.  This functionality extends to allow simple user defined functions: the user can define a new function `my_ytilt(x,y,coeffs)=coeffs[1]+coeffs[2]*y`, which creates a tilted plane with z position set by the first entry in `coeffs` and y tilt set by the second entry.  The user would then pass in `my_ytilt` with no arguments to the function list vector (`funcsList` above).  This is explored with a more complicated and fun surface in [Sideloading a new surface function: caustics in the swimming pool](@ref custSurf).

Here there are three surfaces.  System units are arbitrary and can be prescribed based on the needs of demonstration.  The light will propagate through air (index of refraction n=1.0) to a [`spherical`](@ref) surface with apex at z=1 and radius +6.0, and then refract into a medium of n=1.5.  Then the light will propagate through a glass (n=1.5) to a second surface with apex at z=2 and radius -6.0.  Finally the light will refract back out to air and propagate to a z-plane ([`zplane`](@ref)) at z=8.0.

## Ray geometry 

The least verbose method to trace a bundle of light is by using the square bundle definition [`bundle(x,y,angx,angy,zpos)`](@ref). Other bundles are possible; for example, [`bundle_as_array`](@ref) can make a circular bundle as shown in [Vectorizing over index, 3D analysis](@ref). In this construction, `x,y,zpos` are physical positions, and `angx,angy` are tilts in radians away from the z-axis in the x and y directions respectively (arguments of direction tangents).  By convention, the rays all originate from the same z position given by scalar `zpos`. Arguments `x,y` can be `Vector`s and the rest scalars to make the bundle of rays all parallel or collimated.  These `Vector`s need not be the same length; this sets the count of the rays in either direction. If `x::Vector{T} where T<:Real = [0.]` as below, the rays will all begin with coordinate x=0, and have y values given by the extent of the vector in y.

```@example 1;continued=true
test_bundle = bundle([0.], (-1:0.2:1), 0., 0., 0.)
```

Analogously, another bundle can be traced but have `angx,angy` as `Vector`s and all other arguments as scalars. This would correspond with all rays originating from a point in space with a particular array of angles.  In any event, the [`bundle`](@ref) represents an array of input rays.

## Evaluation and plotting

The output of the method on the instance of [`opticalstack`](@ref) is variable depending on keyword. 

With a bundle passed in and no render optional argument `rend` used the default output is a [`trace`](@ref).

Passed into [`rac`](@ref), the ray lateral coordinate at one position in trajectory can be plotted as a function of ray lateral coordinate at another position.  In this case, the position is related to the ray bundle information.  Index 1 corresponds to the set of ray information before striking any surfaces, and the last index corresponds to the ray information at the last requested surface.  Therefore when traveling through k surfaces, the number of positions along ray trajectory will be (k+1).  The three-surface system is evaluated at ray position four with respect to ray position one.

```@example 1
trace1 = optSta(test_bundle)
p_rac = rac(trace1, 1, 4)
```

The bundle strikes the final surface on a z-plane.  Root-mean-square ray deviation [`rms_spot`](@ref) with no optional keyword `pos` is evaluated at the last surface.

```@example 1
rms = rms_spot(trace1)
```
The [`trace`](@ref) can be evaluated by other means using the evaluation functions explained in [Bundle, surface, and trace options](@ref). The [`bigtrace`](@ref) type stores more information than the [`trace`](@ref) type and can be extracted by setting optional argument `isbigtrace=true` into the `opticalstack` call. 

If the render mode keyword `rend` is used, the output is a plot, here set as a 2D plot of the lens in the YZ plane:

```@example 1
p_lens = optSta(test_bundle; rend = "YZ")

```

Argument options for plotting are `rend="YZ"`, `rend="XZ"`, `rend="3Dcirc"`, `rend="3Dsq"`.  The global plot domain of surfaces can be adjusted by `xdom` and `ydom`.

The various plots can be plotted statically as shown in the docs or as __interactive plots__ by using the plotly backend of Plots with this package. Just add `using Plots; Plots.plotly()` and the plots will be interactive.

```@example 1
p_lens = optSta(test_bundle; rend = "YZ", ydom = -2:0.1:2)

```
Yes, this is the eleventh line of code for those of you keeping count, and I said ten lines. But it's not necessary to repeat the plotting operation. You get the point :)

# *Now, what next?* 

There are other methods and types in GeometricalOptics.jl which can be used, and this package can be interfaced with other packages in Julia.  Other more complicated demonstrations can be constructed as in [Examples (2D and 3D)](@ref). But the core functionality still remains the same as in this tutorial.  The [`opticalstack`](@ref) connects a bundle of rays and with a sequence of surfaces. Simply substitute different functions for surfaces or vectorize over lenses to look at different configurations depending on the behavior you wish to demonstrate.  Other bundle prescriptions are also possible.  The abstract vectorization of Julia allows a surprisingly intuitive extension of this simpler example to other practical demonstrations.
