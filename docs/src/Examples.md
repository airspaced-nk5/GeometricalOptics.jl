# Examples (2D and 3D)

## Vectorizing over bundle angles

The [Quickstart tutorial](@ref) showed tracing a single bundle through a system.  Now suppose there is a need to trace multiple types of bundles through a system all at once.  Consider a system of two parabolic surfaces which have bundles incident at different angles. One could simply do the procedure in the tutorial three different ways, but this would be somewhat verbose.  Instead this can be done more efficiently by defining functions over the variables which are of interest and then calling them on vectors of the different variable values.

Below a setup is implemented very similarly to the tutorial, but a function `bund_ang(angle)` is defined to originate a bundle of fixed spatial origin but variable y-direction-tangent angle.

```@example 1;continued=true
using GeometricalOptics, Plots

surfslist = [aspherical, aspherical, zplane]
coeffslist = [[5., -6., -1.], [1., 2., -1.], [8.]]      # two parabolic surfaces
nlist = [1., -1., 1.]                               # index sign change means reflection

optsta_instance = opticalstack(coeffslist, surfslist, nlist)

bund_ang(angle) = bundle([0.], (-0.5:0.5:0.5).+1, 0., angle, 0.)
``` 

Then vectors to vary for the system are defined. Three different angles are of interest, and for visiblity it is best to plot in different colors and only plot surfaces along with a single set of rays (as the system does not change otherwise). Then a Plots.jl plot pane is initialized so that it can be passed into all of the vectorized cases of analysis and they can be plotted together.

```@example 1; continued=true
angles = [0., 0.028, 0.04]                    # a list of three angles
colors = [:red, :green, :blue]                # each of which is colored differently
issurfplotvec = [true, false, false]          # only plot surfaces along with first bundle
plot_together_in = Plots.plot()               # must put this into the system plotter to get all angles.
```

Finally a function `stack_render` is defined to accept all of the varying parameters, and then the vectors are all passed into the function to return a vector of plots. The vector of plots are linked to the same original `Plots.plot` instance, so calling the last element of the vector of plots will show all three instances together.

```@example 1
stack_render(angle, pc, isp) = optsta_instance(bund_ang(angle); rend = "YZ", color = pc, numsurfpoints = 100, plobj = plot_together_in, issurfplot = isp)

plot_together_out_vec = stack_render.(angles, colors, issurfplotvec)
plot_together_out_vec[3]                    # the last entry in the vector of plots has all bundles
```

And an analogous procedure generates all of the spot plots together. A function is defined to set a now square bundle for a variable y-angle, and then another function plots all of the spots together on the final plane.

```@example 1 
splot_init = Plots.plot()

bund2_ang(angle) = bundle((-0.5:0.2:0.5), (-0.5:0.2:0.5).+1, 0., angle, 0.)

function optsta_render_splot(angle, pc, isp)
    traceout = optsta_instance(bund2_ang(angle))
    splot(traceout, pl_splot = splot_init, color = pc)
end

splot_vec = optsta_render_splot.(angles, colors, issurfplotvec)
splot_vec[3]                                # the last entry in the vector of plots has all bundles
```

## Advantageous start conditions and fixing a broken raytrace

Broadly, there are few ways to break the raytrace analysis. One way is to attempt to violate physics or basic geometry. Another way is due to extremely steep surface or ray angles with respect to z.  This should be tested. Still another ays is to not provide proper start conditions for numerical raytrace methods.  There is a utility to address this latter case.

A few other things should be tested:
- Rays must be known to converge on the surface.
    - Don't choose a function that rolls off before the ray can hit it; that cannot be fixed.
- The surface function must be designed to return a number at the ray incident position.
    - Example: built-in [`zernike`](@ref) is defined to return `NaN` when outside the normalization radius, so it is __designed to fail in certain cases__.
    
Using a modification of the above example, it is possible to "break" the raytrace.  Take the same system for analysis as above, but change the parameters of the first asphere to be positioned and focused differently.  Note that a few rays traced through the first surface are traced erratically.  This is due to improper convergence of numerical raytracing.

```@example 2
using GeometricalOptics, Plots

surfslist = [aspherical, aspherical, zplane]
coeffslist = [[4., -4., -1.], [1., 2., -1.], [8.]]      # two parabolic surfaces
nlist = [1., -1., 1.]                               # index sign change means reflection

optsta_instance=opticalstack(coeffslist, surfslist, nlist)

bund_ang(angle) = bundle([0.], (-0.5:0.5:0.5).+1, 0., angle, 0.)

angles = [0., 0.028, 0.04]                     # a list of three angles
colors = [:red, :green, :blue]                 # each of which is colored differently
issurfplotvec = [true, false, false]           # only plot surfaces along with first bundle
plot_together_in = Plots.plot()                # must put this into the system plotter to get all angles.

stack_render(angle,pc,isp) = optsta_instance(bund_ang(angle); rend = "YZ",color = pc,ydom = -2:0.1:2, plobj = plot_together_in, issurfplot = isp)

plot_together_out_vec = stack_render.(angles, colors, issurfplotvec)
plot_together_out_vec[3]                    # the last entry in the vector of plots has all bundles
title!("broke the trace :(")
```

The spots show the erratic rays as well:

```@example 2
splot_init = Plots.plot()

bund2_ang(angle) = bundle((-0.5:0.2:0.5), (-0.5:0.2:0.5).+1, 0., angle, 0.)

function optsta_render_splot(angle, pc, isp)
    traceout = optsta_instance(bund2_ang(angle))
    splot(traceout, pl_splot = splot_init, color = pc)
end

splot_vec = optsta_render_splot.(angles, colors, issurfplotvec)
splot_vec[3]                                # the last entry in the vector of plots has all bundles
title!("broke the trace :(")
```

This case can be fixed by using the optional argument `init_prop_guesses` to an instance of `opticalstack`.  The optional argument can be passed in as a vector of same length as the number of surfaces, or as `nothing`.  When `nothing`, all rays are first guessed to propagate zero distance to the next surface.  When passed in as a vector, each value of the vector is the initial guess of signed propagation distance to the corresponding surface. Because the error arises in finding intersection of the rays with the second surface, the initial guess for the second surface is altered, and in particular made negative (due to reflection) and given magnitude close to the distance expected for the ray.  Then all of the rays in a bundle will be traced separately, only with newly specified starting guesses for distances traversed between surfaces.

The issue with the raytrace is rectified by adding the vector of values to provide a better initial guess for the negative propagation to the second surface:

```@example 3
using GeometricalOptics, Plots

surfslist = [aspherical, aspherical, zplane]
coeffslist = [[4., -4., -1.], [1., 2., -1.], [8.]]      # two parabolic surfaces
nlist = [1., -1., 1.]                               # index sign change means reflection

optsta_instance = opticalstack(coeffslist, surfslist, nlist)

bund_ang(angle) = bundle([0.], (-0.5:0.5:0.5).+1, 0., angle, 0.)

angles = [0., 0.028, 0.04]                      # a list of three angles
colors = [:red, :green, :blue]                 # each of which is colored differently
issurfplotvec = [true, false, false]            # only plot surfaces along with first bundle
plot_together_in = Plots.plot()               # must put this into the system plotter to get all angles.

init_prop_guesses_choice = [0., -2., 0.]

stack_render(angle, pc, isp) = optsta_instance(bund_ang(angle); rend = "YZ", color = pc, ydom = -2:0.1:2, plobj = plot_together_in, issurfplot = isp, init_prop_guesses = init_prop_guesses_choice)

plot_together_out_vec = stack_render.(angles, colors, issurfplotvec)
plot_together_out_vec[3]                    # the last entry in the vector of plots has all bundles
title!("fixed the trace :)")
```

The spots behave now too:

```@example 3
splot_init = Plots.plot()

bund2_ang(angle) = bundle((-0.5:0.2:0.5), (-0.5:0.2:0.5).+1, 0., angle, 0.)

function optsta_render_splot(angle, pc, isp)
    traceout = optsta_instance(bund2_ang(angle); init_prop_guesses = init_prop_guesses_choice)
    splot(traceout, pl_splot = splot_init, color = pc)
end

splot_vec = optsta_render_splot.(angles, colors, issurfplotvec)
splot_vec[3]                                # the last entry in the vector of plots has all bundles
title!("fixed the trace :)")
```


## Vectorizing over index, 3D analysis 

The same sort of thing that was done for different bundle angles can also be done for different system configurations.  Let's examine a few variations on this, plotting in 3D instead of 2D and using a custom bundle array defined to be a circular bundle footprint.  First define the parameters of the system that are fixed. The bundle is fixed with a circular footprint defined by `bundle_as_array`.  The array arguments are defined such that one array subscript corresponds to a radial distance from the central axis, and another subscript corresponds to an azimuthal angle around the axis.

```@example 4; continued=true
using GeometricalOptics, Plots

surfslist = [aspherical, aspherical, zplane]
coeffslist = [[1., 6., -5.], [2., -6., -5.], [8.]]
# index list is no longer up here, as the material is varied in modeling.

rdom = 0.01:0.2: 1.
thdom = 0.:pi/10:2pi
r = rdom' .* ones(length(thdom))
th = ones(length(rdom))' .* thdom
bund_instance = bundle_as_array(r .* cos.(th), r .* sin.(th), zeros(size(r)), zeros(size(r)), 0.)
```

Then define the variable system parameters.  The selected index values for the second material are in the vector `nsel`, and all three of these will be tested; results will be plotted in different colors on the same plot with surfaces only plotted alongside the first case.

```@example 4; continued=true

nsel = [1.46, 1.5, 1.56]                # a list of three indices,
colors = [:red, :green, :blue]         # each of which is colored differently
issurfplotvec = [true, false, false]    # only plot surfaces with one of the bundles.
plot_lens_init = Plots.plot()       # must put this into the system plotter to get all cases.

```

Now make a function to vary the system over the desired values, and then call over the vectors defined above:

```@example 4
function optsta_render(n_u, pc, isp)
    nlist_chr = [1., n_u, 1.]
    optsta_instance = opticalstack(coeffslist, surfslist, nlist_chr)
    optsta_instance(bund_instance; rend = "3Dcirc", color = pc, numsurfpoints = 100, plobj = plot_lens_init, issurfplot = isp)
end

plot_lens_vec = optsta_render.(nsel, colors, issurfplotvec)
plot_lens_vec[3]                # only need to plot the last of the three elements in the vector; it caught all cases.
```

Using the same configuration, a YZ plot can be made.  This plot, with traces of a circular bundle, is a projection of all of the rays onto the YZ plane.  Thus the YZ plot (properly) appears to have irregular rays because the rays span 3D space and are projected into a 2D plane.  This is why it may be advantageous to plot fans (2D bundles) for lens viewing instead of 3D bundles, as in previous examples.

```@example 4
plot_lens_init3d = Plots.plot()  
function optsta_render(n_u, pc, isp)
    nlist_chr = [1., n_u, 1.]
    optsta_instance = opticalstack(coeffslist, surfslist, nlist_chr)
    optsta_instance(bund_instance; rend = "YZ", color = pc, numsurfpoints = 100, plobj = plot_lens_init3d, issurfplot = isp)
end
plot_lens_vec3d = optsta_render.(nsel, colors, issurfplotvec)
plot_lens_vec3d[3]                 # only need to plot the last of the three elements in the vector; it caught all cases.
```

Finally, an analogous process is used to get the spot plot.  Different index of refraction corresponds to different size of the spot for this configuration.

```@example 4
splot_init = Plots.plot()
function optsta_render_splot(n_u, pc, isp)
    nlist_chr = [1., n_u, 1.]
    optsta_instance = opticalstack(coeffslist, surfslist, nlist_chr)
    traceout = optsta_instance(bund_instance)
    splot(traceout, pl_splot = splot_init, color = pc, markersize = 3)
end
splot_vec = optsta_render_splot.(nsel, colors, issurfplotvec)
splot_vec[3]
```


## [Sideloading a new surface function: caustics in the swimming pool](@id custSurf)

Let's try to model the lines that appear when the sun strikes the surface waves on water, as in a swimming pool.  

Suppose we can approximate the surface waves in a swimming pool by a sinusoid on a polynomial (a better model might exist for this, but this one is qualitatively sufficient). This function can be defined as `wavy_water` below and sideloaded into the raytrace analysis. The function is first defined over arguments `x,y,coeffs`, where `x,y` are lateral coordinates and `coeffs::Vector{T} where T<: Real` has entries which are decided by the user.  Then put the function list and put a corresponding coefficient vector with desired values into `coeffslist`.

The different densities of the ray spots at the analysis plane generates a wavy line pattern as can be seen at the bottom of a swimming pool.

```@example 5
using GeometricalOptics

wavy_water(x, y, coeffs) = coeffs[1] + coeffs[2] * cos(coeffs[3] * (x^2 * 0.05 + 0.2*x +y))

surfslist = [wavy_water, zplane]
coeffslist = [[1., 0.1, 3.], [8.]]
nlist = [1., 1.33]
optsta_instance = opticalstack(coeffslist, surfslist, nlist)

bund_instance = bundle((-2:0.1:2), (-2:0.1:2), 0., 0., 0.)

traceout = optsta_instance(bund_instance)
plot_spot = splot(traceout)
```

To view the raytrace, let's use a less dense bundle than the one used for evaluation.  The rays and surfaces can be viewed from the side. Note that the side plot shows a projection of all rays in the bundle with various different angles.

```@example 5 
bund_instance = bundle((-2:0.5:2), (-2:0.5:2), 0., 0., 0.)
plot_sys = optsta_instance(bund_instance, rend = "YZ")
```

Or the rays and surfaces can be viewed in 3D. The bottom of the pool is at the top of the fixed view pane shown below.  This is better observed when the plot is viewed interactively - try it out!

```@example 5 
plot_sys2 = optsta_instance(bund_instance, rend = "3Dsq")
```

This example illustrates how a user defined surface of arbitrary shape can be directly used in analysis, accepting a vector of coefficients just as for built-in functions.

## [Advanced: Gradient index (GRIN) lens with flat surfaces](@id grin)

The index of refraction can vary continuously in a material.  Light can be focused with a gradient index distribution even in the absence of a curved surface.  This can be modeled by replacing the index entry with a quadratic GRIN function and passing a couple extra arguments to the `opticalstack` for the coefficients of the index function and the discrete ray sampling in the GRIN medium.

First the quadratic GRIN coefficient is calculated to give the approximate desired focal length for a given thickness of the GRIN, and the GRIN function for use is defined. It must be written as `n(x,y,z,coeffs)` where `coeffs` is the coefficient vector mapped into the GRIN function and `x,y,z` are the coordinates being modeled. 

```@example 6;continued=true
using GeometricalOptics 

f = 10.               # some quick physical calculation to get approximate focal length
t = 2.
n_r2 = -1 / 2 / f / t

n_custom(x, y, z, coeffs) = coeffs[1] + coeffs[2] * (x^2 + y^2)
``` 

Now the system `opticalstack` is set up almost exactly the same except now with two more arguments. `dt` sets the GRIN step size and `n_coeffslist` provides a list of coefficient lists accessible to an index function at a corresponding position in `nlist`.  In this way, `n_coeffslist` serves the GRIN analogous to how `coeffslist` serves surfaces in construction.

```@example 6;continued=true

surfslist = [zplane, zplane, zplane]
coeffslist = [[1.], [1 + t], [1 + t + f]]
nlist = [1., n_custom, 1.]

dt = 0.003                                # set the discrete step of the curved raytrace
n_coeffslist = [[0.], [1.5, n_r2], [0.]]     # only entries where values matter are those 
                                        #   corresponding to entries of nlist which are functions

optsta_instance = opticalstack(coeffslist, surfslist, nlist, n_coeffslist, dt)
``` 

Then a y-fan can be traced and the YZ cross-section plotted:

```@example 6
bund_instance = bundle([0.], (-4:0.5:4), 0., 0., 0.)
plot_lens = optsta_instance(bund_instance; rend = "YZ")
```

And the ray position at input and output in the y direction can be compared:

```@example 6
traceout = optsta_instance(bund_instance)
plot_rac = rac(traceout, 1, 4)
```

"Whoa," you say, "that doesn't look right, my input and output ray coordinates are not the same size as this plot would suggest!"  That would be correct.  The reason is that the GRIN implementation traces a bunch of tiny ray segments and stores them in the `trace` or `bigtrace` output.  Unlike for the basic systems with reflective or refractive surfaces, now the propagation subscripts of the rays in the bundle are no longer linked to the number of surfaces, but to the number of surfaces *plus* the number of ray increments in the GRIN.

Let's extract a ray and examine its subscript length.  Just arbitrarily pick the first ray in the trace here:

```@example 6
zvalues_ray11 = trace_extract_ray(traceout, 1, 1; coord = "z")
length(zvalues_ray11)
```

And now it is apparent that the ray had to be traced in a total of 598 segments to meet the GRIN ray increment resolution of `dt=0.003`.  Using this information identifies the propagation position corresponding to the intersection of the ray with the last plane.

```@example 6
traceout = optsta_instance(bund_instance)
plot_rac = rac(traceout, 1, 598)
```

Note that when using this in vectorized form or with dense ray bundles, the GRIN implementation will slow the raytrace for small chosen values of `dt`.  This is a trade of speed and accuracy that should be made carefully.

## [Advanced: Grating/hologram/metasurface lens](@id grating) 
The procedure to raytrace through gratings, holograms, and metasurfaces is the same, as all three contribute a phase perturbation to the wavefront local to the ray intersection with a surface.  

It is worth noting that although the raytrace procedure is the same for gratings, holograms, and metasurfaces, their real use can be quite different.

First a few calculations for the coefficient of the quadratic diffractive are carried out, and the wavelength in system units is defined as 0.00055.  When the system units are taken to be millimeters, this wavelength corresponds to green light.  All physical quantities in this package must be entered with consistent units; choice of this consistent unit scale is flexible.

```@example 7; continued=true
using GeometricalOptics

lambda_wl = 0.00055 # note: must keep units consistent!
f = 7.
quad_opl_cff = -1 / 2 / f
quad_phase_cff = 2 * pi / lambda_wl * quad_opl_cff

phase_custom(x, y, z, coeffs) = coeffs[1] * (x^2 + y^2)
```

The diffraction grating is defined by a volume phase function in a list of functions (`difflist` below) sampled by a surface in `surflist`. A volume phase function `phase_custom(x,y,z,coeffs)` is defined with a quadratic profile to accept the calculated coefficient value in the coefficient vector `coeffs`.  The volume phase function ``\phi(x,y,z)`` is used to get the real phase contribution ``\phi_{surf}(x,y)`` on an arbitrarily shaped surface ``z(x,y)`` by

```math
\phi_{surf}(x, y) = \phi(x, y, z(x, y)).
```

There is no restriction on the shape of the surface. This enables a phase function on a freeform surface to be modeled.

Where there are no diffraction gratings, the entry in `difflist` should be zero.  A vector of integers (`diffmlist` below) and a vector of vectors (`diffcoeffslist` below) are passed in alongside the diffraction function vector to pass in the desired diffraction order and the vector of coefficients for all diffractive gratings.  The entries to these later vectors are only important when the entry in `difflist` is a function and nonzero.

Important note: the index of refraction before and after the diffractive surface must be equal. If the desired physical system has two different materials on either side of the diffractive, then the system should be modeled as a sequence of a diffractive surface and an identically located and defined reflective or refractive surface.

```@example 7

difflist = [phase_custom, 0.]
diffmlist = [1, 0]
diffcoeffslist = [[quad_phase_cff], [1.]]

surfslist = [zplane, zplane]
coeffslist = [[1.], [8.]]
nlist = [1., 1.]

optsta_instance = opticalstack(coeffslist, surfslist, nlist, difflist, diffcoeffslist, diffmlist, lambda_wl)

bund_instance = bundle([0.], (-4:0.5:4), 0., 0., 0.)
plot_lens = optsta_instance(bund_instance; rend = "YZ")
```

The plot of the y coordinates of the YZ fan at input and output positions in propagation are then given

```@example 7
traceout = optsta_instance(bund_instance)
plot_rac = rac(traceout, 1, 3)
```

Noting the vertical scale, it is observed how raytraced diffractives theoretically promise tight or even "perfect" focusing.  The imperfect performance of this diffractive is appreciated only when the wave nature of light is considered, which is outside the scope of pure ray-based analysis.

