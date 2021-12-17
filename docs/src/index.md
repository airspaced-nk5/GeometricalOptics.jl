# GeometricalOptics.jl



## Overview 

GeometricalOptics.jl is a simple, compact, and extensible tool for optical raytrace evaluation in the Julia Programming Language.  

The language of this package can be defined by the sequence:

Bundle -> Optical System (Stack) -> Trace

Where a Bundle represents a set of input rays, an Optical System or Stack receives and operates on the Bundle, and the returned object is a Trace, which contains the information about how the input Bundle traversed the Optical System.

This package was designed with simplicity and extensibility as high priorities (but, of course, the user community will determine whether this goal was met).  As is seen in the [Quickstart tutorial](@ref), an optical system can be defined and analyzed in very few lines of code and minimal package-specific language where it can be avoided.  Then [Examples (2D and 3D)](@ref) shows how the general framework extends and how other functionality can be added to the system under analysis.  Examples include [implementing a custom surface](@ref custSurf), [implementing a custom gradient index](@ref grin), and [implementing a custom grating/hologram/metasurface](@ref grating).

This tool is not a substitute for, but rather an accompaniment to, education in optics and photonics. Minimal jargon is used in the docs, though optical concepts will be apparent in spite of this.  The developer hopes this is a useful educational tool, and that it also serves as a good toolbox for exploring different functional forms of raytraced optical functions.

## Installation

Add GeometricalOptics.jl using `Pkg`:

```
using Pkg; Pkg.add("GeometricalOptics")
```

Though not strictly necessary for basic analyses, it is also highly recommended to add Plots.jl.  This will appear in some of the examples.


## Roadmap

The package will be periodically updated with the following objectives in mind, in order of decreasing priority:

1. Technical function bug fixes
2. Improved interoperability between functions and types of GeometricalOptics.jl
3. Improved consistency with base Julia (generalizing or narrowing type declarations, for example)
4. Improved intraoperability with other packages
5. Generalization of existing functions
6. New functions


## How to support

Support this project by
- Making contact - indicating support, comments, or issues on Github
- Citing the package
- Contributing to the package

## License

This package is provided for use under the MIT License.