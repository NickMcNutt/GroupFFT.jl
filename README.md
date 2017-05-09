# GroupFFT.jl

[![Build Status](https://travis-ci.org/NickMcNutt/GroupFFT.jl.svg?branch=master)](https://travis-ci.org/NickMcNutt/GroupFFT.jl)
[![codecov.io](http://codecov.io/github/NickMcNutt/GroupFFT.jl/coverage.svg?branch=master)](http://codecov.io/github/NickMcNutt/GroupFFT.jl?branch=master)

## Overview

* A Julia package for performing [fast Fourier transforms](https://en.wikipedia.org/wiki/Fast_Fourier_transform) over arbitrary [compact groups](https://en.wikipedia.org/wiki/Compact_group)

* Currently supports the rotation group, [SO(3)](https://en.wikipedia.org/wiki/Rotation_group_SO(3)), and the rotation/reflection group, [O(3)](https://en.wikipedia.org/wiki/Orthogonal_group).

* Licensed under the [MIT License](https://opensource.org/licenses/MIT)

## How to install

In Julia:
```julia
Pkg.clone("https://github.com/NickMcNutt/GroupFFT.jl")
```

## Examples

Generate an element of the [orthogonal group](https://en.wikipedia.org/wiki/Orthogonal_group) O(3) using the Euler ZYZ parameterization:

```julia
using GroupFFT

α, β, γ = π/2, π/2, π/2
g = O3(α, β, γ)
```

Compute the unitary irreducible representations of element `g`. For the group O(3), these are [Wigner-D matrices](https://en.wikipedia.org/wiki/Wigner_D-matrix).

```julia
unitary_irreps = [unitary_irrep(g, i) for i in 0:2]
display.(round.(unitary_irreps, 3))
```
```julia
1×1 Array{Complex{Float64},2}:
 1.0+0.0im

3×3 Array{Complex{Float64},2}:
 -0.5+0.0im     0.0+0.707im   0.5+0.0im  
 -0.0-0.707im   0.0+0.0im     0.0-0.707im
  0.5+0.0im    -0.0+0.707im  -0.5-0.0im  

5×5 Array{Complex{Float64},2}:
   0.25-0.0im  -0.0-0.5im  -0.612+0.0im  0.0+0.5im    0.25+0.0im
    0.0+0.5im   0.5-0.0im     0.0+0.0im  0.5+0.0im     0.0-0.5im
 -0.612+0.0im  -0.0-0.0im    -0.5-0.0im  0.0-0.0im  -0.612-0.0im
   -0.0-0.5im   0.5+0.0im    -0.0+0.0im  0.5+0.0im    -0.0+0.5im
   0.25+0.0im  -0.0+0.5im  -0.612-0.0im  0.0-0.5im    0.25+0.0im
```

If we don't want complex numbers, we can generate real irreducible representations instead:

```julia
orthogonal_irreps = [orthogonal_irrep(g, i) for i in 0:2]
display.(round.(orthogonal_irreps, 3))
```

```julia
1×1 Array{Float64,2}:
 1.0

3×3 Array{Float64,2}:
 -1.0  -0.0  0.0
  0.0  -0.0  1.0
 -0.0   1.0  0.0

5×5 Array{Float64,2}:
 -0.0  -1.0  -0.0     0.0  -0.0  
 -1.0   0.0   0.0    -0.0   0.0  
  0.0  -0.0  -0.5    -0.0   0.866
 -0.0   0.0  -0.0     1.0   0.0  
  0.0  -0.0   0.866   0.0   0.5
```

## Functionality

GroupFFT.jl specifies an abstract base type `Group` from which all group types derive. Currently, support is available for the following groups:

```julia
OrthogonalGroup{T, 3}
SpecialOrthogonalGroup{T, 3}
SymmetricGroup{T, N}
```

For convenience, type aliases are provided for groups of low dimension:

```julia
O3{T}
SO3{T}
```

## Contributing

GroupFFT.jl aims to be a comprehensive package that provides base types and linear representations for all of the most widely used groups.
We welcome the support of anyone who wishes to contribute to this project.

#### Completed

* Base types and irreducible representations for O(3), SO(3), and S<sub>n</sub>

#### Todo

* Base types and representations for O(n), U(n), SU(n), SL(n), SE(n), PSL(n)
* Support for other kinds of group representations (standard, faithful, etc.)
* Support for fields other than the complex and real numbers
