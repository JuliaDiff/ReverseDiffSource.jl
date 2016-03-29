#########################################################################
#
#   Base rules
#
#########################################################################

# derivation neutral functions
ReverseDiffSource.@deriv_rule colon(x,y)   x     0.
ReverseDiffSource.@deriv_rule colon(x,y)   y     0.

ReverseDiffSource.@deriv_rule length(x)    x     0.

ReverseDiffSource.@deriv_rule size(x)      x     0.
ReverseDiffSource.@deriv_rule size(x,y)    x     0.
ReverseDiffSource.@deriv_rule size(x,y)    y     0.

ReverseDiffSource.@deriv_rule fill(x,y)    x     0.
ReverseDiffSource.@deriv_rule fill(x,y)    y     0.

ReverseDiffSource.@deriv_rule similar(x,y) x     0.
ReverseDiffSource.@deriv_rule similar(x,y) y     0.

ReverseDiffSource.@deriv_rule zeros(x)     x     0.

ReverseDiffSource.@deriv_rule ones(x)      x     0.

ReverseDiffSource.@deriv_rule cell(x)      x     0.

ReverseDiffSource.@deriv_rule sign(x)      x     0.

ReverseDiffSource.@deriv_rule reverse(x)   x     0.

#  tuple  TODO : specialized macro for this kind of function
ReverseDiffSource.@deriv_rule tuple(x)        x     ds[1]
ReverseDiffSource.@deriv_rule tuple(x,y)      x     ds[1]
ReverseDiffSource.@deriv_rule tuple(x,y)      y     ds[2]
ReverseDiffSource.@deriv_rule tuple(x,y,z)    x     ds[1]
ReverseDiffSource.@deriv_rule tuple(x,y,z)    y     ds[2]
ReverseDiffSource.@deriv_rule tuple(x,y,z)    z     ds[3]
ReverseDiffSource.@deriv_rule tuple(x,y,z,t)  x     ds[1]
ReverseDiffSource.@deriv_rule tuple(x,y,z,t)  y     ds[2]
ReverseDiffSource.@deriv_rule tuple(x,y,z,t)  x     ds[3]
ReverseDiffSource.@deriv_rule tuple(x,y,z,t)  t     ds[4]

#  vcat
ReverseDiffSource.@deriv_rule vcat(x)        x     ds[1]
ReverseDiffSource.@deriv_rule vcat(x,y)      x     ds[1]
ReverseDiffSource.@deriv_rule vcat(x,y)      y     ds[2]
ReverseDiffSource.@deriv_rule vcat(x,y,z)    x     ds[1]
ReverseDiffSource.@deriv_rule vcat(x,y,z)    y     ds[2]
ReverseDiffSource.@deriv_rule vcat(x,y,z)    z     ds[3]
ReverseDiffSource.@deriv_rule vcat(x,y,z,t)  x     ds[1]
ReverseDiffSource.@deriv_rule vcat(x,y,z,t)  y     ds[2]
ReverseDiffSource.@deriv_rule vcat(x,y,z,t)  x     ds[3]
ReverseDiffSource.@deriv_rule vcat(x,y,z,t)  t     ds[4]

# reshape
ReverseDiffSource.@deriv_rule reshape(x::AbstractArray, a, b)        x    reshape(ds, size(x))
ReverseDiffSource.@deriv_rule reshape(x::AbstractArray, a, b)        a    0.
ReverseDiffSource.@deriv_rule reshape(x::AbstractArray, a, b)        b    0.
ReverseDiffSource.@deriv_rule reshape(x::AbstractArray, d::Tuple)    x    reshape(ds, size(x))
ReverseDiffSource.@deriv_rule reshape(x::AbstractArray, d::Tuple)    d    0.


# square root
ReverseDiffSource.@deriv_rule sqrt(x::Real)              x     0.5 * x ^ (-0.5) * ds
ReverseDiffSource.@deriv_rule sqrt(x::AbstractVector)    x     0.5 .* x .^ (-0.5) .* ds

# addition
ReverseDiffSource.@deriv_rule +(x::Real         , y::Real )            x     ds
ReverseDiffSource.@deriv_rule +(x::AbstractArray, y::AbstractArray)    x     ds
ReverseDiffSource.@deriv_rule +(x::Real         , y::Real )            y     ds
ReverseDiffSource.@deriv_rule +(x::AbstractArray, y::AbstractArray)    y     ds

ReverseDiffSource.@deriv_rule +(x::Real         , y::AbstractArray)    x     sum(ds)
ReverseDiffSource.@deriv_rule +(x::AbstractArray, y       )            x     ds
ReverseDiffSource.@deriv_rule +(x::AbstractArray, y::Real )            y     sum(ds)
ReverseDiffSource.@deriv_rule +(x               , y::AbstractArray)    y     ds

# dot addition
ReverseDiffSource.@deriv_rule .+(x::Real         , y::Real )            x     ds
ReverseDiffSource.@deriv_rule .+(x::Real         , y::AbstractArray)    x     sum(ds)
ReverseDiffSource.@deriv_rule .+(x::AbstractArray, y       )            x     ds
ReverseDiffSource.@deriv_rule .+(x::Real         , y::Real )            y     ds
ReverseDiffSource.@deriv_rule .+(x::AbstractArray, y::Real )            y     sum(ds)
ReverseDiffSource.@deriv_rule .+(x               , y::AbstractArray)    y     ds

# unary substraction
ReverseDiffSource.@deriv_rule -(x::Real )                              x     -ds
ReverseDiffSource.@deriv_rule -(x::AbstractArray)                      x     -ds

# binary substraction
ReverseDiffSource.@deriv_rule -(x::Real         , y::Real )            x     ds
ReverseDiffSource.@deriv_rule -(x::AbstractArray, y::AbstractArray)    x     ds
ReverseDiffSource.@deriv_rule -(x::Real         , y::Real )            y     -ds
ReverseDiffSource.@deriv_rule -(x::AbstractArray, y::AbstractArray)    y     -ds

# dot binary substraction
ReverseDiffSource.@deriv_rule .-(x::Real         , y::Real )            x     ds
ReverseDiffSource.@deriv_rule .-(x::Real         , y::AbstractArray)    x     sum(ds)
ReverseDiffSource.@deriv_rule .-(x::AbstractArray, y       )            x     ds
ReverseDiffSource.@deriv_rule .-(x::Real         , y::Real )            y     -ds
ReverseDiffSource.@deriv_rule .-(x::AbstractArray, y::Real )            y     -sum(ds)
ReverseDiffSource.@deriv_rule .-(x               , y::AbstractArray)    y     -ds

# sum()
ReverseDiffSource.@deriv_rule sum(x::Real )                            x     ds
ReverseDiffSource.@deriv_rule sum(x::AbstractArray)                    x     ones(size(x)).*ds

# dot()
ReverseDiffSource.@deriv_rule dot(x::Real         , y::Real )          x     y * ds
ReverseDiffSource.@deriv_rule dot(x::Real         , y::Real )          y     x * ds

ReverseDiffSource.@deriv_rule dot(x::AbstractArray, y::AbstractArray)  x     y.*ds
ReverseDiffSource.@deriv_rule dot(x::AbstractArray, y::AbstractArray)  y     x.*ds

# log() and exp()
ReverseDiffSource.@deriv_rule log(x::Real )                            x     ds / x
ReverseDiffSource.@deriv_rule log(x::AbstractArray)                    x     ds ./ x

ReverseDiffSource.@deriv_rule exp(x::Real )                            x     exp(x) * ds 
ReverseDiffSource.@deriv_rule exp(x::AbstractArray)                    x     exp(x) .* ds

# sin() and cos()
ReverseDiffSource.@deriv_rule sin(x::Real )                            x     cos(x) * ds
ReverseDiffSource.@deriv_rule sin(x::AbstractArray)                    x     cos(x) .* ds

# func = colon ; args= [(:x,Any), (:y,Any)] ; dv = :x ; diff = 0.
# deriv_rule(sin, [(:x, Any)], :x, :(cos(x)*ds))

ReverseDiffSource.@deriv_rule cos(x::Real )                            x     -sin(x) * ds
ReverseDiffSource.@deriv_rule cos(x::AbstractArray)                    x     -sin(x) .* ds

# abs, max(), min()
ReverseDiffSource.@deriv_rule abs(x::Real )                            x     sign(x) * ds
ReverseDiffSource.@deriv_rule abs(x::AbstractArray)                    x     sign(x) .* ds

ReverseDiffSource.@deriv_rule max(x::Real         , y::Real )          x     (x > y) * ds
ReverseDiffSource.@deriv_rule max(x::Real         , y::AbstractArray)  x     sum((x .> y) .* ds)
ReverseDiffSource.@deriv_rule max(x::AbstractArray, y::Real )          x     (x .> y) .* ds
ReverseDiffSource.@deriv_rule max(x::AbstractArray, y::AbstractArray)  x     (x .> y) .* ds

ReverseDiffSource.@deriv_rule max(x::Real         , y::Real )          y     (x < y) * ds
ReverseDiffSource.@deriv_rule max(x::Real         , y::AbstractArray)  y     (x .< y) .* ds
ReverseDiffSource.@deriv_rule max(x::AbstractArray, y::Real )          y     sum((x .< y) .* ds)
ReverseDiffSource.@deriv_rule max(x::AbstractArray, y::AbstractArray)  y     (x .< y) .* ds

ReverseDiffSource.@deriv_rule min(x::Real         , y::Real )          x     (x < y) * ds
ReverseDiffSource.@deriv_rule min(x::Real         , y::AbstractArray)  x     sum((x .< y) .* ds)
ReverseDiffSource.@deriv_rule min(x::AbstractArray, y::Real )          x     (x .< y) .* ds
ReverseDiffSource.@deriv_rule min(x::AbstractArray, y::AbstractArray)  x     (x .< y) .* ds

ReverseDiffSource.@deriv_rule min(x::Real         , y::Real )          y     (x > y) * ds
ReverseDiffSource.@deriv_rule min(x::Real         , y::AbstractArray)  y     (x .> y) .* ds
ReverseDiffSource.@deriv_rule min(x::AbstractArray, y::Real )          y     sum((x .> y) .* ds)
ReverseDiffSource.@deriv_rule min(x::AbstractArray, y::AbstractArray)  y     (x .> y) .* ds

# maximum, minimum
ReverseDiffSource.@deriv_rule maximum(x::Real         )     x     ds
ReverseDiffSource.@deriv_rule maximum(x::AbstractArray)     x     (x .== maximum(x)) .* ds

ReverseDiffSource.@deriv_rule minimum(x::Real         )     x     ds
ReverseDiffSource.@deriv_rule minimum(x::AbstractArray)     x     (x .== minimum(x)) .* ds


# multiplication
ReverseDiffSource.@deriv_rule *(x::Real         , y::Real )            x     y * ds
ReverseDiffSource.@deriv_rule *(x::Real         , y::AbstractArray)    x     sum(y .* ds)
ReverseDiffSource.@deriv_rule *(x::AbstractArray, y::Real )            x     y .* ds
# ReverseDiffSource.@deriv_rule *(x::AbstractArray, y::Vector)          x     gemm!('N', 'T', 1., ds, reshape(y, length(y), 1), 1., dx)  # reshape needed 
ReverseDiffSource.@deriv_rule *(x::AbstractArray, y::AbstractArray)    x     ds * y'
 
ReverseDiffSource.@deriv_rule *(x::Real         , y::Real )            y     x * ds
ReverseDiffSource.@deriv_rule *(x::Real         , y::AbstractArray)    y     x .* ds
ReverseDiffSource.@deriv_rule *(x::AbstractArray, y::Real )            y     sum(x .* ds)
# ReverseDiffSource.@deriv_rule *(x::AbstractArray, y::Vector)          y     gemm!('T', 'N', 1., x, reshape(ds, length(ds), 1), 1., dy)
ReverseDiffSource.@deriv_rule *(x::AbstractArray, y::AbstractArray)    y     x' * ds

# dot multiplication
ReverseDiffSource.@deriv_rule .*(x::Real         , y::Real )           x     y .* ds
ReverseDiffSource.@deriv_rule .*(x::Real         , y::AbstractArray)   x     sum(y .* ds)
ReverseDiffSource.@deriv_rule .*(x::AbstractArray, y::Real )           x     y .* ds
ReverseDiffSource.@deriv_rule .*(x::AbstractArray, y::AbstractArray)   x     y .* ds

ReverseDiffSource.@deriv_rule .*(x::Real         , y::Real )           y     x * ds
ReverseDiffSource.@deriv_rule .*(x::Real         , y::AbstractArray)   y     x .* ds
ReverseDiffSource.@deriv_rule .*(x::AbstractArray, y::Real )           y     sum(x .* ds)
ReverseDiffSource.@deriv_rule .*(x::AbstractArray, y::AbstractArray)   y     x .* ds

# power  (both args reals)
ReverseDiffSource.@deriv_rule ^(x::Real, y::Real)                      x     y * x ^ (y-1) * ds
ReverseDiffSource.@deriv_rule ^(x::Real, y::Real)                      y     log(x) * x ^ y * ds

# dot power
ReverseDiffSource.@deriv_rule .^(x::Real         , y::Real )           x     y * x ^ (y-1) * ds
ReverseDiffSource.@deriv_rule .^(x::Real         , y::AbstractArray)   x     sum(y .* x .^ (y-1) .* ds)
ReverseDiffSource.@deriv_rule .^(x::AbstractArray, y::Real )           x     y * x .^ (y-1) .* ds
ReverseDiffSource.@deriv_rule .^(x::AbstractArray, y::AbstractArray)   x     y .* x .^ (y-1) .* ds

ReverseDiffSource.@deriv_rule .^(x::Real         , y::Real )           y     log(x) * x ^ y * ds
ReverseDiffSource.@deriv_rule .^(x::AbstractArray, y::Real )           y     sum( log(x) .* x .^ y .* ds)
ReverseDiffSource.@deriv_rule .^(x::Real         , y::AbstractArray)   y     log(x) .* x .^ y .* ds
ReverseDiffSource.@deriv_rule .^(x::AbstractArray, y::AbstractArray)   y     log(x) .* x .^ y .* ds

# division
ReverseDiffSource.@deriv_rule /(x::Real          , y::Real )           x     ds / y
ReverseDiffSource.@deriv_rule /(x::AbstractArray , y::Real )           x     ds ./ y

ReverseDiffSource.@deriv_rule /(x::Real          , y::Real )           y     -x * ds / (y * y)
ReverseDiffSource.@deriv_rule /(x::AbstractArray , y::Real )           y     sum(-x .* ds) / (y * y)

# dot division
ReverseDiffSource.@deriv_rule ./(x::Real         , y::Real )           x     ds / y
ReverseDiffSource.@deriv_rule ./(x::Real         , y::AbstractArray)   x     sum(ds ./ y)
ReverseDiffSource.@deriv_rule ./(x::AbstractArray, y::Real )           x     ds ./ y
ReverseDiffSource.@deriv_rule ./(x::AbstractArray, y::AbstractArray)   x     ds ./ y

ReverseDiffSource.@deriv_rule ./(x::Real         , y::Real )           y     -x * ds / (y * y)
ReverseDiffSource.@deriv_rule ./(x::Real         , y::AbstractArray)   y     -x * ds ./ (y .* y)
ReverseDiffSource.@deriv_rule ./(x::AbstractArray, y::Real )           y     -sum(x .* ds) / (y * y)
ReverseDiffSource.@deriv_rule ./(x::AbstractArray, y::AbstractArray)   y     -x .* ds ./ (y .* y)

# mod2pi
ReverseDiffSource.@deriv_rule mod2pi(x::Real)                          x     ds

# transpose
ReverseDiffSource.@deriv_rule transpose(x::Real )                      x     ds
ReverseDiffSource.@deriv_rule transpose(x::AbstractArray)              x     transpose(ds)

