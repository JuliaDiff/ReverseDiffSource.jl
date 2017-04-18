#########################################################################
#
#   Base rules
#
#########################################################################

# derivation neutral functions
@deriv_rule colon(x,y)   x     0.
@deriv_rule colon(x,y)   y     0.

@deriv_rule length(x)    x     0.

@deriv_rule size(x)      x     0.
@deriv_rule size(x,y)    x     0.
@deriv_rule size(x,y)    y     0.

@deriv_rule fill(x,y)    x     0.
@deriv_rule fill(x,y)    y     0.

@deriv_rule similar(x,y) x     0.
@deriv_rule similar(x,y) y     0.

@deriv_rule zeros(x)     x     0.

@deriv_rule ones(x)      x     0.

@deriv_rule cell(x)      x     0.

@deriv_rule sign(x)      x     0.

@deriv_rule reverse(x)   x     0.

#  tuple  TODO : specialized macro for this kind of function
@deriv_rule tuple(x)        x     ds[1]
@deriv_rule tuple(x,y)      x     ds[1]
@deriv_rule tuple(x,y)      y     ds[2]
@deriv_rule tuple(x,y,z)    x     ds[1]
@deriv_rule tuple(x,y,z)    y     ds[2]
@deriv_rule tuple(x,y,z)    z     ds[3]
@deriv_rule tuple(x,y,z,t)  x     ds[1]
@deriv_rule tuple(x,y,z,t)  y     ds[2]
@deriv_rule tuple(x,y,z,t)  z     ds[3]
@deriv_rule tuple(x,y,z,t)  t     ds[4]

#  vcat
@deriv_rule vcat(x)        x     ds[1]
@deriv_rule vcat(x,y)      x     ds[1]
@deriv_rule vcat(x,y)      y     ds[2]
@deriv_rule vcat(x,y,z)    x     ds[1]
@deriv_rule vcat(x,y,z)    y     ds[2]
@deriv_rule vcat(x,y,z)    z     ds[3]
@deriv_rule vcat(x,y,z,t)  x     ds[1]
@deriv_rule vcat(x,y,z,t)  y     ds[2]
@deriv_rule vcat(x,y,z,t)  z     ds[3]
@deriv_rule vcat(x,y,z,t)  t     ds[4]

# reshape
@deriv_rule reshape(x::AbstractArray, a, b)        x    reshape(ds, size(x))
@deriv_rule reshape(x::AbstractArray, a, b)        a    0.
@deriv_rule reshape(x::AbstractArray, a, b)        b    0.
@deriv_rule reshape(x::AbstractArray, d::Tuple)    x    reshape(ds, size(x))
@deriv_rule reshape(x::AbstractArray, d::Tuple)    d    0.


# square root
@deriv_rule sqrt(x::Real)              x     0.5 * x ^ (-0.5) * ds
@deriv_rule sqrt(x::AbstractVector)    x     0.5 .* x .^ (-0.5) .* ds

# addition
@deriv_rule +(x::Real         , y::Real )            x     ds
@deriv_rule +(x::AbstractArray, y::AbstractArray)    x     ds
@deriv_rule +(x::Real         , y::Real )            y     ds
@deriv_rule +(x::AbstractArray, y::AbstractArray)    y     ds

@deriv_rule +(x::Real         , y::AbstractArray)    x     sum(ds)
@deriv_rule +(x::AbstractArray, y       )            x     ds
@deriv_rule +(x::AbstractArray, y::Real )            y     sum(ds)
@deriv_rule +(x               , y::AbstractArray)    y     ds

# dot addition
@deriv_rule .+(x::Real         , y::Real )            x     ds
@deriv_rule .+(x::Real         , y::AbstractArray)    x     sum(ds)
@deriv_rule .+(x::AbstractArray, y       )            x     ds
@deriv_rule .+(x::Real         , y::Real )            y     ds
@deriv_rule .+(x::AbstractArray, y::Real )            y     sum(ds)
@deriv_rule .+(x               , y::AbstractArray)    y     ds

# unary substraction
@deriv_rule -(x::Real )                              x     -ds
@deriv_rule -(x::AbstractArray)                      x     -ds

# binary substraction
@deriv_rule -(x::Real         , y::Real )            x     ds
@deriv_rule -(x::AbstractArray, y::AbstractArray)    x     ds
@deriv_rule -(x::Real         , y::Real )            y     -ds
@deriv_rule -(x::AbstractArray, y::AbstractArray)    y     -ds

# dot binary substraction
@deriv_rule .-(x::Real         , y::Real )            x     ds
@deriv_rule .-(x::Real         , y::AbstractArray)    x     sum(ds)
@deriv_rule .-(x::AbstractArray, y       )            x     ds
@deriv_rule .-(x::Real         , y::Real )            y     -ds
@deriv_rule .-(x::AbstractArray, y::Real )            y     -sum(ds)
@deriv_rule .-(x               , y::AbstractArray)    y     -ds

# sum()
@deriv_rule sum(x::Real )                            x     ds
@deriv_rule sum(x::AbstractArray)                    x     ones(size(x)).*ds

# dot()
@deriv_rule dot(x::Real         , y::Real )          x     y * ds
@deriv_rule dot(x::Real         , y::Real )          y     x * ds

@deriv_rule dot(x::AbstractArray, y::AbstractArray)  x     y.*ds
@deriv_rule dot(x::AbstractArray, y::AbstractArray)  y     x.*ds

# log() and exp()
@deriv_rule log(x::Real )                            x     ds / x
@deriv_rule log(x::AbstractArray)                    x     ds ./ x

@deriv_rule exp(x::Real )                            x     exp(x) * ds
@deriv_rule exp(x::AbstractArray)                    x     exp(x) .* ds

@deriv_rule log1p(x::Real)                           x     ds  / (1 + x)
@deriv_rule log1p(x::AbstractArray)                  x     ds ./ (1 + x)

@deriv_rule expm1(x::Real)                           x     (1. + expm1(x))  * ds
@deriv_rule expm1(x::AbstractArray)                  x     (1. + expm1(x)) .* ds
# note : derivative uses expm1() and not exp() to reuse the
#   already calculated expm1()

# trig functions
@deriv_rule sin(x::Real )                            x     cos(x) * ds
@deriv_rule sin(x::AbstractArray)                    x     cos(x) .* ds

@deriv_rule cos(x::Real )                            x     -sin(x) * ds
@deriv_rule cos(x::AbstractArray)                    x     -sin(x) .* ds

@deriv_rule tan(x::Real )                            x     (1. + tan(x)  * tan(x))  * ds
@deriv_rule tan(x::AbstractArray)                    x     (1. + tan(x) .* tan(x)) .* ds

@deriv_rule sinh(x::Real )                           x     cosh(x) * ds
@deriv_rule sinh(x::AbstractArray)                   x     cosh(x) .* ds

@deriv_rule cosh(x::Real )                           x     sinh(x) * ds
@deriv_rule cosh(x::AbstractArray)                   x     sinh(x) .* ds

@deriv_rule tanh(x::Real )                           x     (1. - tanh(x)  * tanh(x))  * ds
@deriv_rule tanh(x::AbstractArray)                   x     (1. - tanh(x) .* tanh(x)) .* ds

@deriv_rule asin(x::Real )                           x     ds  / sqrt(1 - x *x)
@deriv_rule asin(x::AbstractArray)                   x     ds ./ sqrt(1 - x.*x)

@deriv_rule acos(x::Real )                           x     -ds  / sqrt(1 - x *x)
@deriv_rule acos(x::AbstractArray)                   x     -ds ./ sqrt(1 - x.*x)

@deriv_rule atan(x::Real )                           x     ds  / (1 + x *x)
@deriv_rule atan(x::AbstractArray)                   x     ds ./ (1 + x.*x)


# round, floor, ceil, trunc, mod2pi
@deriv_rule round(x::Real)                           x     0.
@deriv_rule round(x::AbstractArray)                  x     0.

@deriv_rule floor(x::Real)                           x     0.
@deriv_rule floor(x::AbstractArray)                  x     0.

@deriv_rule ceil(x::Real)                            x     0.
@deriv_rule ceil(x::AbstractArray)                   x     0.

@deriv_rule trunc(x::Real)                           x     0.
@deriv_rule trunc(x::AbstractArray)                  x     0.

@deriv_rule mod2pi(x::Real)                          x     ds


# abs, max(), min()
@deriv_rule abs(x::Real )                            x     sign(x) * ds
@deriv_rule abs(x::AbstractArray)                    x     sign(x) .* ds

@deriv_rule max(x::Real         , y::Real )          x     (x > y) * ds
@deriv_rule max(x::Real         , y::AbstractArray)  x     sum((x .> y) .* ds)
@deriv_rule max(x::AbstractArray, y::Real )          x     (x .> y) .* ds
@deriv_rule max(x::AbstractArray, y::AbstractArray)  x     (x .> y) .* ds

@deriv_rule max(x::Real         , y::Real )          y     (x < y) * ds
@deriv_rule max(x::Real         , y::AbstractArray)  y     (x .< y) .* ds
@deriv_rule max(x::AbstractArray, y::Real )          y     sum((x .< y) .* ds)
@deriv_rule max(x::AbstractArray, y::AbstractArray)  y     (x .< y) .* ds

@deriv_rule min(x::Real         , y::Real )          x     (x < y) * ds
@deriv_rule min(x::Real         , y::AbstractArray)  x     sum((x .< y) .* ds)
@deriv_rule min(x::AbstractArray, y::Real )          x     (x .< y) .* ds
@deriv_rule min(x::AbstractArray, y::AbstractArray)  x     (x .< y) .* ds

@deriv_rule min(x::Real         , y::Real )          y     (x > y) * ds
@deriv_rule min(x::Real         , y::AbstractArray)  y     (x .> y) .* ds
@deriv_rule min(x::AbstractArray, y::Real )          y     sum((x .> y) .* ds)
@deriv_rule min(x::AbstractArray, y::AbstractArray)  y     (x .> y) .* ds

# maximum, minimum
@deriv_rule maximum(x::Real         )     x     ds
@deriv_rule maximum(x::AbstractArray)     x     (x .== maximum(x)) .* ds

@deriv_rule minimum(x::Real         )     x     ds
@deriv_rule minimum(x::AbstractArray)     x     (x .== minimum(x)) .* ds


# multiplication
@deriv_rule *(x::Real         , y::Real )            x     y * ds
@deriv_rule *(x::Real         , y::AbstractArray)    x     sum(y .* ds)
@deriv_rule *(x::AbstractArray, y::Real )            x     y .* ds
@deriv_rule *(x::AbstractArray, y::AbstractArray)    x     ds * y'

@deriv_rule *(x::Real         , y::Real )            y     x * ds
@deriv_rule *(x::Real         , y::AbstractArray)    y     x .* ds
@deriv_rule *(x::AbstractArray, y::Real )            y     sum(x .* ds)
@deriv_rule *(x::AbstractArray, y::AbstractArray)    y     x' * ds

# dot multiplication
@deriv_rule .*(x::Real         , y::Real )           x     y .* ds
@deriv_rule .*(x::Real         , y::AbstractArray)   x     sum(y .* ds)
@deriv_rule .*(x::AbstractArray, y::Real )           x     y .* ds
@deriv_rule .*(x::AbstractArray, y::AbstractArray)   x     y .* ds

@deriv_rule .*(x::Real         , y::Real )           y     x * ds
@deriv_rule .*(x::Real         , y::AbstractArray)   y     x .* ds
@deriv_rule .*(x::AbstractArray, y::Real )           y     sum(x .* ds)
@deriv_rule .*(x::AbstractArray, y::AbstractArray)   y     x .* ds

# power  (both args reals)
@deriv_rule ^(x::Real, y::Real)                      x     y * x ^ (y-1) * ds
@deriv_rule ^(x::Real, y::Real)                      y     log(x) * x ^ y * ds

# dot power
@deriv_rule .^(x::Real         , y::Real )           x     y * x ^ (y-1) * ds
@deriv_rule .^(x::Real         , y::AbstractArray)   x     sum(y .* x .^ (y-1) .* ds)
@deriv_rule .^(x::AbstractArray, y::Real )           x     y * x .^ (y-1) .* ds
@deriv_rule .^(x::AbstractArray, y::AbstractArray)   x     y .* x .^ (y-1) .* ds

@deriv_rule .^(x::Real         , y::Real )           y     log(x) * x ^ y * ds
@deriv_rule .^(x::AbstractArray, y::Real )           y     sum( log(x) .* x .^ y .* ds)
@deriv_rule .^(x::Real         , y::AbstractArray)   y     log(x) .* x .^ y .* ds
@deriv_rule .^(x::AbstractArray, y::AbstractArray)   y     log(x) .* x .^ y .* ds

# division
@deriv_rule /(x::Real          , y::Real )           x     ds / y
@deriv_rule /(x::AbstractArray , y::Real )           x     ds ./ y

@deriv_rule /(x::Real          , y::Real )           y     -x * ds / (y * y)
@deriv_rule /(x::AbstractArray , y::Real )           y     sum(-x .* ds) / (y * y)

# dot division
@deriv_rule ./(x::Real         , y::Real )           x     ds / y
@deriv_rule ./(x::Real         , y::AbstractArray)   x     sum(ds ./ y)
@deriv_rule ./(x::AbstractArray, y::Real )           x     ds ./ y
@deriv_rule ./(x::AbstractArray, y::AbstractArray)   x     ds ./ y

@deriv_rule ./(x::Real         , y::Real )           y     -x * ds / (y * y)
@deriv_rule ./(x::Real         , y::AbstractArray)   y     -x * ds ./ (y .* y)
@deriv_rule ./(x::AbstractArray, y::Real )           y     -sum(x .* ds) / (y * y)
@deriv_rule ./(x::AbstractArray, y::AbstractArray)   y     -x .* ds ./ (y .* y)

# mod2pi
ReverseDiffSource.@deriv_rule mod2pi(x::Real)                          x     ds

# transpose
@deriv_rule transpose(x::Real )                      x     ds
@deriv_rule transpose(x::AbstractArray)              x     transpose(ds)

# erf, erfc, gamma, beta, lbeta, lgamma
@deriv_rule erf(x::Real)                       x     2/sqrt(π) * exp(-x  * x)  * ds
@deriv_rule erf(x::AbstractArray)              x     2/sqrt(π) * exp(-x .* x) .* ds

@deriv_rule erfc(x::Real)                      x     -2/sqrt(π) * exp(-x  * x)  * ds
@deriv_rule erfc(x::AbstractArray)             x     -2/sqrt(π) * exp(-x .* x) .* ds

@deriv_rule gamma(x::Real)                     x     polygamma(0,x)  * gamma(x)  * ds
@deriv_rule gamma(x::AbstractArray)            x     polygamma(0,x) .* gamma(x) .* ds

@deriv_rule lgamma(x::Real)                    x     polygamma(0,x)  * ds
@deriv_rule lgamma(x::AbstractArray)           x     polygamma(0,x) .* ds

@deriv_rule beta(x::Real         , y::Real)            x   beta(x,y)  * (digamma(x)-digamma(x+y))  * ds
@deriv_rule beta(x::AbstractArray, y::AbstractArray)   x   beta(x,y) .* (digamma(x)-digamma(x+y)) .* ds
@deriv_rule beta(x::Real         , y::Real)            y   beta(x,y)  * (digamma(y)-digamma(x+y))  * ds
@deriv_rule beta(x::AbstractArray, y::AbstractArray)   y   beta(x,y) .* (digamma(y)-digamma(x+y)) .* ds

@deriv_rule lbeta(x::Real         , y::Real)            x   (polygamma(0,x)-polygamma(0,x+y))  * ds
@deriv_rule lbeta(x::AbstractArray, y::AbstractArray)   x   (polygamma(0,x)-polygamma(0,x+y)) .* ds
@deriv_rule lbeta(x::Real         , y::Real)            y   (polygamma(0,y)-polygamma(0,x+y))  * ds
@deriv_rule lbeta(x::AbstractArray, y::AbstractArray)   y   (polygamma(0,y)-polygamma(0,x+y)) .* ds
