using SpecialFunctions
using HypergeometricFunctions
using Memoize
using FastGaussQuadrature
using LinearAlgebra
using TaylorSeries
using QuadGK

# To adapt to use a different (non-isochrone) zeroth-order, the
# following quantities must be changed: phi0l, rho0l, Omegal, mu0l
# Additionally, the variable substitution in Ikl (which is x = pi*s for the
# isochrone model) should be adjusted for optimal performance.

# Zeroth-order potential and density for isochrone model (see Eq. (5.1))
function phi00(r)
    return 1.0 ./ (1.0 .+ sqrt.(1.0 .+ r.^2))
end

function phi0l(l,r)
    return r.^l .* phi00(r).^(1 + 2l)
end

function rho0l(l,r)
    return -(((1 + 2*l)*r.^l.*(1 .+ sqrt.(1 .+ r.^2)).^(-3 - 2*l).*(3 .+ 2*r.^2 .+ 3*sqrt.(1 .+ r.^2) .+ 2*l*(1 .+ r.^2 .+ sqrt.(1 .+ r.^2))))./(1 .+ r.^2).^(3/2))
end

# Mellin-space weight function corresponding to isochrone potential/density
# This comes from ωl[l,s] in the Mathematica script isochrone_expressions.m
function Omegal(l::Int64,s::Float64)::Float64
    return abs( (gamma(l + 3/2 + im*s) * gamma(l/2 + 1/4 + im*s/2))/gamma(3*l/2 + 7/4 + im*s/2))^2
end

# Un-reduced incomplete beta function
function incbeta(x,a,b)
    if (b>0)
        return beta(a,b)*beta_inc(a,b,x)[1]
    end
    return x^a/a * _₂F₁(a,1-b,a+1,x)
end

# The recurrence coeff for the orthogonal polynomials/basis function (Eq. (5.2))
@memoize function betajl(j::Int64,l::Int64)::Float64
    if (j == -1)
        return 0
    end

    if (j == 0)
        # This could be found numerically, but for the isochrone we know the exact value (see Eq. (G.4))
        mu0l = (2*l+1)/(24*pi) * exp(lgamma(2*l+2) + lgamma(l+0.5) - lgamma(3*l+1.5)) * (4.0^(-l) - 2*(2*l+1)*incbeta(0.5,3*l+2.5,-l-0.5))
        return mu0l
    end

    return Ikl(j,l)/Ikl(j-1,l)
end

# The orthogonal polynomials (in Mellin-space)
### @memoize function pnl(n::Int64,l::Int64,s::Float64)::Float64
###     if (n == -1)
###         return 0
###     end
###     if (n == 0)
###         return 1
###     end
###     return s*pnl(n-1,l,s) - betajl(n-1,l)*pnl(n-2,l,s)
### end

function pnl(n::Int64,l::Int64,s::Float64)::Float64
    ret = Array{Float64}(undef, n+1)

    ret[1] = 1.0

    if (n==0)
        return ret[1]
    end
    
    ret[2] = s*ret[1]

    if (n==1)
        return ret[2]
    end

    for j in (3:(n+1))
        ret[j] = s*ret[j-1] - betajl(j-2,l)*ret[j-2]
    end

    return ret[n+1]
end


# 'Normalisation' constants I_{kl} (Eq. (5.3) & (5.6))
# 
# Note the variable substitution s = x/pi, which comes about because
# we determined that Omegal(l,s) asymptotically decays as
# exp(-pi*s). For a different weight function (& for optimal
# performance) this variable substitution should be tuned to match the
# correct exponential decay.
@memoize function Ikl(k::Int64,l::Int64)::Float64
    f(x) = exp(x + log(Omegal(l,x/pi)) + 2*log(abs(pnl(k,l,x/pi))))
    x, w = gausslaguerre(max(k+l,20))
    return exp(log(2/pi) + log(dot(w,f.(x))) + 2*log(2l+1) - (2l+6)*log(2) - 2*log(pi))
end

# A vector of Ikl(k,l) values
@memoize function Ikl_vec(k,l)
    return map(n -> Ikl(n,l), 0:k)
end

# Return a vector of the (0...n)th derivatives of f
function nthDerivs(f,x,n)
    t = Taylor1(n)
    facts = map(i -> gamma(i+1), 0:n)
    if (length(x) == 1)
        series = f(x + t).coeffs
        return (series .* facts)'
    end
    
    series = map(s -> s.coeffs',f(x .+ t))
    return reduce(vcat, series) .* facts'
end

# Repeated differentiation of zeroth-order potential/density functions
function phi0lderivs(p,l,r)
    f = s -> phi0l(l,exp.(s))
    return nthDerivs(f,log.(r),p)
end

function rho0lderivs(p,l,r)
    f = s -> rho0l(l,exp.(s))
    return nthDerivs(f,log.(r),p)
end

# matrix coeffs of transformation from phi0lderivs to phinl (Eq. (5.14))
@memoize function Anjl(n,j,l)
    if (n<0 || j<0 || n<j)
        return 0
    end

    if (n==0 && j==0)
        return 1
    end

    if (j==0)
        return 0.5*Anjl(n-1,0,l) + betajl(n-1,l)*Anjl(n-2,0,l)
    end

    return Anjl(n-1,j-1,l) + 0.5*Anjl(n-1,j,l) + betajl(n-1,l)*Anjl(n-2,j,l)
end

# matrix coeffs of transformation from rho0lderivs to rhonl
@memoize function Bnjl(n,j,l)
    if (n<0 || j<0 || n<j)
        return 0
    end

    if (n==0 && j==0)
        return 1
    end

    if (j==0)
        return 2.5*Bnjl(n-1,0,l) + betajl(n-1,l)*Bnjl(n-2,0,l)
    end

    return Bnjl(n-1,j-1,l) + 2.5*Bnjl(n-1,j,l) + betajl(n-1,l)*Bnjl(n-2,j,l)
end

# matrices constructed from the matrix coeffs above
@memoize function Amat(nmax,l)
    m = zeros(nmax+1,nmax+1)

    for n in (0:nmax), j in (0:nmax)
        m[n+1,j+1] = (-1)^n * Anjl(n,j,l)
    end

    return m
end

@memoize function Bmat(nmax,l)
    m = zeros(nmax+1,nmax+1)

    for n in (0:nmax), j in (0:nmax)
        m[n+1,j+1] = (-1)^n * Bnjl(n,j,l)
    end

    return m
end

# The actual basis functions
function phinl(nmax::Int64,l::Int64,r::Float64)
    x = r/Rbasis # Dimensionless radius
    return -1.0 * (Amat(nmax,l) * phi0lderivs(nmax,l,x)') ./ sqrt.(Ikl_vec(nmax,l)) * (sqrt(G/Rbasis))
end

function rhonl(nmax::Int64,l::Int64,r::Float64)
    x = r/Rbasis # Dimensionless radius
    return -1.0 * (Bmat(nmax,l) * rho0lderivs(nmax,l,x)') ./ sqrt.(Ikl_vec(nmax,l)) / ((sqrt(G)*Rbasis^(5/2)*(4*pi)))
end

# Alternative implementations of phinl and rhonl that appear to be slower
function phinl2(n::Int64,l::Int64,r::Float64)
    t = Taylor1(n)
    x = log(r/Rbasis)

    # the change of variables x = log(r/R) means that r d/dr = d/dx,
    # which reduces the total number of required algebraic operations
    # on the Taylor series

    ret = Array{Taylor1{Float64}}(undef, n+1)

    ret[1] = phi0l(l,exp(x+t))

    if (n==0)
        return ret
    end
    
    ret[2] = derivative(ret[1]) + 0.5*ret[1]

    if (n==1)
        return ret
    end

    for j in (3:n+1)
        ret[j] = derivative(ret[j-1]) + 0.5*ret[j-1] + betajl(j-2,l)*ret[j-2]
    end
    
    sign_correction = map(j -> (-1)^j, (1:(n+1)))

    return map(y -> y.coeffs[1], ret) .* sign_correction ./ sqrt.(Ikl_vec(n,l)) * (sqrt(G/Rbasis))
end

function rhonl2(n::Int64,l::Int64,r::Float64)
    t = Taylor1(n)
    x = log(r/Rbasis)

    ret = Array{Taylor1{Float64}}(undef, n+1)

    ret[1] = rho0l(l,exp(x+t))

    if (n==0)
        return ret
    end

    ret[2] = derivative(ret[1]) + 2.5*ret[1]

    if (n==1)
        return ret
    end

    for j in (3:n+1)
        ret[j] = derivative(ret[j-1]) + 2.5*ret[j-1] + betajl(j-2,l)*ret[j-2]
    end
    
    sign_correction = map(j -> (-1)^j, (1:(n+1)))

    return map(y -> y.coeffs[1], ret) .* sign_correction ./ sqrt.(Ikl_vec(n,l)) / ((sqrt(G)*Rbasis^(5/2)*(4*pi)))
end



# Check that the potential & density basis functions are mutually
# orthogonal, for n=0...nmax (at fixed l).
# Returns a matrix that should have 1.0 along the diagonal (and ~0 for
# all other entries) if the functions are correctly orthogonal.
function checkOrthogIso(l,nmax)
    f = r -> phinl(nmax,l,r)
    g = r -> rhonl(nmax,l,r)
    integrand = r -> (r == Inf) ? 0 : r^2*f(r)*g(r)'
    normalisation = -1.0
    result = quadgk(integrand, 0, Inf)
    return result[1] ./ normalisation
end



########################################
# The following functions are for interfacing with JBF's code
# Note that for consistency with JBF's code,
# (1) "np" is off-by-1 from the "n" used in the functions above
# (2) the order of n & l is swapped

function UlnpIso(l::Int64,np::Int64,r::Float64)
    return phinl(np-1,l,r)[np]
end

function DlnpIso(l::Int64,np::Int64,r::Float64)
    return rhonl(np-1,l,r)[np]
end

function tabUlnpIso!(l::Int64,r::Float64,tabUlnp::Array{Float64,1})
    result = phinl(nradial-1,l,r)
    for n=1:nradial
        tabUlnp[n] = result[n]
    end
end


function profile_test(f,x)
    for r = 0.0:0.1:x
        f(r);
    end
end

