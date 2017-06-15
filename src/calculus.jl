"""
    trapezint(x::AbstractVector, y::AbstractVector)
    trapezint(step::Real, y::AbstractVector)
Calculate the integral of `y` versus `x` using the trapezoidal rule.
The second method assumes equally-spaced `y`.
"""
function trapezint{X<:Real, Y<:Real}(x::AbstractVector{X}, y::AbstractVector{Y})
  integral::Float64 = 0.0
  @inbounds for i in 1:(length(x)-1)
    integral += (x[i+1] - x[i])*(y[i+1] + y[i])
  end
  integral *= 0.5
end
function trapezint{Y<:Real}(step::Real, y::AbstractVector{Y})
  integral::Float64 = 0.0
  @inbounds for i in 1:(length(y)-1)
    integral += step*(y[i+1] + y[i])
  end
  integral *= 0.5
end

"""
    fdderiv(x::AbstractVector, y::AbstractVector)
    fdderiv(step::Real, y::AbstractVector)
Calculate the derivative of `y` versus `x` using the forward difference approximation.
The second method assumes equally-spaced `y`. Notice that the last element of the
derivative cannot be computed and therefore it is equalled to the second-to-last.
"""
function fdderiv(x::AbstractVector, y::AbstractVector)
  der = zeros(y)
  @inbounds @simd for i in 1:length(y)-1
    der[i] = (y[i+1]-y[i])/(x[i+1] - x[i])
  end
  der[end] = der[end-1]
  return der
end
function fdderiv(x::Real, y::AbstractVector)
  der = zeros(y)
  @inbounds @simd for i in 1:length(y)-1
    der[i] = (y[i+1]-y[i])/(x)
  end
  der[end] = der[end-1]
  return der
end

"""
    fdderiv!(der::AbstractVector, x::AbstractVector, y::AbstractVector)
    fdderiv!(der::AbstractVector, step::Real, y::AbstractVector)
Calculate the derivative of `y` versus `x` using the forward difference approximation
and write it in-place for `der`.
The second method assumes equally-spaced `y`. Notice that the last element of
the derivative
cannot be computed and therefore it is equalled to the second-to-last.
"""
function fdderiv!(der::AbstractVector, x::AbstractVector, y::AbstractVector)
  @inbounds @simd for i in 1:length(y)-1
    der[i] = (y[i+1]-y[i])/(x[i+1] - x[i])
  end
  der[end] = der[end-1]
  return der
end
function fdderiv!(der::AbstractVector, x::Real, y::AbstractVector)
  @inbounds @simd for i in 1:length(y)-1
    der[i] = (y[i+1]-y[i])/(x)
  end
  der[end] = der[end-1]
  return der
end

"""Return the order of magnitude of a number"""
mag(x::Real) = ceil(log10(x))
