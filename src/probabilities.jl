export renewal_counting
"""
```julia
renewal_counting(n::Int, t, pdf, first_pdf = pdf)
```
Recursively calculates the p.d.f. of the `n`-th event happening at time `t`, as
well as the probability to have `n` events occuring up to time `t`.
# Inputs:
1. `n::Int` : The number of events you care about. Notice that everything
  starting from `k=1` up to `k=n` will be calculated (and returned).
2. `t::Vector` : The time vector versus which the `pdf` is calculated. Has to have a
  `step()` method defined (i.e. be equally-spaced).
3. `pdf::Vector` : The vector containing the histogram (prob. density function)
  of the
  interarrival times of the events. It should be integral-normalized versus `t`,
  and if it is not, this is taken care internally. It is expected that the
  histogram/p.d.f. is calculated with **left-closed intervals**! It is also
  expected that the p.d.f. has finite support, i.e. it is 0 for any time ∉ `t`.
4. `first_pdf::Vector` : [optional] The p.d.f. of the interarrival times of the
  *first* event, in
  case it is different from the general p.d.f.
# Outputs:
1. `realt::Vector` : Time vector of length: `L = l*(n+1) - n` where `l = length(t)`
2. `Qn::Matrix` : A `L×n` matrix. Each column `j` is the **p.d.f.** of the
  `j`-th event versus `realt`,  for `j in 1:n`
3. `Pn::Matrix` : A `L×(n+1)` matrix. Each column `j` is the **probability** to
  have a total of `j` events up to time `realt` for `j in 0:n`.

Notice that the output `Qn` represents probability density function, while the
output `Pn` represents commulative probability function (*direct* probability).
"""
function renewal_counting(n::Int, t::AbstractVector, pdf::AbstractVector,
  first_pdf = pdf)

  if length()
  l = length(t)
  tmax = maximum(t)
  integral = trapezint(t, pdf)

  if !isapprox(integral, 1)
    warn("Integral of p.d.f. versus t was not normalized,
using a re-normalized copy instead...")
    f = deepcopy(pdf)
    f ./= integral
  else
    f = pdf
  try
    stepp = step(t)
  catch
    throw(ArgumentError("Input `t` must be equally spaced (have `step()`)"))
  end


  #this is the time vector my functions will be functions of
  realt = 0:stepp:ceil((n+1)*tmax, -Int(mag(stepp)))
  L = length(realt) #should be l*n+1 - n because you only count 0.0 once
  append!(f, zeros(L - l)) # extend f for realt
  # Calculate average event time:
  κ = trapezint(stepp, realt.*f)
  W = zeros(realt)
  # Calculate general survival function W:
  for i in 1:l
    W[i] = 1 - trapezint(t[1:i], f[1:i])
  end

  # Define the matrix P: each column a P_n, from n=0
  # P_n[i] is probability to have n events up to time t[i]
  P = zeros(eltype(t), L, n+1)
  # P_0, probability to have 0 events up to time t (identical to
  # survival function of first event)
  W = zeros(realt)
  # Calculate general survival function W:
  for i in 1:l
    P[i, 1] = 1 - trapezint(t[1:i], first_pdf[1:i])
  end
  # Define the matrix Q: each column a Q_n from n=1
  # Qn is the p.d.f. of the n-th event happening at time t
  Q = zeros(eltype(t), L, n)
  # Calculate Q1: probability that 1st event happens at time t
  # equivalent to p.d.f. of the first event or to general p.d.f.
  # if they are the same
  Q1[1:l] .= first_pdf
  # Calculate Qn, Pn
  for j in 1:n
    # Qn:
    if j != 1
      reverse_Qnm1 = reverse(Q[:,j-1])
      for t_idx in 2:L #no reason to start from 1, gives 0 anyway
        Q[t_idx, j] =
        trapezint(realt[1:t_idx], f[1:t_idx].*reverse_Qnm1[(L-t_idx+1):end])
      end
      qint = trapezint(realt, Q[:, j])
      if !isapprox(qint, 1; rtol = 1e-14)
        Q[:, j] ./= qint
      end
    end

    # Pn:
    reverse_Qn = reverse(Q[:,j])
    for t_idx in 2:L #t_idx is the timepoint t I want to find the Pn(t)
      P[t_idx, j+1] =
      trapezint(realt[1:t_idx], W[1:t_idx].*reverse_Qn[(L-t_idx+1):end])
    end
  end
  return realt, Q, P
end
