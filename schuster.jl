"""
    schuster(x, y, s, f)

1-dimesional periodogram of unevenly sampled data using basic periodogram

x = independant variable
y = dependant variable
s = RMS error
f = frequency vector
"""
function schuster(x::AbstractArray, y::AbstractArray, s::AbstractArray, f::AbstractArray)
    scv  = sincos.(2pi.*f*x')
    sqrt.((last.(scv)*(y./s)).^2 .+ (first.(scv)*(y./s)).^2)
end
