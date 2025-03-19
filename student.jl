"""
    student(x, y, f)

1-dimensional periodogram of unevenly sampled data using the Student-T distribution

x = independant variable
y = dependant variable
f = frequency vector
"""
function student(x, y, f)
    M = length(x)
    scv = sincos.(2*pi.*f*x')
    C = ((last.(scv)*y).^2 .+ (first.(scv)*y).*2)/M
    (1 .- 2 .*C./(y'*y)).^((2-M)/2)
end
