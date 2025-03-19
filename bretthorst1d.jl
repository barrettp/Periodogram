"""
    bretthorst1d(x, y, s, f)

1-dimensional periodogram of unevenly sampled data using Bretthorst's Bayesian algorithm

x = independant variable
y = dependant variable
s = RMS error
f = frequency vector
"""
function brett1d(x::AbstractArray, y::AbstractArray, s::AbstractArray, freq::AbstractArray; zloge=0.0)
    M, N = length(x), length(freq)
    B, h2bar, st = zeros(eltype(y), (N, 2)), zeros(eltype(y), N), zeros(eltype(y), N)
    stloge, sig, phat = zeros(eltype(y), N), zeros(eltype(y), N), zeros(eltype(y), N)

    Threads.@threads for j=1:N
        Gj = hcat(cos.(2pi*freq[j].*x), sin.(2pi*freq[j].*x))
        lj, ej = eigen(Hermitian(Gj'*Gj))
        # lj, ej = eigen2(Gj'*Gj)
        hj = ej*Gj'*y./sqrt.(lj)
        B[j,:] = ej*hj./sqrt.(lj)
        h2bar[j] = hj'*hj/2
        stl = log(1 - (hj'*hj)./(y'*y)) .* ((2 - M)/2)
        stloge[j] = stl
        st[j] = abs(zloge) != 0 ? exp(stl - zloge) : 0.0
        sig[j] = sqrt((y'*y - hj'*hj)/(M-4))
        phat[j] = hj'*hj * st[j]
    end
    (B, h2bar, st, stloge, sig, phat)

end

function eigen2(m::Matrix)
    a, b, c, d = reshape(m, :)
    Δ = (a - d)/2
    T, S = (a + d)/2, sqrt(Δ^2 + b*c)
    θ = acos(Δ/S)/2
    if θ >= π/4 θ += π end
    y, x = sincos(θ)
    (values=[T - S, T + S], vectors=[y -x; -x -y])
end