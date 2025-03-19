"""
    bretthorst2d(x, y, z, f; zloge=0)

2-dimensional periodogram of unevenly sampled data using Bretthorst's Bayesian algorithm

x = independant variable
y = dependant variable
s = RMS error
f = frequency vector
"""
function brett2d(x::AbstractArray, y::AbstractArray, s::AbstractArray, freq::AbstractArray; zloge=0.0)
    M, N = length(x), length(freq)
    B, h2bar, st = zeros(eltype(y), (N, N, 4)), zeros(eltype(y), (N, N)), zeros(eltype(y), (N, N))
    stloge, sig, phat = zeros(eltype(y), (N, N)), zeros(eltype(y), (N, N)), zeros(eltype(y), (N, N))

    for j=1:N
        f1 = freq[j]
        Threads.@threads for k=1:N
            f2 = (freq[k] == f1 ? 1.01 : 1.0)*freq[k]
            Gjk = hcat(cos.(2pi*f1.*x), sin.(2pi*f1.*x), cos.(2pi*f2.*x), sin.(2pi*f2.*x))
            F = eigen(Hermitian(Gjk'*Gjk))
            ljk, ejk = F.values, F.vectors
            hjk = ejk*Gjk'*y./sqrt.(ljk)
            B[j,k,:] = ejk*hjk./sqrt.(ljk)
            h2bar[j,k] = hjk'*hjk/2
            # stl = log(1 - (hj'*hj)./(y'*y)) .* ((2 - M)/2)
            stloge[j,k] = 0.0
            st[j,k] = abs(zloge) != 0 ? exp(stl - zloge) : 0.0
            # sig[j,k] = sqrt((y'*y - hjk'*hjk)/(M-4))
            # phat[j,k] = hjk'*hjk * st[j,k]
        end
    end
    (B, h2bar, st, stloge, sig, phat)
end
