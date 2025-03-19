
function eigen2(m::Matrix)
    a, b, c, d = reshape(m, :)
    Δ = (a - d)/2
    T, S = (a + d)/2, sqrt(Δ^2 + b*c)
    θ = acos(Δ/S)/2
    if θ >= π/4 θ += π end
    y, x = sincos(θ)
    (values=[T - S, T + S], vectors=[y -x; -x -y])
end
