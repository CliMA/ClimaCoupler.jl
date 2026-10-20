#=
Control-volume weights for the exchange-grid scatter. Each SE node owns a box in
reference coordinates; the boxes of an element tile it, and each box's area is the node's
quadrature weight `(Jw)_n`. The weight of polygon `Ω_k` at node `n` is the area of `Ω_k`
inside the box, `B̃_{kn} = ∫_{Ω_k ∩ box(n)} dA ≥ 0`.
=#

"""
    ControlVolumes(space)

Mass-matched control volumes of every element of the spectral-element `space`: node `(i, j)`
of element `e` owns `[ξ[i, j, e], ξ[i + 1, j, e]] × [η[j, e], η[j + 1, e]]` in reference
coordinates. The row edges `η` give each row of boxes its row's total weight, and the column
edges `ξ` within each row give each box its weight `(Jw)_{ij}`; each edge is found by
bisection on the exact area under the interpolated Jacobian.
"""
struct ControlVolumes
    η::Matrix{Float64}     # (Nq + 1) × n_elem
    ξ::Array{Float64, 3}   # (Nq + 1) × Nq × n_elem
end

function ControlVolumes(space)
    ξs, ws = CC.Quadratures.quadrature_points(Float64, CC.Spaces.quadrature_style(space))
    Nq = length(ξs)
    Jw = get_ConservativeRegriddingCCExt().se_node_weights(space)
    n_elem = length(Jw) ÷ Nq^2
    # Φ(x)[p] = ∫_{-1}^x ϕ_p, exact by GLL quadrature on [-1, x]
    Φ(x) =
        (x + 1) / 2 *
        (CC.Quadratures.interpolation_matrix(@.((x + 1) * (ξs + 1) / 2 - 1), ξs)' * ws)
    η = zeros(Nq + 1, n_elem)
    ξ = zeros(Nq + 1, Nq, n_elem)
    for e in 1:n_elem
        W = reshape(Jw[((e - 1) * Nq ^ 2 + 1):(e * Nq ^ 2)], Nq, Nq)
        J = W ./ (ws * ws') # nodal Jacobian
        row_weight = J' * ws
        η[1, e], η[Nq + 1, e] = -1, 1
        for j in 1:(Nq - 1)
            η[j + 1, e] =
                _increasing_root(y -> sum(row_weight .* Φ(y)), sum(W[:, 1:j]), η[j, e], 1.0)
        end
        for j in 1:Nq
            JΔΦ = J * (Φ(η[j + 1, e]) - Φ(η[j, e]))
            ξ[1, j, e], ξ[Nq + 1, j, e] = -1, 1
            for i in 1:(Nq - 1)
                ξ[i + 1, j, e] =
                    _increasing_root(x -> sum(Φ(x) .* JΔΦ), sum(W[1:i, j]), ξ[i, j, e], 1.0)
            end
        end
    end
    return ControlVolumes(η, ξ)
end

# Root of an increasing function `f(x) = target` on [lo, hi], by bisection
function _increasing_root(f, target, lo, hi)
    for _ in 1:200
        mid = (lo + hi) / 2
        (mid == lo || mid == hi) && break
        f(mid) < target ? (lo = mid) : (hi = mid)
    end
    return (lo + hi) / 2
end

"""
    control_volume_weights(cv::ControlVolumes, space, elem, polygon) -> Matrix{Float64}

Control-volume weights of one exchange-grid polygon inside SE element `elem`: an
`Nq × Nq` block with `B̃[i, j] = ∫_{polygon ∩ box(i, j)} J dξ dη ≥ 0`, where the
polygon's vertices are mapped to reference coordinates and joined by straight edges
there. The outermost box edges lie at `±2` to keep vertices that invert to just
outside `[-1, 1]`.
"""
function control_volume_weights(cv::ControlVolumes, space, elem, polygon)
    CRExt = get_ConservativeRegriddingCCExt()
    GI, GO = CRExt.GI, CRExt.GO
    Nq = size(cv.ξ, 2)
    bary, wq = CR.TriangleQuadrature.reference_rule(2 * (Nq - 1))
    vertices = [
        CRExt.inverse_element_map(space, elem, (GI.x(v), GI.y(v), GI.z(v))) for
        v in GI.getpoint(GI.getexterior(polygon))
    ]
    reference_polygon = GI.Polygon([GI.LinearRing(vertices)])
    edge(edges, i) = i == 1 ? -2.0 : i == Nq + 1 ? 2.0 : edges[i]
    B = zeros(Nq, Nq)
    for j in 1:Nq, i in 1:Nq
        ξ0, ξ1 = edge(view(cv.ξ, :, j, elem), i), edge(view(cv.ξ, :, j, elem), i + 1)
        η0, η1 = edge(view(cv.η, :, elem), j), edge(view(cv.η, :, elem), j + 1)
        box =
            GI.Polygon([GI.LinearRing([(ξ0, η0), (ξ1, η0), (ξ1, η1), (ξ0, η1), (ξ0, η0)])])
        piece = GO.intersection(GO.ConvexConvexSutherlandHodgman(), reference_polygon, box)
        piece_vertices = collect(GI.getpoint(GI.getexterior(piece)))[1:(end - 1)]
        for t in CR.TriangleQuadrature.fan_triangulate(piece_vertices)
            A = CR.TriangleQuadrature.planar_triangle_area(t)
            for (λ, w) in zip(bary, wq)
                ξ = λ[1] * t[1][1] + λ[2] * t[2][1] + λ[3] * t[3][1]
                η = λ[1] * t[1][2] + λ[2] * t[2][2] + λ[3] * t[3][2]
                B[i, j] += 2A * w * CRExt.element_jacobian_at(space, elem, ξ, η)
            end
        end
    end
    return B
end
