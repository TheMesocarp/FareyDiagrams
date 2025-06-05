module FareyDiagrams

using Luxor
using Colors
using SimpleWeightedGraphs
using Graphs
using GLMakie
using GeometryTypes

# A module providing functionality to create and manipulate Farey diagrams and topographs.

# ---------------------------------
# Graph Construction Utilities
# ---------------------------------

"""Create an empty weighted graph and associated vertex-value dictionary."""
function create_custom_graph()
    g = SimpleWeightedGraph(0)
    vertex_value = Dict{Int, Set{Rational{Int}}}()
    return g, vertex_value
end

"""Create an empty simple graph (no weights) and associated vertex-point dictionary."""
function create_custom_Tgraph()
    g = SimpleGraph(0)
    vertex_value = Dict{Int, Set{Luxor.Point}}()
    return g, vertex_value
end

"""Add a vertex to a simple graph, storing a Luxor.Point in the vertex dictionary."""
function add_topovertex!(g::SimpleGraph, vertex_value::Dict{Int, Set{Luxor.Point}}, point::Luxor.Point)
    v_index = nv(g) + 1
    add_vertices!(g, 1)
    vertex_value[v_index] = Set([point])
    return v_index
end

"""Add a vertex to a weighted graph, storing a rational set in the vertex dictionary."""
function add_farey_vertex!(g::SimpleWeightedGraph, vertex_value::Dict{Int, Set{Rational{Int}}}, values::Set{Rational{Int}})
    v_index = nv(g) + 1
    add_vertices!(g, 1)
    vertex_value[v_index] = values
    return v_index
end

"""Add an edge between two vertices in a weighted Farey graph."""
function add_farey_edge!(g::SimpleWeightedGraph, v1::Int, v2::Int, weight::Float64)
    add_edge!(g, v1, v2, weight)
end

"""Update the weight of an existing edge in a weighted graph."""
function update_edge_weight!(g::SimpleWeightedGraph, v1::Int, v2::Int, new_weight::Float64)
    if has_edge(g, v1, v2)
        g.weights[v1, v2] = new_weight
        g.weights[v2, v1] = new_weight  # Ensure symmetry
    end
end

# ---------------------------------
# Farey Graph Construction
# ---------------------------------

"""
Build the Farey graph of order n. Returns a tuple (Graph, vertex_values, triangles).
- Graph: SimpleWeightedGraph representing the Farey diagram
- vertex_values: Dict mapping vertex index to Set{Rational{Int}} representing the fraction(s)
- triangles: Vector of integer triples representing new triangles introduced at this level
"""
function FareyGraph(n::Int)
    if n == 1
        F_1, Vval = create_custom_graph()
        # Base vertices: 1/1, infinity (±1/0), -1/1, 0/1
        v1 = add_farey_vertex!(F_1, Vval, Set(Rational(1, 1)))
        v2 = add_farey_vertex!(F_1, Vval, Set([Rational(1, 0), Rational(-1, 0)]))  # Infinity
        v3 = add_farey_vertex!(F_1, Vval, Set(Rational(-1, 1)))
        v4 = add_farey_vertex!(F_1, Vval, Set(Rational(0, 1)))

        # Connect base edges with initial weight 0.1
        add_farey_edge!(F_1, v1, v2, 0.1)
        add_farey_edge!(F_1, v2, v3, 0.1)
        add_farey_edge!(F_1, v3, v4, 0.1)
        add_farey_edge!(F_1, v1, v4, 0.1)

        triangles = [[v1, v2, v4], [v2, v3, v4]]
        return F_1, Vval, triangles
    else
        # Recursively build (n-1)-graph
        Graph_prev, Value_prev, triangles_prev = FareyGraph(n - 1)
        # Copy graph to modify
        CopyG = SimpleWeightedGraph(Graph_prev)
        for e in edges(CopyG)
            p1, p2 = src(e), dst(e)
            w = weight(e)
            f_1 = first(Value_prev[p1])
            f_2 = first(Value_prev[p2])
            # Handle sign adjustments for infinity
            if f_1 in [Rational(1, 0), Rational(-1, 0)] || f_2 in [Rational(1, 0), Rational(-1, 0)]
                if f_1 in [Rational(1, 0), Rational(-1, 0)]
                    if f_1 == Rational(1, 0)
                        f_2 < 0 ? f_1 = -f_1 : f_1 = f_1
                    else
                        f_2 > 0 ? f_1 = -f_1 : f_1 = f_1
                    end
                else
                    if f_2 == Rational(1, 0)
                        f_1 < 0 ? f_2 = -f_2 : f_2 = f_2
                    else
                        f_1 > 0 ? f_2 = -f_2 : f_2 = f_2
                    end
                end
            end
            # Compute mediant
            new_value = Rational(numerator(f_1) + numerator(f_2), denominator(f_1) + denominator(f_2))
            if w ≈ 0.1
                v_child = add_farey_vertex!(Graph_prev, Value_prev, Set(new_value))
                add_farey_edge!(Graph_prev, v_child, p1, 0.1)
                add_farey_edge!(Graph_prev, v_child, p2, 0.1)
                update_edge_weight!(Graph_prev, p1, p2, 1.0)
                push!(triangles_prev, [v_child, p1, p2])
            end
        end
        return Graph_prev, Value_prev, triangles_prev
    end
end

"""
Build an extended Farey graph variant (FareyGraph2) that introduces extra vertices.
"""
function FareyGraph2(n::Int)
    if n == 1
        F_1, Vval = create_custom_graph()
        # Base vertices plus additional ones: ±2/1, ±1/2
        v1 = add_farey_vertex!(F_1, Vval, Set(Rational(1, 1)))
        v2 = add_farey_vertex!(F_1, Vval, Set([Rational(1, 0), Rational(-1, 0)]))
        v3 = add_farey_vertex!(F_1, Vval, Set(Rational(-1, 1)))
        v4 = add_farey_vertex!(F_1, Vval, Set(Rational(0, 1)))
        v5 = add_farey_vertex!(F_1, Vval, Set(Rational(2, 1)))
        v6 = add_farey_vertex!(F_1, Vval, Set(Rational(-2, 1)))
        v7 = add_farey_vertex!(F_1, Vval, Set(Rational(1, 2)))
        v8 = add_farey_vertex!(F_1, Vval, Set(Rational(-1, 2)))
        # Initial edges
        add_farey_edge!(F_1, v1, v3, 0.1)
        add_farey_edge!(F_1, v2, v7, 0.1)
        add_farey_edge!(F_1, v4, v5, 0.1)
        add_farey_edge!(F_1, v4, v6, 0.1)
        add_farey_edge!(F_1, v2, v8, 0.1)
        return F_1, Vval
    else
        Graph_prev, Value_prev = FareyGraph2(n - 1)
        CopyG = SimpleWeightedGraph(Graph_prev)
        for e in edges(CopyG)
            p1, p2 = src(e), dst(e)
            w = weight(e)
            f_1 = first(Value_prev[p1])
            f_2 = first(Value_prev[p2])
            # Handle sign adjustments for infinity
            if f_1 in [Rational(1, 0), Rational(-1, 0)] || f_2 in [Rational(1, 0), Rational(-1, 0)]
                if f_1 in [Rational(1, 0), Rational(-1, 0)]
                    if f_1 == Rational(1, 0)
                        f_2 < 0 ? f_1 = -f_1 : f_1 = f_1
                    else
                        f_2 > 0 ? f_1 = -f_1 : f_1 = f_1
                    end
                else
                    if f_2 == Rational(1, 0)
                        f_1 < 0 ? f_2 = -f_2 : f_2 = f_2
                    else
                        f_1 > 0 ? f_2 = -f_2 : f_2 = f_2
                    end
                end
            end
            if w ≈ 0.1
                # Compute extended mediants
                new_value1 = Rational(numerator(f_1) + 2 * numerator(f_2), denominator(f_1) + 2 * denominator(f_2))
                new_value2 = Rational(2 * numerator(f_1) + numerator(f_2), 2 * denominator(f_1) + denominator(f_2))
                new_value3 = Rational(numerator(f_1) - 2 * numerator(f_2), denominator(f_1) - 2 * denominator(f_2))
                new_value4 = Rational(-2 * numerator(f_1) + numerator(f_2), -2 * denominator(f_1) + denominator(f_2))
                v_child1 = add_farey_vertex!(Graph_prev, Value_prev, Set(new_value1))
                v_child2 = add_farey_vertex!(Graph_prev, Value_prev, Set(new_value2))
                v_child3 = add_farey_vertex!(Graph_prev, Value_prev, Set(new_value3))
                v_child4 = add_farey_vertex!(Graph_prev, Value_prev, Set(new_value4))
                add_farey_edge!(Graph_prev, v_child2, p1, 0.1)
                add_farey_edge!(Graph_prev, v_child1, p2, 0.1)
                add_farey_edge!(Graph_prev, v_child4, p1, 0.1)
                add_farey_edge!(Graph_prev, v_child3, p2, 0.1)
                update_edge_weight!(Graph_prev, p1, p2, 1.0)
            end
        end
        return Graph_prev, Value_prev
    end
end

# ---------------------------------
# Triangle and Vertex Utilities
# ---------------------------------

"""Return the triples of new triangles (vertex indices) introduced at level n."""
function triangleF(n::Int)
    if n == 1
        _, _, base_triangles = FareyGraph(1)
        return base_triangles
    else
        _, _, M = FareyGraph(n)
        # For each triangle, list the vertex and its neighbors
        vertex_neighbors = [[v; neighbors(FareyGraph(n)[1], v)] for v in M]
        return vertex_neighbors
    end
end

"""Compute the continued fraction series of a rational number."""
function continued_fraction(r::Rational)
    series = Int[]
    num = numerator(r)
    den = denominator(r)
    while den != 0
        q = div(num, den)
        push!(series, q)
        num, den = den, num - q * den
    end
    return series
end

# ---------------------------------
# Angle Computations
# ---------------------------------

"""Compute the standard Farey diagram angle for a rational set {q}."""
function angle(a::Set{Rational{Int}})
    v = first(a)
    if v ∈ [Rational(-1, 0), Rational(1, 0)]
        return π
    elseif v >= -1 && v <= 1
        return -float(v) * π / 2
    else
        return (float(1 / v) + 2 * sign(v) * 1) * π / 2
    end
end

"""Compute the continued-fraction-based angle perturbation for equidistant layout."""
function theta(q::Rational)
    if q == 0//1
        return 0
    elseif q == 1//1
        return π / 2
    end
    cf = continued_fraction(q)
    angle_val = 0.0
    for (idx, a_i) in enumerate(cf[2:end], 2)
        partial_sum = sum(cf[2:idx])
        angle_val += (-1) ^ idx * ((1 / 2) ^ partial_sum) * π
    end
    return angle_val
end

"""Compute the equidistant Farey diagram angle for a rational set {q}."""
function equi_angle(a::Set{Rational{Int}})
    v = first(a)
    if v ∈ [Rational(-1, 0), Rational(1, 0)]
        return π
    elseif v >= 0 && v <= 1
        return mod(theta(v), 2π)
    elseif v >= -1 && v < 0
        return mod(-theta(-v), 2π)
    elseif v > 1
        return mod(π - theta(1 // v), 2π)
    else
        return mod(π + theta(-1 // v), 2π)
    end
end

# ---------------------------------
# Topograph Construction
# ---------------------------------

# Global radius for geometric computations
const my_radius = 600.0

"""Compute the centroid of a triangle (given by vertex indices) on the Farey diagram."""
function centertriangle(triangle::Vector{Int}, Farey::Tuple)
    DictG = Farey[2]
    # Compute coordinates for each vertex
    coords = [
        Luxor.Point(my_radius * cos(equi_angle(DictG[idx])), my_radius * sin(equi_angle(DictG[idx])))
        for idx in triangle
    ]
    center = Luxor.Point((coords[1].x + coords[2].x + coords[3].x) / 3,
                          (coords[1].y + coords[2].y + coords[3].y) / 3)
    return center
end

"""Find an existing vertex corresponding to a triangle centroid or add a new one."""
function get_or_add_vertex!(G::SimpleGraph, vertex_dict::Dict{Int, Set{Luxor.Point}}, triangle::Vector{Int}, Farey::Tuple)
    center = centertriangle(triangle, Farey)
    for (idx, points) in vertex_dict
        if any(≈(center), points)
            return idx
        end
    end
    return add_topovertex!(G, vertex_dict, center)
end

"""Compute centroids for every triangle up to level m in the Farey graph."""
function equitopograph(m::Int, G::Tuple)
    DictG = G[2]
    T = Luxor.Point[]
    for i in 1:m-1
        for t in triangleF(i)
            push!(T, centertriangle(t, G))
        end
    end
    return T
end

"""Recursively build the topograph up to level m."""
function graphtopo(m::Int)
    if m == 1
        topo_graph, vertex_dict = create_custom_Tgraph()
        Ffarey = FareyGraph(1)
        # Base centroids for triangles of level 1
        T = equitopograph(1, Ffarey)
        v1 = add_topovertex!(topo_graph, vertex_dict, T[1])
        v2 = add_topovertex!(topo_graph, vertex_dict, T[2])
        add_edge!(topo_graph, v1, v2)
        return topo_graph, vertex_dict
    else
        Farey = (FareyGraph(m)[1], FareyGraph(m)[2])
        prev_graph, prev_dict = graphtopo(m - 1)
        T_new = equitopograph(m, Farey)
        A = triangleF(m)
        B = triangleF(m - 1)
        # Connect new and existing centroids where triangles share an edge
        for t in B
            for p in A
                if length(intersect(t, p)) == 2
                    c_i = centertriangle(t, Farey)
                    c_j = centertriangle(p, Farey)
                    v_i = get_or_add_vertex!(prev_graph, prev_dict, t, Farey)
                    v_j = get_or_add_vertex!(prev_graph, prev_dict, p, Farey)
                    add_edge!(prev_graph, v_i, v_j)
                end
            end
        end
        return prev_graph, prev_dict
    end
end

"""Compute centroids using standard (non-equi) angles up to level m."""
function topograph(m::Int, G::Tuple)
    DictG = G[2]
    T = Luxor.Point[]
    for i in 1:m-1
        for t in triangleF(i)
            # Coordinates using standard angle
            coords = [
                Luxor.Point(my_radius * cos(angle(DictG[idx])), my_radius * sin(angle(DictG[idx])))
                for idx in t
            ]
            center = Luxor.Point((coords[1].x + coords[2].x + coords[3].x) / 3,
                                  (coords[1].y + coords[2].y + coords[3].y) / 3)
            push!(T, center)
        end
    end
    return T
end

# ---------------------------------
# Plotting Utilities
# ---------------------------------

# Helper: Compute line-arc intersection center for drawing arcs
function arc_center_radius(coord1::Tuple{Float64, Float64}, coord2::Tuple{Float64, Float64})
    x1, y1 = coord1
    x2, y2 = coord2
    if y1 == 0
        y1 = 1e-10
    end
    cx = (y1 * x2^2 - y2 * x1^2 + (y2 - y1) * y2 * y1) / (y1 * x2 - y2 * x1)
    cy = -(x1 / y1) * cx + ((x1^2) / y1) + y1
    r = sqrt((x1 - cx)^2 + (y1 - cy)^2)
    return (cx, cy), r
end

"""Draw an arc between two points with given center."""
function arc2r(ax, center::Tuple{Float64, Float64}, pt1::Tuple{Float64, Float64}, pt2::Tuple{Float64, Float64}; color=:black)
    x_c, y_c = center
    x1, y1 = pt1
    x2, y2 = pt2
    θ1 = atan(y1 - y_c, x1 - x_c)
    θ2 = atan(y2 - y_c, x2 - x_c)
    if θ2 < θ1
        θ2 += 2π
    end
    Makie.arc!(ax, Point2f(x_c, y_c), sqrt((x1 - x_c)^2 + (y1 - y_c)^2), θ1, θ2; color=color)
end

"""Plot the Farey graph with arcs (standard angle)."""
function plot_farey_graph(ax, m::Int, radius::Float64; variant::Int=1)
    if variant == 1
        Graph, DictG = (FareyGraph(m)[1], FareyGraph(m)[2])
        c = :blue
    else
        Graph, DictG = (FareyGraph2(m)[1], FareyGraph2(m)[2])
        c = :green
    end
    Edges = edges(Graph)
    # Draw main circle
    θ = range(0, 2π, length=200)
    circle_x = radius .* cos.(θ)
    circle_y = radius .* sin.(θ)
    lines!(ax, circle_x, circle_y, color=:black)
    # Draw vertices
    for (k, pts) in DictG
        angle_rad = angle(pts)
        coord_x = radius * cos(angle_rad)
        coord_y = radius * sin(angle_rad)
        scatter!(ax, [coord_x], [coord_y]; color=:red, markersize=8)
        # Label
        frac = first(pts)
        label = if denominator(frac) == 1
            string(numerator(frac))
        elseif frac ∈ [Rational(-1, 0), Rational(1, 0)]
            "∞"
        else
            string(numerator(frac), "/", denominator(frac))
        end
        text!(ax, label, position=(coord_x * 1.1, coord_y * 1.1), align=(:center, :center), fontsize=12)
    end
    # Draw edges as arcs
    for e in Edges
        p1, p2 = src(e), dst(e)
        angle1 = angle(DictG[p1])
        angle2 = angle(DictG[p2])
        coord1 = (radius * cos(angle1), radius * sin(angle1))
        coord2 = (radius * cos(angle2), radius * sin(angle2))
        center, _ = arc_center_radius(coord1, coord2)
        f_1 = first(DictG[p1]); f_2 = first(DictG[p2])
        if f_1 in [Rational(1, 0), Rational(-1, 0)] || f_2 in [Rational(1, 0), Rational(-1, 0)]
            if f_1 in [Rational(1, 0), Rational(-1, 0)]
                if f_2 < 0
                    arc2r(ax, center, coord2, coord1; color=c)
                else
                    arc2r(ax, center, coord1, coord2; color=c)
                end
            else
                if f_1 < 0
                    arc2r(ax, center, coord2, coord1; color=c)
                else
                    arc2r(ax, center, coord1, coord2; color=c)
                end
            end
        elseif f_1 < f_2
            arc2r(ax, center, coord1, coord2; color=c)
        else
            arc2r(ax, center, coord2, coord1; color=c)
        end
    end
    hideticks!(ax)
    return ax
end

"""Plot the Farey graph with equidistant layout."""
function plot_farey_graph_equi(ax, m::Int, radius::Float64; variant::Int=1)
    if variant == 1
        Graph, DictG = (FareyGraph(m)[1], FareyGraph(m)[2])
        c = :blue
    else
        Graph, DictG = (FareyGraph2(m)[1], FareyGraph2(m)[2])
        c = :green
    end
    Edges = edges(Graph)
    θ = range(0, 2π, length=200)
    circle_x = radius .* cos.(θ)
    circle_y = radius .* sin.(θ)
    lines!(ax, circle_x, circle_y, color=:black)
    for (k, pts) in DictG
        angle_rad = equi_angle(pts)
        coord_x = radius * cos(angle_rad)
        coord_y = radius * sin(angle_rad)
        scatter!(ax, [coord_x], [coord_y]; color=:red, markersize=8)
        frac = first(pts)
        label = if denominator(frac) == 1
            string(numerator(frac))
        elseif frac ∈ [Rational(-1, 0), Rational(1, 0)]
            "∞"
        else
            string(numerator(frac), "/", denominator(frac))
        end
        text!(ax, label, position=(coord_x * 1.1, coord_y * 1.1), align=(:center, :center), fontsize=12)
    end
    for e in Edges
        p1, p2 = src(e), dst(e)
        angle1 = equi_angle(DictG[p1])
        angle2 = equi_angle(DictG[p2])
        coord1 = (radius * cos(angle1), radius * sin(angle1))
        coord2 = (radius * cos(angle2), radius * sin(angle2))
        center, _ = arc_center_radius(coord1, coord2)
        f_1 = first(DictG[p1]); f_2 = first(DictG[p2])
        if f_1 in [Rational(1, 0), Rational(-1, 0)] || f_2 in [Rational(1, 0), Rational(-1, 0)]
            if f_1 in [Rational(1, 0), Rational(-1, 0)]
                if f_2 < 0
                    arc2r(ax, center, coord2, coord1; color=c)
                else
                    arc2r(ax, center, coord1, coord2; color=c)
                end
            else
                if f_1 < 0
                    arc2r(ax, center, coord2, coord1; color=c)
                else
                    arc2r(ax, center, coord1, coord2; color=c)
                end
            end
        elseif f_1 < f_2
            arc2r(ax, center, coord2, coord1; color=c)
        else
            arc2r(ax, center, coord1, coord2; color=c)
        end
    end
    hideticks!(ax)
    return ax
end

"""Plot the topograph edges (for visualization purposes)."""
function plot_topograph(ax, m::Int)
    topo_graph, topo_dict = graphtopo(m)
    edges_topo = edges(topo_graph)
    for e in edges_topo
        v1, v2 = src(e), dst(e)
        pts1 = first(topo_dict[v1])
        pts2 = first(topo_dict[v2])
        LuDrawLine(pts1, pts2)
    end
    return ax
end

# Interactive plotting function (requires a running Makie display)
"""Launch an interactive Farey graph viewer with sliders for m and variant."""
function interactive_farey_graph()
    fig = Figure(resolution = (1200, 800))
    ax = Axis(fig[1, 1]; aspect = DataAspect())
    m_slider = Slider(fig[2, 1], range = 1:10, startvalue = 1, width = 400)
    variant_toggle = Toggle(fig[3, 1], active = false)
    on(m_slider.value) do m
        if variant_toggle.active[]
            plot_farey_graph_equi(ax, m, 500.0; variant=2)
        else
            plot_farey_graph(ax, m, 500.0; variant=1)
        end
    end
    on(variant_toggle.active) do active
        if active
            plot_farey_graph_equi(ax, m_slider.value[], 500.0; variant=2)
        else
            plot_farey_graph(ax, m_slider.value[], 500.0; variant=1)
        end
    end
    return fig
end

end # module FareyDiagrams