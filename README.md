# FareyDiagrams.jl

**FareyDiagrams.jl** is a Julia library for constructing, exploring, and visualizing **Farey diagrams** and their associated **topographs**. It supports plotting using Makie.jl as well as computing triangle centers using both standard and equidistant angular embeddings.

> This library can generate both the classic Farey diagram and a variant that reveals extra structure when `n = 2`.

---

## Features

- Constructs Farey graphs of any integer order `n = {1,2}`.
- Supports two graph types:
  - `FareyGraph(n)`: the classic Farey diagram built via mediants.
  - `FareyGraph2(n)`: an extended variant that adds "extra" rational points (e.g. ±2/1, ±1/2) and reveals new symmetries at `n = 2`.
- Triangle generation and tracking: identify new triangles formed at each level.
- Computes **centroids** of triangles (used in topographs).
- Embeds vertices on the unit circle using:
  - **Standard angle mapping** (`angle`): based on rational number geometry.
  - **Equidistant angle mapping** (`equi_angle`): based on continued fractions.
- Arc-based plotting using Luxor and Makie.
- **interactive viewer** with sliders for level and variant.

---

## TO-DO 
- [ ] look for the cases for n>2 
- [ ] ...


---

##  Example

```julia
using FareyDiagrams

# Build the Farey graph of order 2
G, vertex_dict, triangles = FareyGraph(2)

# Or use the extended variant
G2, vertex_dict2 = FareyGraph2(2)

# Plot (requires GLMakie)
using GLMakie

fig = Figure(resolution=(800, 800))
ax = Axis(fig[1, 1], aspect=DataAspect())
plot_farey_graph(ax, 2, 500.0; variant=1)  # change to variant=2 for FareyGraph2
display(fig)
