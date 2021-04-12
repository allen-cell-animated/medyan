### A Pluto.jl notebook ###
# v0.14.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 4e9478cb-11d3-4d00-8456-8e1e1e20178e
begin
	using Plots
	using PlutoUI
	using StaticArrays
	using StatsBase
	
	gr()
	
	src_dir = joinpath(@__DIR__, "../../build")
	work_dir = @__DIR__
end

# ╔═╡ 1ef4e7be-b0ca-4d37-9b87-416495802473
struct VertexData
    cell_area :: Float64

    coord :: SVector{3, Float64}
    species_count :: SVector{3, Float64}
end

# ╔═╡ 78e4be62-f030-4e3b-8e5c-a1bff309e735
struct FrameData
    time :: Float64
    vertices :: Vector{VertexData}
end

# ╔═╡ 00f87b79-d447-4732-9ebf-4378dd804905
function read_vertices(filepath::AbstractString)
    frames = Vector{FrameData}(undef, 0)
    open(filepath) do f
        for line ∈ eachline(f)
            splitted = split(line)
            if splitted[1] == "FRAME"
                push!(frames, FrameData(parse(Float64, splitted[2]), Vector{VertexData}(undef, 0)))
            else
                vertices = frames[end].vertices
                push!(
                    vertices,
                    VertexData(
                        parse(Float64, splitted[2]),
                        SVector{3}([parse(Float64, splitted[3]), parse(Float64, splitted[4]), parse(Float64, splitted[5])]),
                        SVector{3}([parse(Int, splitted[6]), parse(Int, splitted[7]), parse(Int, splitted[8])])
                    )
                )
            end
        end
    end
    frames
end

# ╔═╡ 18500426-5217-4784-b6e3-9e04030270bc
function read_triangles(filepath::AbstractString)
	open(filepath) do f
		num_triangles = parse(Int, readline(f))
		# The file uses 0-based indexing, and must be converted to 1-based indexing here.
		inds = map(x->parse(Int, x) + 1, f |> readline |> split)
		vertex_indices = reshape(inds, 3, num_triangles)
	end
end

# ╔═╡ 1b83f02a-dc6a-4124-9b7f-eaaef5edb95d
begin
	frames = read_vertices(joinpath(src_dir, "side-surface-mesh-diffusion.txt"))
	triangle_vertex_indices = read_triangles(joinpath(src_dir, "side-surface-mesh-diffusion-triangle.txt"))

    max_copy = 10000
    β = 1 / 4.1

    total_area = 480 * 480
    areas = [v.cell_area for v ∈ frames[1].vertices]
    xs = [v.coord[1] for v ∈ frames[1].vertices]
    ys = [v.coord[2] for v ∈ frames[1].vertices]

	energy_at_coord(coord) = -5 * sin(2π / 1000 * coord[1])

    # calculate expected copy numbers
    expected_uniform = [max_copy * v.cell_area / total_area for v ∈ frames[1].vertices]
    expected_sin = [v.cell_area * exp(-β * energy_at_coord(v.coord)) for v ∈ frames[1].vertices]
    expected_sin *= max_copy / sum(expected_sin)
	
	expected_conc = [expected_uniform, expected_uniform, expected_sin]
end

# ╔═╡ f7e0efbd-2370-48af-bb4a-929d7fcb9d5d
md"
Species index $(@bind species_index Slider(1:3; show_value=true))
"

# ╔═╡ 2e4b60fd-e59c-42e0-b68e-9da50f73c09c

function make_avg_plot()
	frame_range = 10:length(frames)
	conc_avg = zeros(length(frames[1].vertices))

	for frame_index ∈ frame_range
		conc_avg += [v.species_count[species_index] / v.cell_area for v ∈ frames[frame_index].vertices]
	end
	conc_avg /= length(frame_range)

	# avg copy numbers
	conc_plot = scatter3d(
		xs, ys, conc_avg;
		label = "average frame $frame_range",
		xlabel = "x (nm)",
		ylabel = "y (nm)",
		zlabel = "concentration (nm^-2)",
	)

	scatter3d!(conc_plot, xs, ys, expected_conc[species_index] ./ areas; label = "expected")
	conc_plot
end

# ╔═╡ 5b4da84b-f765-4e49-9251-0c2426985ff0
begin
	avg_plot = make_avg_plot()
	savefig(avg_plot, joinpath(work_dir, "surface-mesh-diffusion-avg.png"))
	avg_plot
end

# ╔═╡ 61623a45-1616-44d0-b468-0658a13fb49c
begin
	anim = @animate for (fi, frame) ∈ enumerate(frames[1:20])
		conc_plot = scatter3d(
			xs, ys, [v.species_count[species_index] / v.cell_area for v ∈ frame.vertices];
			label = "frame $fi",
			xlabel = "x (nm)",
			ylabel = "y (nm)",
			zlabel = "concentration (nm^-2)",
		)

		scatter3d!(conc_plot, xs, ys, expected_conc[species_index] ./ areas; label = "expected")
	end
	gif(anim, joinpath(work_dir, "surface-mesh-diffusion-anim.gif"); fps=4)
end

# ╔═╡ 5e0f3cd3-b0b6-46c5-a812-f1edc0c56e9e
begin
	# Make mesh graph with energy
	p = plot(;
		xlabel = "x (nm)",
		ylabel = "y (nm)",
		xlim = (0, 500),
		ylim = (0, 500),
		aspect_ratio = :equal,
	)
	
	# Energy part
	heat_xrange = 0:5:500
	heat_yrange = 0:100:500
	heatmap!(p,
		heat_xrange, heat_yrange,
		(x, y) -> energy_at_coord(@SVector [x, y, 0]);
		color = :copper,
	)

	# Mesh part
	for ti ∈ 1:size(triangle_vertex_indices)[2]
		vs = triangle_vertex_indices[[1,2,3,1], ti]
		xs = map(vi -> frames[1].vertices[vi].coord[1], vs)
		ys = map(vi -> frames[1].vertices[vi].coord[2], vs)
		plot!(p, xs, ys; legend=:none, color=:white, colorbar_title="energy (pN nm)")
	end

	savefig(p, joinpath(work_dir, "surface-mesh-diffusion-setup.png"))
	p
end

# ╔═╡ Cell order:
# ╠═4e9478cb-11d3-4d00-8456-8e1e1e20178e
# ╠═1ef4e7be-b0ca-4d37-9b87-416495802473
# ╠═78e4be62-f030-4e3b-8e5c-a1bff309e735
# ╟─00f87b79-d447-4732-9ebf-4378dd804905
# ╟─18500426-5217-4784-b6e3-9e04030270bc
# ╠═1b83f02a-dc6a-4124-9b7f-eaaef5edb95d
# ╟─f7e0efbd-2370-48af-bb4a-929d7fcb9d5d
# ╟─2e4b60fd-e59c-42e0-b68e-9da50f73c09c
# ╠═5b4da84b-f765-4e49-9251-0c2426985ff0
# ╠═61623a45-1616-44d0-b468-0658a13fb49c
# ╠═5e0f3cd3-b0b6-46c5-a812-f1edc0c56e9e
