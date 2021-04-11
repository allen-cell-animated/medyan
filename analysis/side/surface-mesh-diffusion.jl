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

# ╔═╡ 1b83f02a-dc6a-4124-9b7f-eaaef5edb95d
begin
	frames = read_vertices(joinpath(@__DIR__, "../../build/side-surface-mesh-diffusion.txt"))
    max_copy = 10000
    β = 1 / 4.1

    total_area = 480 * 480
    areas = [v.cell_area for v ∈ frames[1].vertices]
    xs = [v.coord[1] for v ∈ frames[1].vertices]
    ys = [v.coord[2] for v ∈ frames[1].vertices]

    # calculate expected copy numbers
    expected_uniform = [max_copy * v.cell_area / total_area for v ∈ frames[1].vertices]
    expected_sin = [v.cell_area * exp(β * 5 * sin(2π / 1000 * v.coord[1])) for v ∈ frames[1].vertices]
    expected_sin *= max_copy / sum(expected_sin)
	
	expected_conc = [expected_uniform, expected_uniform, expected_sin]
end

# ╔═╡ f7e0efbd-2370-48af-bb4a-929d7fcb9d5d
md"
Species index $(@bind species_index Slider(1:3; show_value=true))
"

# ╔═╡ 2e4b60fd-e59c-42e0-b68e-9da50f73c09c

begin
	function make_plot()
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

	p = make_plot()
end

# ╔═╡ Cell order:
# ╠═4e9478cb-11d3-4d00-8456-8e1e1e20178e
# ╠═1ef4e7be-b0ca-4d37-9b87-416495802473
# ╠═78e4be62-f030-4e3b-8e5c-a1bff309e735
# ╟─00f87b79-d447-4732-9ebf-4378dd804905
# ╠═1b83f02a-dc6a-4124-9b7f-eaaef5edb95d
# ╠═f7e0efbd-2370-48af-bb4a-929d7fcb9d5d
# ╠═2e4b60fd-e59c-42e0-b68e-9da50f73c09c
