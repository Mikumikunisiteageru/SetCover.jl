module SetCover

using DelimitedFiles
using GLPK
using JuMP

function index(list::Vector)
	array = unique(list)
	n = length(array)
	dict = Dict(array .=> 1:n)
	n, array, dict
end

struct SCP
	n_c::Int
	counties::Vector
	county_index::Dict
	n_s::Int
	species::Vector
	species_index::Dict
	matrix::BitMatrix
	costs::Vector{<:Real}
	total_cost::T where {T <: Real}
	function SCP(occurrences, areas=nothing)
		m = size(occurrences, 1)
		@assert size(occurrences) == (m, 2)
		n_c, counties, county_index = index(occurrences[:, 1])
		n_s, species, species_index = index(occurrences[:, 2])
		matrix = falses(n_c, n_s)
		for (c, s) = eachrow(occurrences)
			i = county_index[c]
			j = species_index[s]
			matrix[i, j] = true
		end
		if isnothing(areas)
			costs = ones(n_c)
		else
			@assert allunique(areas[:, 1])
			dict_area = Dict(areas[:, 1] .=> areas[:, 2])
			@assert all(haskey(dict_area, c) for c = counties)
			costs = [dict_area[c] for c = counties]
		end
		total_cost = sum(costs)
		new(n_c, counties, county_index, n_s, species, species_index, matrix, 
			costs, total_cost)
	end
end

function Base.show(io::IO, scp::SCP)
	density = sum(scp.matrix) / length(scp.matrix)
	println(io, "Set covering problem")
	println(io, "    # Counties : ", scp.n_c)
	if scp.n_c > 5
		println(io, "      Counties : ", join(scp.counties[1:5], ", "), "...")
	else
		println(io, "      Counties : ", join(scp.counties, ", "))
	end
	println(io, "    # Species  : ", scp.n_s)
	if scp.n_s > 4
		println(io, "      Species  : ", join(scp.species[1:4], ", "), "...")
	else
		println(io, "      Species  : ", join(scp.species, ", "))
	end
	println(io, "    Density    : ", density)
	println(io, "    Total area : ", scp.total_cost)
end

struct SCPSolution
	n_c::Int
	counties::Vector
	solution::BitVector
	target::T where {T <: Real}
	total_cost::T where {T <: Real}
	function SCPSolution(scp::SCP, solution)
		n_c = scp.n_c
		counties = scp.counties
		total_cost = scp.total_cost
		target = sum(scp.costs[solution])
		new(n_c, counties, solution, target, total_cost)
	end
end

function Base.show(io::IO, sol::SCPSolution)
	println(io, "Solution of a set covering problem")
	println(io, "    Selected num  : ", sum(sol.solution), " / ", sol.n_c)
	println(io, "    Selected area : ", sol.target, " / ", sol.total_cost)
end

function complementary(scp::SCP)
	matrix = scp.matrix
	costs = scp.costs
	n, m = size(matrix)
	row_unselected = trues(n)
	col_uncovered = trues(m)
	sum_row = sum(matrix, dims=2)[:] .- costs ./ (maximum(costs) + 1)
	num_col_uncovered = m
	while num_col_uncovered > 0
		i = findmax(sum_row)[2]
		js = findall(matrix[i, :] .& col_uncovered)
		sum_row .-= sum(matrix[:, js], dims=2)[:]
		row_unselected[i] = false
		col_uncovered[js] .= false
		num_col_uncovered -= length(js)
	end
	SCPSolution(scp, .!row_unselected)
end

function greedy(scp::SCP)
	matrix = scp.matrix
	costs = scp.costs
	n, m = size(matrix)
	row_unselected = trues(n)
	col_uncovered = trues(m)
	sum_row = sum(matrix, dims=2)[:]
	num_col_uncovered = m
	while num_col_uncovered > 0
		i = findmax(sum_row ./ costs)[2]
		js = findall(matrix[i, :] .& col_uncovered)
		sum_row .-= sum(matrix[:, js], dims=2)[:]
		row_unselected[i] = false
		col_uncovered[js] .= false
		num_col_uncovered -= length(js)
	end
	SCPSolution(scp, .!row_unselected)
end

function lagrangian(scp::SCP; _smallnumber=1e-5)
	# THE IMPLEMENTATION OF THIS FUNCTION HAS REFERRED TO  SetCoverPy ,
	# A PYTHON PACKAGE BY  Guangtun Ben Zhu
	matrix = scp.matrix
	costs = scp.costs
	s = any(matrix[:, sum(matrix, dims=1)[:] .== 1], dims=2)[:]
	iuncovered = .!any(matrix[s, :], dims=1)[:]
	u = zeros(size(matrix, 2))
	cost_inv = sum(matrix, dims=2)[:] ./ costs
	t = maximum(cost_inv .* matrix[:, iuncovered], dims=1)[:]
	u[iuncovered] = 1. ./ max.(t, minimum(cost_inv))
	while any(iuncovered)
		id = findall(.!s)
		mu = max.(_smallnumber, sum(matrix[id, iuncovered], dims=2)[:])
		gamma = costs[id] .- (matrix * u)[id]
		foo(ga, mu) = ga .>= 0 ? ga / mu : ga * mu
		inewcolumn = id[findmin(foo.(gamma, mu))[2]]
		s[inewcolumn] = true
		new = matrix[inewcolumn, :]
		iuncovered[new] .= false
		u[new] .= 0
	end
	SCPSolution(scp, s)
end

function ilp_solver(scp::SCP; optimizer=GLPK.Optimizer())
	matrix = scp.matrix
	costs = scp.costs
	n, m = size(matrix)
	model = direct_model(optimizer)
	@variable(model, x[1:n], binary=true)
	@constraint(model, matrix' * x .>= 1)
	@objective(model, Min, costs' * x)
	optimize!(model)
	SCPSolution(scp, Bool.(value.(x)))
end

function write_compact(file_name::AbstractString, sol::SCPSolution)
	writedlm(file_name, sol.counties[sol.solution])
end

function write_full(file_name::AbstractString, sol::SCPSolution, dlm='\t')
	writedlm(file_name, [sol.counties Int.(sol.solution)], dlm)
end

export SCP, SCPSolution
export complementary, greedy, lagrangian, ilp_solver
export target, write_compact, write_full

end # module
