using Symbolics, Groebner, SymbolicUtils
include("univar.jl")
include("coeffs.jl")
include("nemo_stuff.jl")

function get_and_sub_factors(subs, filtered_expr, subbed_factors)
    factors = []
    @variables I
    if any(isequal(I, var) for var in collect(keys(subs)))
        u, factors = factor_use_nemo_and_split(filtered_expr, I)
        delete!(subs, I)
    else
        factors = deepcopy(subbed_factors)
    end


    # sub into factors 
    for i = 1:length(factors)
        factors[i] = Symbolics.substitute(factors[i], subs, fold=false)
    end

    return factors
end

struct RootsOf
    poly::Num
end
Base.show(io::IO, r::RootsOf) = print(io, "roots_of(", r.poly, ")")


function solve(expression, x, mult=false)
    args = []
    mult_n = 1
    try
        exp = Symbolics.unwrap(simplify(expression))
        args = unsorted_arguments(exp)
        operation = SymbolicUtils.operation(exp)
        if isequal(operation, ^) && args[2] isa Int64
            expression = Symbolics.wrap(args[1])
            mult_n = args[2]
        end
    catch e
    end
    expression = simplify(expression)

    subs, filtered_expr = filter_poly(expression, x)
    filtered_expr = simplify(real(filtered_expr))

    degree = Symbolics.degree(filtered_expr, x)
    u, subbed_factors = factor_use_nemo(filtered_expr)
    subbed_factors = convert(Vector{Any}, subbed_factors)
    factors = get_and_sub_factors(subs, filtered_expr, subbed_factors)

    arr_roots = []

    if degree < 5 && length(factors) == 1
        arr_roots = get_roots(expression, x)

        # multiplicities
        if mult
            for i = 1:(mult_n-1)
                try
                    push!(arr_roots, arr_roots[1])    
                catch e
                end
            end
        end

        return arr_roots
    end

    if length(factors) != 1
        @assert isequal(expand(filtered_expr - u*expand(prod(subbed_factors))), 0)

        for factor in factors
            roots = solve(factor, x, mult)
            if isequal(typeof(roots), RootsOf)
                push!(arr_roots, roots)
            else
                append!(arr_roots, roots)
            end
        end
    end


    if isequal(arr_roots, [])
        return RootsOf(Symbolics.wrap(expression))
    end

    return arr_roots
end

# You can compute the GCD between a system of polynomials by doing the following:
# Get the GCD between the first two polys,
# and get the GCD between this result and the following index,
# say: solve([x^2 - 1, x - 1, (x-1)^20], x)
# the GCD between the first two terms is obviously x-1,
# now we call gcd_use_nemo() on this term, and the following,
# gcd_use_nemo(x - 1, (x-1)^20), which is again x-1.
# now we just need to solve(x-1, x) to get the common root in this
# system of equations.
function solve(polys::Vector, x::Num, mult=false)
    polys = unique(polys)

    if length(polys) < 1
        throw("No expressions entered")
    end
    if length(polys) == 1
        return solve(polys[1], x, mult)
    end

    gcd = gcd_use_nemo(polys[1], polys[2])

    for i = 3:length(polys)
        gcd = gcd_use_nemo(gcd, polys[i])
    end
    
    if isequal(gcd, 1)
        @info "Nemo gcd is 1."
        return []
    end

    return solve(gcd, x, mult)
end

function contains(var, vars)
    for variable in vars
        if isequal(var, variable)
            return true
        end
    end
    return false
end

function add_sol(solutions, new_sols, var, index)
    sol_used = solutions[index]
    deleteat!(solutions, index)
    for new_sol in new_sols
        sol_used[var] = new_sol
        push!(solutions, deepcopy(sol_used))
    end
    return solutions
end

function add_sol_to_all(solutions, new_sols, var)
    existing_solutions = deepcopy(solutions)
    solutions = []
    for new_sol in new_sols
        copy_sol = deepcopy(existing_solutions)
        for i = 1:length(copy_sol)
            copy_sol[i][var] = new_sol
        end
        append!(solutions, copy_sol)
    end
    return solutions
end

function solve(eqs::Vector{Num}, vars::Vector{Num}, mult=false)
    # do the trick
    @variables HAT
    old_len = length(vars)
    push!(vars, HAT)
    new_eqs = []
    generating = true
    while generating
        new_eqs = deepcopy(eqs)
        eq = HAT
        for i = 1:(old_len)
            eq -= rand(1:10)*vars[i]
        end
        push!(new_eqs, eq)
        new_eqs = convert(Vector{Any}, Symbolics.groebner_basis(new_eqs, ordering=Lex(vars)))

        if length(new_eqs) <= length(vars) 
            generating &= false
        end

        for i = 2:length(new_eqs)
            generating |= all(Symbolics.degree(var) > 1 for var in Symbolics.get_variables(new_eqs[i]))
        end
    end

    solutions = []

    # handle "unsolvable" cases
    if isequal(1, new_eqs[1])
        return solutions
    end
    if length(new_eqs) < length(vars)
        throw("Infinite number of solutions")
    end


    # first, solve any single variable equations
    i = 1
    while !(i > length(new_eqs))
            present_vars = Symbolics.get_variables(new_eqs[i])
        for var in vars
            if size(present_vars, 1) == 1 && isequal(var, present_vars[1])
                new_sols = solve(Symbolics.wrap(new_eqs[i]), var, mult)

                if length(solutions) == 0
                    append!(solutions, [Dict{Num, Any}(var => sol) for sol in new_sols])
                else
                    solutions = add_sol_to_all(solutions, new_sols, var)
                end

                deleteat!(new_eqs, i)
                i = i - 1
                break
            end
        end
        i = i + 1
    end


    # second, iterate over eqs and sub each found solution
    # then add the roots of the remaining unknown variables 
    for eq in new_eqs
        solved = false
        present_vars = Symbolics.get_variables(eq)
        size_of_sub = length(solutions[1])

        if size(present_vars, 1) <= (size_of_sub + 1)
            while !solved 
                subbed_eq = eq
                for (var, root) in solutions[1]
                    subbed_eq = Symbolics.substitute(subbed_eq, Dict([var => root]), fold=false)
                end

                var_tosolve = Symbolics.get_variables(subbed_eq)[1]
                new_var_sols = solve(subbed_eq, var_tosolve, mult)
                solutions = add_sol(solutions, new_var_sols, var_tosolve, 1)

                solved = all(x -> length(x) == size_of_sub+1, solutions)
            end
        end
    end

    for roots in solutions
        delete!(roots, HAT)
    end
    return solutions
end
    

#@variables x y z
#solve(x^4 + sqrt(complex(-8//1)), x)
@variables x
exp = x^3 + Symbolics.term(sqrt, Symbolics.term(complex, -2//1))*x + 2
get_roots_deg3(exp, x)