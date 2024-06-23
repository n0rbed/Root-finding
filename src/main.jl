using Symbolics, Groebner, SymbolicUtils
include("univar.jl")
include("coeffs.jl")
include("nemo_stuff.jl")


function solve(expression, x)
    args = []
    m = 1
    try
        args = Symbolics.arguments(expression.val)
        if isequal(SymbolicUtils.operation(expression.val), ^) && SymbolicUtils.arguments(expression.val)[2] isa Int64
            expression = Symbolics.wrap(args[1])
            m = args[2]
        end
    catch e
    end
    

    expression = expand(expression)
    expression = simplify.(expression)
    degree = Symbolics.degree(expression, x)

    subs, filtered_expression = filter_poly(expression, x)
    u, factors = factor_use_nemo(filtered_expression)

    # sub into factors 
    for i = 1:length(factors)
        for (var, sub) in subs 
            factors[i] = Symbolics.substitute(factors[i], Dict([var => sub]), fold=false)
        end
    end

    arr_roots = []

    if degree < 5 && length(factors) == 1
        arr_roots = get_roots(expression, x)
        #sub_roots(arr_roots, subs)
        for i = 1:(m-1)
            try
                push!(arr_roots, arr_roots[1])    
            catch e
            end
        end
        return arr_roots
    end

    if length(factors) != 1
        @assert isequal(expand(expression - u*expand(prod(factors))), 0)

        for factor in factors
            append!(arr_roots, solve(factor, x))
        end
    end


    if isequal(arr_roots, [])
        throw("This expression does not have an exact solution, use a numerical method instead.")
    end

    # is this necessary?
    #sub_roots(arr_roots, subs)
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
function solve(polys::Vector, x::Num)
    polys = unique(polys)

    if length(polys) < 1
        throw("No expressions entered")
    end
    if length(polys) == 1
        return solve(polys[1], x)
    end

    gcd = gcd_use_nemo(polys[1], polys[2])

    for i = 3:length(polys)
        gcd = gcd_use_nemo(gcd, polys[i])
    end
    
    if isequal(gcd, 1)
        @info "Nemo gcd is 1."
        return []
    end
    return solve(gcd, x)
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

function solve(eqs::Vector{Num}, vars::Vector{Num})
    # do the trick
    @variables HAT
    generated = false
    old_vars = deepcopy(vars)
    push!(vars, HAT)
    new_eqs = []
    while !generated
        new_eqs = deepcopy(eqs)
        eq = HAT
        for i = 1:(length(old_vars))
            eq -= rand(1:10)*vars[i]
        end
        push!(new_eqs, eq)
        new_eqs = convert(Vector{Any}, Symbolics.groebner_basis(new_eqs, ordering=Lex(vars)))

        if length(new_eqs) <= length(vars) 
            generated = true
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
                new_sols = solve(Symbolics.wrap(new_eqs[i]), var)

                if length(solutions) == 0
                    append!(solutions, [Dict(var => sol) for sol in new_sols])
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
    j = 1
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
                subbed_eq = Symbolics.wrap(subbed_eq)


                var_tosolve = Symbolics.get_variables(subbed_eq)[1]
                new_var_sols = solve(subbed_eq, var_tosolve)
                solutions = add_sol(solutions, new_var_sols, var_tosolve, 1)

                solved = all(x -> length(x) == j+1, solutions)
            end
        end
        if solved
            j = j + 1
        end
    end

    for roots in solutions
        delete!(roots, HAT)
    end
    return solutions
end
    

@variables x y z
eqs = [x+y^2+z, z*x*y, z+3x+y]
solve(eqs, [x,y,z])
