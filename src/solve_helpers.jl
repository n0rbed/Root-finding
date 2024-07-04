# general helpers
function clean_f(filtered_expr)
    filtered_expr = simplify_fractions(simplify(real(filtered_expr)))
    unwrapped_f = Symbolics.unwrap(filtered_expr)
    oper = 0
    try
        oper = operation(unwrapped_f)
    catch e
        return filtered_expr
    end
    if oper === (/)
        args = unsorted_arguments(unwrapped_f)
        filtered_expr = Symbolics.wrap(args[1])
        @info args[2] != 0
    end
    return filtered_expr
end

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



# multivar stuff
function contains_var(var, vars)
    for variable in vars
        if isequal(var, variable)
            return true
        end
    end
    return false
end

function add_sol!(solutions, new_sols, var, index)
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

