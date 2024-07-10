# general helpers
function ssubs(expr, dict)
    if haskey(dict, expr)
        return dict[expr]
    end
    expr = Symbolics.unwrap(expr)
    if iscall(expr)
        op = ssubs(Symbolics.operation(expr), dict)
        args = map(Symbolics.arguments(expr)) do x
            ssubs(x, dict)
        end
        op(Symbolics.wrap.(args)...)
    else
        return expr
    end
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

