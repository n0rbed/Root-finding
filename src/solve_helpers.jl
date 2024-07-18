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

function ssqrt(n)
    n = Symbolics.unwrap(n)

    if n isa Real
        n > 0 && return sqrt(n)
        return sqrt(complex(n))
    end

    if n isa Complex 
        return sqrt(n)
    end

    if n isa SymbolicUtils.BasicSymbolic{Real}
        return Symbolics.term(ssqrt, n)
    end

end


function scbrt(n)
    n = Symbolics.unwrap(n)

    if n isa Real
        return cbrt(n)
    end

    if n isa Complex
        return (n)^(1/3)
    end

    if n isa SymbolicUtils.BasicSymbolic{Real}
        return Symbolics.term(scbrt, n)
    end

end

function slog(n)
    n = Symbolics.unwrap(n)

    if n isa Real
        return n > 0 ? log(n) : log(complex(n))
    end

    if n isa Complex
        return log(n)
    end

    if n isa SymbolicUtils.BasicSymbolic{Real}
        return Symbolics.term(slog, n)
    end
end


struct RootsOf
    poly::Num
    var::Num
end

Base.show(io::IO, r::RootsOf) = print(io, "roots_of(", r.poly, ", ", x, ")")


function check_expr_validity(expr)
    type_expr = typeof(expr)

    if type_expr == Num || type_expr == SymbolicUtils.BasicSymbolic{Real} || 
    type_expr == Complex{Num} || type_expr == Symbolics.ComplexTerm{Real}
        return true
    end
    @assert isequal(expr, 0) "Invalid input"
end



function f_numbers(n)
    if n isa Num || n isa SymbolicUtils.BasicSymbolic ||
     n isa Complex{Num} || n isa Symbolics.ComplexTerm || n isa Float64
        return n
    end

    if n isa Integer
        n = abs(n) > 1000 && n isa Integer ? BigInt(n) : n
        return n
    end

    if n isa Complex
        real_part = f_numbers(n.re)
        im_part = f_numbers(n.im)
        return real_part + im_part*im
    end

    if n isa Rational && n isa Real
        top = numerator(n) > 100 ? BigInt(numerator(n)) : numerator(n)
        bottom = denominator(n) > 100 ? BigInt(denominator(n)) : denominator(n)

        return top//bottom
    end

end

function comp_rational(x,y)
    try
        r = x//y
        return r
    catch e
        r = nothing
        if x isa ComplexF64
            real_p = real(x)
            imag_p = imag(x)
            r = Rational{BigInt}(real_p)//y
            if !isequal(imag_p, 0)
                r += (Rational{BigInt}(imag_p)//y)*im
            end
        elseif x isa Float64 
            r = Rational{BigInt}(x)//y
        end
        return r
    end
end



### multivar stuff ###
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

