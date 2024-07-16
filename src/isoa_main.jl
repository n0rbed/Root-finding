using Symbolics

function isolate(lhs, var)
    rhs = Vector{Any}([0])
    lhs = Symbolics.unwrap(lhs)
    while !isequal(lhs, var)
        subs, poly = filter_poly(lhs, var)
        if check_poly_inunivar(poly, var)
            roots = []
            for i = 1:length(rhs)
                append!(roots, solve(Symbolics.wrap(lhs-rhs[i]), var))
            end
            return roots
        end

        oper = Symbolics.operation(lhs)
        args = arguments(lhs)

        if oper === (+)
            for arg in args
                vars = Symbolics.get_variables(arg)
                if any(isequal(x, var) for x in vars)
                    continue
                end
                lhs = lhs - arg
                rhs = map(sol -> sol - arg, rhs)
            end

        elseif oper === (*)
            for arg in args
                vars = Symbolics.get_variables(arg)
                if any(isequal(x, var) for x in vars)
                    continue
                end
                lhs = lhs / arg
                rhs = map(sol -> sol/arg, rhs)
            end

        # TODO: add / oper
        elseif oper === (/)
            var_innumerator = any(isequal(x, var) for x in Symbolics.get_variables(args[1]))
            if var_innumerator
                # x / 2 = y
                lhs = args[1]
                rhs = map(sol -> sol*args[2], rhs)
            else
                # 2 / x = y
                lhs = args[2]
                rhs = map(sol -> args[1]//sol, rhs)
            end



        elseif oper === (^)
            if any(isequal(x, var) for x in Symbolics.get_variables(args[1])) && all(!isequal(x, var) for x in Symbolics.get_variables(args[2]))
                lhs = args[1]
                rhs = map(sol -> Symbolics.term(^, sol, 1//args[2]), rhs)
                if args[2] isa Int64 && iseven(args[2])
                    append!(rhs, map(sol -> -sol, rhs))
                end
            else
                lhs = args[2]
                rhs = map(sol -> Symbolics.term(/, Symbolics.term(slog, sol), Symbolics.term(slog, args[1])), rhs)
            end

        elseif oper === (log)
            lhs = args[1]
            rhs = map(sol -> Symbolics.term(^, Base.MathConstants.e, sol), rhs)

        elseif oper === (log2)
            lhs = args[1]
            rhs = map(sol -> Symbolics.term(^, 2, sol), rhs)

        elseif oper === (log10)
            lhs = args[1]
            rhs = map(sol -> Symbolics.term(^, 10, sol), rhs)

        elseif oper === (sqrt)
            lhs = args[1]
            rhs = map(sol -> Symbolics.term(^, sol, 2), rhs)

        elseif oper === (cbrt)
            lhs = args[1]
            rhs = map(sol -> Symbolics.term(^, sol, 3), rhs)

        elseif oper === (sin) || oper === (cos) || oper === (tan)
            rev_oper = Dict(sin=>asin, cos=>acos, tan=>atan)
            lhs = args[1]
            new_var = gensym()
            new_var = (@variables $(new_var))[1]
            rhs = map(sol -> Symbolics.term(rev_oper[oper], sol) + Symbolics.term(*, Base.MathConstants.pi, 2*new_var), rhs)
            @info string(new_var) * " ϵ" * " Ζ: e.g. 0, 1, 2..."

        elseif oper === (asin)
            lhs = args[1]
            rhs = map(sol -> Symbolics.term(sin, sol), rhs)

        elseif oper === (acos)
            lhs = args[1]
            rhs = map(sol -> Symbolics.term(cos, sol), rhs)

        elseif oper === (atan)
            lhs = args[1]
            rhs = map(sol -> Symbolics.term(tan, sol), rhs)
        elseif oper === (exp)
           lhs = args[1] 
           rhs = map(sol -> Symbolics.term(slog, sol), rhs)
        end

        lhs = simplify(lhs)
    end
    return map(postprocess_root, rhs)
end


function attract(lhs, var)
    if n_func_occ(simplify(lhs), var) <= n_func_occ(lhs, var)
        lhs = simplify(lhs)
    end

    if detect_exponential(lhs, var)
        lhs = attract_exponential(lhs, var)
    end
    if detect_addlogs(lhs, var)
        lhs = attract_logs(lhs, var)
    end
    lhs = attract_trig(lhs, var)

    n_func_occ(lhs, var) == 1 && return isolate(lhs, var)

    lhs, sub = turn_to_poly(lhs, var)
    if (isequal(sub, Dict()) || n_func_occ(lhs, collect(keys(sub))[1]) != 1) 
        throw("This system cannot be solved with the methods available to nl_solve. Try \
        a numerical method instead.")
    end

    new_var = collect(keys(sub))[1]
    new_var_val = collect(values(sub))[1]

    roots = isolate(lhs, new_var)
    new_roots = []
    for root in roots
        new_sol = isolate(new_var_val - root, var)
        push!(new_roots, new_sol)
    end
    new_roots = collect(Iterators.flatten(new_roots))

    return new_roots
end

function nl_solve(lhs, var)
    nx = n_func_occ(lhs, var)
    if nx == 0
        throw("Var not present in given expression.")
    elseif nx == 1
        return isolate(lhs, var)
    elseif nx > 1
        return attract(lhs, var)
    end

end

# @variables x y 
# nl_solve(lhs, x)
# attract(9^x + 3^x + 1, x)