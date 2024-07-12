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
                rhs = map(sol -> log(sol)//log(args[1]), rhs)
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

        elseif oper === (sin)
            lhs = args[1]
            rhs = map(sol -> Symbolics.term(asin, sol), rhs)

        elseif oper === (cos)
            lhs = args[1]
            rhs = map(sol -> Symbolics.term(acos, sol), rhs)

        elseif oper === (tan)
            lhs = args[1]
            rhs = map(sol -> Symbolics.term(atan, sol), rhs)

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
           rhs = map(sol -> Symbolics.term(log, sol), rhs)
        end

        lhs = simplify(lhs)
    end
    return map(postprocess_root, rhs)
end


function attract(lhs, var)
    if n_func_occ(simplify(lhs), var) < n_func_occ(lhs, var)
        lhs = simplify(lhs)
    end

    if detect_exponential(lhs, var)
        lhs = attract_exponential(lhs, var)
    end
    if detect_addlogs(lhs, var)
        lhs = attract_logs(lhs, var)
    end

    n_func_occ(lhs, var) == 1 && return lhs

    throw("This system cannot be solved with the methods available to nl_solve. Try \
a numerical method instead.")

end

function nl_solve(lhs, var)
    nx = n_func_occ(lhs, var)
    if nx == 0
        throw("Var not present in given expression.")
    elseif nx == 1
        return isolate(lhs, var)
    elseif nx > 1
        return isolate(attract(lhs, var), var)
    end

end

# @variables x y 
# nl_solve(exp(2x)*exp(x^2 + 3) + 3, x)
# nl_solve(2sin(x+1)cos(x+1) + 1, x)
# nl_solve(x/5 + 3x^2, x)
# nl_solve(x + 2, x)
# nl_solve(expr, x)
# nl_solve(2^(x+1) + 5^(x+3), x)
 
