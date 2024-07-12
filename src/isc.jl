using Symbolics

function isolate(lhs, var)
    rhs = [0]
    lhs = Symbolics.unwrap(lhs)
    while !isequal(lhs, var)
        try
            subs, poly = filter_poly(lhs-rhs, var)
            if check_polynomial(Symbolics.wrap(poly)) && n_occurrences(poly, var) > 1
                return solve(Symbolics.wrap(lhs-rhs), var)
            end
        catch e
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


function attract_logs(lhs, var)
    contains_var(arg) = occursin(string(var), string(arg))

    r_addlogs = Vector{Any}() 
    push!(r_addlogs, @acrule log(~x::(contains_var)) + log(~y::(contains_var)) => log(~x * ~y))
    push!(r_addlogs, @acrule ~z*log(~x::(contains_var)) + log(~y::(contains_var)) => log((~x)^(~z) * ~y))
    push!(r_addlogs, @acrule ~z*log(~x::(contains_var)) + ~h*log(~y::(contains_var)) => log((~x)^(~z) * (~y)^(~h)))


    for r in r_addlogs
        lhs = SymbolicUtils.Fixpoint(r)(Symbolics.unwrap(lhs))
    end

    return lhs
end


function attract_exponential(lhs, var)
    contains_var(arg) = occursin(string(var), string(arg))

    r_addexpon = Vector{Any}()
    push!(r_addexpon, @acrule (~b)^(~f::(contains_var)) + (~d)^(~g::(contains_var)) => ~f*Symbolics.term(log, ~b) - ~g*Symbolics.term(log, ~d) + Symbolics.term(log, Symbolics.term(complex, -1)))
    push!(r_addexpon, @acrule (~a)*(~b)^(~f::(contains_var)) + (~d)^(~g::(contains_var)) => Symbolics.term(log, -~a) + ~f*Symbolics.term(log, ~b) - ~g*Symbolics.term(log, ~d))
    push!(r_addexpon, @acrule (~a)*(~b)^(~f::(contains_var)) + (~c)*(~d)^(~g::(contains_var)) => Symbolics.term(log, -(~a)/(~c)) + ~f*Symbolics.term(log, ~b) - ~g*Symbolics.term(log, ~d))

    for r in r_addexpon
        lhs = SymbolicUtils.Fixpoint(r)(Symbolics.unwrap(lhs))
    end

    return expand(lhs)
end


function attract(lhs, var)
    unwrapped_lhs = Symbolics.unwrap(lhs)

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

function n_func_occ(expr, var)
    expr = Symbolics.unwrap(expr)
    !iscall(expr) && return n_occurrences(expr, var)
    args, cur_oper = Symbolics.arguments(expr), Symbolics.operation(expr)
    counted_ops = [sin, log, log2, log10, cos, tan, asin, acos, atan, exp]
    n = 0


    if cur_oper === (*) || cur_oper === (+)

        outside = false
        for arg in args
            n_occurrences(arg, var) == 0 && continue

            if !iscall(arg) && isequal(var, Symbolics.get_variables(expr)[1]) && !outside
                n += 1
                outside = true
                continue
            end
            !iscall(arg) && continue
            oper = Symbolics.operation(arg)

            args_arg = Symbolics.arguments(arg)
            oper_arg = Symbolics.operation(arg)
            case_1_pow = oper_arg === (^) && n_occurrences(args_arg[2], var) == 0  && n_occurrences(args_arg[1], var) != 0 && !check_poly_inunivar(args_arg[1], var)
            case_2_pow = oper_arg === (^) && n_occurrences(args_arg[2], var) != 0  && n_occurrences(args_arg[1], var) == 0  
            is_var_outside(arg) = check_poly_inunivar(arg, var) && !outside && n_occurrences(arg, var) != 0

            # any transcedental operation and the case:  (weird_transcedental_f(x))^(something)
            if any(isequal(oper, op) for op in counted_ops) || case_1_pow
                n += n_func_occ(args_arg[1], var)
            
            # the case (some constant)^(f(x))
            elseif case_2_pow
                n += n_func_occ(args_arg[2], var) 
            
                
            # var is outside 'x'+1
            elseif is_var_outside(arg)
                n += 1
                outside = true

            # n(2 / x) = 1; n(x/x^2) = 2?
            elseif oper_arg === (/)
                n += n_func_occ(args_arg[1], var)
                n += n_func_occ(args_arg[2], var)

            # multiplication cases
            elseif oper_arg === (*)
                args_arg = arguments(arg)

                for sub_arg in args_arg
                    # x*log(2)
                    if is_var_outside(sub_arg)
                        n += 1
                        outside = true
                    # log(x)*y
                    elseif !check_poly_inunivar(sub_arg, var)
                        n += n_func_occ(sub_arg, var)
                    end
                end
            end
        end

    else
        for arg in args
            n += n_func_occ(arg, var)
        end
    end
    
    return n
end

function n_occurrences(expr, var)
    n = 0
    !iscall(Symbolics.unwrap(expr)) && any(isequal(var, x) for x in Symbolics.get_variables(expr)) && return 1
    !iscall(Symbolics.unwrap(expr)) && return 0 

    args = Symbolics.arguments(Symbolics.unwrap(expr))

    for arg in args
        n += traverse(arg, var)
    end

    return n
end

function traverse(argument, var)
    args = []
    try
        args = Symbolics.arguments(argument)
    catch e
        if isequal(argument, var)
            return 1 
        end
        return 0
    end

    n = 0

    for arg in args
        n += traverse(arg, var)
    end
    return n
end



@variables x y 
nl_solve(log(x) - 1, x)
# nl_solve(x/5 + 3x^2, x)
# nl_solve(x + 2, x)
# nl_solve(expr, x)
# nl_solve(2^(x+1) + 5^(x+3), x)
 
