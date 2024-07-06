using Symbolics

function isolate(lhs, var)
    rhs = 0
    lhs = Symbolics.unwrap(lhs)
    while !isequal(lhs, var)
        try
            subs, poly = filter_poly(lhs-rhs, var)
            if check_polynomial(Symbolics.wrap(poly))
                return solve(Symbolics.wrap(lhs-rhs), var)
            end
        catch e
        end

        oper = Symbolics.operation(lhs)
        args = unsorted_arguments(lhs)

        if oper === (+)
            for arg in args
                vars = Symbolics.get_variables(arg)
                if any(isequal(x, var) for x in vars)
                    continue
                end
                lhs = lhs - arg
                rhs = rhs - arg
            end

        elseif oper === (*)
            for arg in args
                vars = Symbolics.get_variables(arg)
                if any(isequal(x, var) for x in vars)
                    continue
                end
                lhs = lhs / arg
                rhs = rhs / arg
            end

        elseif oper === (^)
            if any(isequal(x, var) for var in Symbolics.get_variables(args[1]))
                lhs = Symbolics.term(^, lhs, (1/args[1]))
                rhs = Symbolics.term(^, rhs, (1/args[1]))
            else
                lhs = args[2]
                rhs = log(rhs) / log(args[1])
            end

        elseif oper === (log)
            lhs = args[1]
            rhs = Symbolics.term(^, Base.MathConstants.e, rhs)

        elseif oper === (log2)
            lhs = args[1]
            rhs = Symbolics.term(^, 2, rhs)

        elseif oper === (log10)
            lhs = args[1]
            rhs = Symbolics.term(^, 10, rhs)

        elseif oper === (sqrt)
            lhs = args[1]
            rhs = Symbolics.term(^, rhs, 2)

        elseif oper === (cbrt)
            lhs = args[1]
            rhs = Symbolics.term(^, rhs, 3)
        elseif oper === (sin)
            lhs = args[1]
            rhs = Symbolics.term(asin, rhs)
        elseif oper === (cos)
            lhs = args[1]
            rhs = Symbolics.term(acos, rhs)
        elseif oper === (tan)
            lhs = args[1]
            rhs = Symbolics.term(atan, rhs)
        elseif oper === (asin)
            lhs = args[1]
            rhs = Symbolics.term(sin, rhs)
        elseif oper === (acos)
            lhs = args[1]
            rhs = Symbolics.term(cos, rhs)
        elseif oper === (atan)
            lhs = args[1]
            rhs = Symbolics.term(tan, rhs)
        end

        lhs = simplify(lhs)
    end
    return Symbolics.wrap(rhs)
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
    push!(r_addexpon, @acrule (~b)^(~f::(contains_var)) + (~d)^(~g::(contains_var)) => ~f*Symbolics.term(log, ~b) - ~g*Symbolics.term(log, ~d))
    push!(r_addexpon, @acrule (~a)*(~b)^(~f::(contains_var)) + (~d)^(~g::(contains_var)) => Symbolics.term(log, -~a) + ~f*Symbolics.term(log, ~b) - ~g*Symbolics.term(log, ~d))
    push!(r_addexpon, @acrule (~a)*(~b)^(~f::(contains_var)) + (~c)*(~d)^(~g::(contains_var)) => Symbolics.term(log, -(~a)/(~c)) + ~f*Symbolics.term(log, ~b) - ~g*Symbolics.term(log, ~d))

    for r in r_addexpon
        lhs = SymbolicUtils.Fixpoint(r)(Symbolics.unwrap(lhs))
    end

    return expand(lhs)
end


function attract_collect(lhs, var)
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

function nl_solve(lhs::Num, var)
    nx = n_func_occ(lhs, var)
    if nx == 0
        throw("Var not present in given expression.")
    elseif nx == 1
        return isolate(lhs, var)
    elseif nx > 1
        return isolate(attract_collect(lhs, var), var)
    end

end

function n_func_occ(expr, var)
    expr = Symbolics.unwrap(expr)
    !iscall(expr) && return n_occurrences(expr, var)
    args, cur_oper = Symbolics.arguments(expr), Symbolics.operation(expr)
    counted_ops = [sin, log, log2, log10, cos, tan, asin, acos, atan]
    n = 0


    if cur_oper === (*) || cur_oper === (+)

        outside = false
        for arg in args
            n_occurrences(arg, var) == 0 && continue

            if !iscall(arg) && any(isequal(var, x) for x in Symbolics.get_variables(arg)) && !outside
                n += 1
                outside = true
            end
            !iscall(arg) && continue
            oper = Symbolics.operation(arg)

            args_arg = Symbolics.arguments(arg)
            oper_arg = Symbolics.operation(arg)
            case_1_pow = oper_arg === (^) && n_occurrences(args_arg[2], var) == 0  && n_occurrences(args_arg[1], var) != 0
            case_2_pow = oper_arg === (^) && n_occurrences(args_arg[2], var) != 0  && n_occurrences(args_arg[1], var) == 0  

            if any(isequal(oper, op) for op in counted_ops)
                n += n_func_occ(args_arg[1], var)
            elseif case_2_pow
                n += n_func_occ(args_arg[2], var) 
            elseif check_poly_inunivar(arg, var) && !outside 
                n += 1
                outside = true
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
    !iscall(expr) && any(isequal(var, x) for x in Symbolics.get_variables(expr)) && return 1
    !iscall(expr) && return 0 

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

# @variables x
# nl_solve(2log(x+1) + log(x-1), x)

@variables x y 
# nl_solve(2^(x+1) + 5^(x+3), x)

n_func_occ(expand((1 + x)*Symbolics.term(log, 2) - (3 + x)*Symbolics.term(log, 5)), x)