using Symbolics

function isolate(lhs::Num, var)
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
        end

        lhs = simplify(lhs)
        rhs = simplify(rhs)
    end
    return Symbolics.wrap(rhs)
end


function attract_collect(lhs, var)
    @variables x y
    r_addlogs = @rule log(~x) + log(~y) => log(~x * ~y)
    unwrapped_lhs = Symbolics.unwrap(lhs)
    args = unsorted_arguments(unwrapped_lhs)
    operation = Symbolics.operation(unwrapped_lhs)

    i = 1
    j = i + 1
    while j <= length(args)
        if Symbolics.operation(args[i]) === (log) && n_occurrences(args[i], var) != 0 && 
            Symbolics.operation(args[j]) === (log) && n_occurrences(args[j], var) != 0
            lhs = expand(simplify(lhs, RuleSet([r_addlogs])))
        end
        i += 1
        j += 1
    end
    if n_occurrences(lhs, var) == 1
        return lhs
    end

end
function nl_solve(lhs::Num, var)
    nx = n_occurrences(lhs, var)
    if nx == 0
        throw("Var not present in given expression." )
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
            if !iscall(arg) && any(isequal(var, x) for x in Symbolics.get_variables(arg)) && !outside
                n += 1
                outside = true
            end
            !iscall(arg) && continue
            oper = operation(arg)

            if any(isequal(oper, op) for op in counted_ops) 
                n += n_func_occ(arguments(arg)[1], var)
            elseif n_occurrences(arg, var) > 0 && !outside
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
        args = unsorted_arguments(argument)
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
