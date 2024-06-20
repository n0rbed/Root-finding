using Symbolics
function nl_solve(lhs::Num, var)
    rhs = 0
    lhs = Symbolics.unwrap(lhs)
    while !isequal(lhs, var)
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
        end
        if oper === (*)
            for arg in args
                vars = Symbolics.get_variables(arg)
                if any(isequal(x, var) for x in vars)
                    continue
                end
                lhs = lhs / arg
                rhs = rhs / arg
            end
        end
        if oper === (^)
            for arg in args
                vars = Symbolics.get_variables(arg)
                if any(isequal(x, var) for x in vars)
                    continue
                end
                lhs = lhs^(1/arg) 
                rhs = rhs^(1/arg)
            end
        end
        if oper === (log)
            lhs = args[1]
            rhs = Base.MathConstants.e^rhs
        end
        if oper === (sqrt)
            lhs = args[1]
            rhs = rhs^2
        end
        if oper === (cbrt)
            lhs = args[1]
            rhs = rhs^3
        end
        lhs = simplify(lhs)
        rhs = simplify(rhs)
    end
    return Symbolics.wrap(rhs)
end

function n_occurrences(expr, var)
    n = 0
    args = unsorted_arguments(Symbolics.unwrap(expr))

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
@variables x
n_occurrences(sin(x+log(x+1)) + 2x, x)