using Symbolics
include("coeffs.jl")
include("nemo_stuff.jl")

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
            for arg in args
                vars = Symbolics.get_variables(arg)
                if any(isequal(x, var) for x in vars)
                    continue
                end
                lhs = lhs^(1/arg) 
                rhs = rhs^(1/arg)
            end

        elseif oper === (log)
            lhs = args[1]
            rhs = Base.MathConstants.e^rhs

        elseif oper === (log2)
            lhs = args[1]
            rhs = 2^rhs

        elseif oper === (log10)
            lhs = args[1]
            rhs = 10^rhs

        elseif oper === (sqrt)
            lhs = args[1]
            rhs = rhs^2

        elseif oper === (cbrt)
            lhs = args[1]
            rhs = rhs^3

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
    println(lhs)


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
#nl_solve(log(x+1)+log(x-1), x)
