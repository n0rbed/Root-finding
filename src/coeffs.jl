using Symbolics 
function filter_poly(expr)
    expr = Symbolics.unwrap(expr)
    if isequal(Symbolics.get_variables(expr)[1], expr)
        return (Dict(), Symbolics.wrap(expr))
    end

    args = unsorted_arguments(expr)
    subs = Dict()
    for (i, arg) in enumerate(args)
        # handle constants
        vars = Symbolics.get_variables(arg)
        if isequal(vars, [])
            if typeof(arg) == Int64
                continue
            end
            var = Symbolics.variables("c"*string(i))[1]
            subs[var] = arg
            args[i] = var.val
            continue
        end

        # handle "x" as an argument
        if isequal(Symbolics.get_variables(arg)[1], arg)
            continue
        end
        
        oper = Symbolics.operation(arg)
        if oper === (^)
            continue
        end

        monomial = unsorted_arguments(arg)
        for (j, x) in enumerate(monomial)
            if !isequal(Symbolics.get_variables(x), [])  || isequal(typeof(x), Int64)
                continue
            end
            var = Symbolics.variables("c"*string(i))[1]
            subs[var] = x
            monomial[j] = var.val
        end
    end
    return (subs, Symbolics.wrap(expr))
end
