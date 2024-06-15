using Symbolics 
function filter_poly(expr)
    expr = Symbolics.unwrap(expr)
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
            args[i] = var
            continue
        end

        # handle "x" as an argument
        if isequal(Symbolics.get_variables(arg)[1], arg)
            continue
        end
        
        monomial = unsorted_arguments(arg)
        for (j, x) in enumerate(monomial)
            if !isequal(Symbolics.get_variables(x), [])  || isequal(typeof(x), Int64)
                continue
            end
            var = Symbolics.variables("c"*string(i))[1]
            subs[var] = x
            monomial[j] = var
        end
    end
    return (subs, expr)
end
@variables x
filter_poly(x^2 - (Symbolics.term(sqrt, 8))^2)