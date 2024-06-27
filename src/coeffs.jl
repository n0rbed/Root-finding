using Symbolics 

function sub(sub_counter, subs, place_to_sub)
    sub_var = Symbolics.variables("c"*string(sub_counter))[1]
    subs[sub_var] = deepcopy(place_to_sub)
    place_to_sub = sub_var.val
    sub_counter += 1
    return sub_counter, place_to_sub
end


function filter_stuff(expr)
    type_expr = typeof(expr)
    if isequal(type_expr, Int64) || isequal(type_expr, Rational{Int64})
        return Dict(), expr
    else
        # to do
        sub_var = Symbolics.variables("c1")[1]
        if isequal(expr, true)
            expr = 1
        end
        subs = Dict{Any, Any}(sub_var=>expr)
        return subs, Symbolics.wrap(sub_var)
    end
end

function merge_filtered_exprs(subs1, expr1, subs2, expr2)
    sub_counter = length(subs1)+1
    og_subs2 = deepcopy(subs2)
    for (var, sub) in og_subs2
        sub_var = Symbolics.variables("c"*string(sub_counter))[1]
        expr2 = Symbolics.substitute(expr2, Dict(var=>sub_var), fold=false)
        pop!(subs2, var)
        subs2[sub_var] = sub
    end
    merged_subs = Dict(subs1..., subs2...)
    return merged_subs, expr2
end

function filter_poly(og_expr, var)
    expr = deepcopy(og_expr)
    expr = Symbolics.unwrap(expr)
    vars = Symbolics.get_variables(expr)
    if !isequal(vars, []) && isequal(Symbolics.get_variables(expr)[1], expr)
        return (Dict(), Symbolics.wrap(expr))
    end
    if isequal(vars, [])
        return filter_stuff(expr)
    end

    args = unsorted_arguments(expr)
    if isequal(typeof(expr), Symbolics.ComplexTerm{Real})
        subs1, subs2 = Dict(), Dict()
        expr1, expr2 = 0, 0
        if !isequal(expr.re, false)
            subs1, expr1 = filter_poly(expr.re, x)
        end
        if !isequal(expr.im, false)
            subs2, expr2 = filter_poly(expr.im, x)
        end
        
        merged_subs, expr2 = merge_filtered_exprs(subs1, expr1, subs2, expr2)
        @variables I
        merged_subs[I] = im
        expr = expr1 + I*expr2
        return merged_subs, Symbolics.wrap(expr)
    end
    subs = Dict{Any, Any}()
    sub_counter = 1
    for (i, arg) in enumerate(args)
        # handle constants
        vars = Symbolics.get_variables(arg)
        type_arg = typeof(arg)
        if isequal(vars, [])
            if type_arg == Int64 || type_arg == Rational{Int64}
                continue
            end
            sub_counter, args[i] = sub(sub_counter, subs, args[i])
            continue
        end

        # handle "x" as an argument
        if length(vars) == 1
            if isequal(arg, var)
                continue
            elseif isequal(vars[1], arg)
                sub_counter, args[i] = sub(sub_counter, subs, args[i])
                continue
            end
        end
        
        oper = Symbolics.operation(arg)
        if oper === (^)
            monomial = unsorted_arguments(arg)
            if any(arg -> isequal(arg, var), monomial) 
                continue
            end
            sub_counter, args[i] = sub(sub_counter, subs, args[i])
            continue
        end

        monomial = unsorted_arguments(args[i])
        for (j, x) in enumerate(monomial)
            type_x = typeof(x)
            vars = Symbolics.get_variables(x)
            if (!isequal(vars, []) && isequal(vars[1], var))  || isequal(type_x, Int64) || isequal(type_x, Rational{Int64})
                continue
            end
            sub_counter, monomial[j] = sub(sub_counter, subs, monomial[j])
        end
    end
    return (subs, Symbolics.wrap(expr))
end


function lead_term(expr, var)
    subs, expr = filter_poly(expr, var)
    coeffs, constant = polynomial_coeffs(expr, [var])
    degree = Symbolics.degree(expr, var)
    lead_term = coeffs[var^degree]*var^degree
    for (var, sub) in subs
        lead_term = Symbolics.substitute(lead_term, Dict([var => sub]), fold=false)
    end

    return lead_term
end

function lead_coeff(expr, var)
    degree = Symbolics.degree(expr, var)
    lead_coeff = lead_term(expr, var) / (var^degree)
    return lead_coeff
end


@variables x
filter_poly(im*(x+3), x)