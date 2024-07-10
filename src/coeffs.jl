using Symbolics 

function sub(sub_counter, subs, place_to_sub)
    sub_var = Symbolics.variables("c"*string(sub_counter))[1]
    subs[sub_var] = deepcopy(place_to_sub)
    place_to_sub = sub_var.val
    sub_counter += 1
    return sub_counter, place_to_sub
end

function clean_f(filtered_expr, var, subs)
    filtered_expr = simplify_fractions(simplify(real(filtered_expr)))
    unwrapped_f = Symbolics.unwrap(filtered_expr)
    !iscall(unwrapped_f) && return filtered_expr
    oper = operation(unwrapped_f)

    if oper === (/)
        args = arguments(unwrapped_f)
        if any(isequal(var, x) for x in Symbolics.get_variables(args[2]))
            return filtered_expr
        end
        filtered_expr = args[1]
        @info substitute(args[2], subs, fold=false) != 0
    end
    return filtered_expr
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
        return subs, sub_var.val
    end
end

function merge_filtered_exprs(subs1, expr1, subs2, expr2)
    sub_counter = length(subs1)+1
    og_subs2 = deepcopy(subs2)
    if length(subs1) == 0 
        return subs2, expr2
    end
    for (var, sub) in og_subs2
        sub_var = Symbolics.variables("n"*string(sub_counter))[1]
        expr2 = Symbolics.substitute(expr2, Dict(var=>sub_var), fold=false)
        pop!(subs2, var)
        subs2[sub_var] = sub
        sub_counter += 1
    end
    og_subs2 = deepcopy(subs2)
    for (var, sub) in og_subs2
        sub_var = Symbolics.variables("c"*string(var)[2:end])[1]
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
    if !isequal(vars, []) && isequal(vars[1], expr)
        return (Dict{Any, Any}(), expr)
    elseif isequal(vars, [])
        return filter_stuff(expr)
    end

    args = arguments(expr)
    if isequal(typeof(expr), Symbolics.ComplexTerm{Real})
        subs1, subs2 = Dict(), Dict()
        expr1, expr2 = 0, 0
        if !isequal(expr.re, false)
            subs1, expr1 = filter_poly(expr.re, var)
        end
        if !isequal(expr.im, false)
            subs2, expr2 = filter_poly(expr.im, var)
        end
        
        merged_subs, expr2 = merge_filtered_exprs(subs1, expr1, subs2, expr2)
        vars = union!(collect(keys(subs1)), collect(keys(subs2)),
        Symbolics.get_variables(expr1), Symbolics.get_variables(expr2))
        c = 0
        i_var = Symbolics.variables("I"*string(c))[1]
        while any(isequal(i_var, present_var) for present_var in vars)
            c += 1
            i_var = Symbolics.variables("I"*string(c))[1]
        end
        merged_subs[i_var] = im
        expr = expr1 + i_var*expr2
        return merged_subs, expr
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
            if isequal(arg, var) || isequal(vars[1], arg)
                continue
            end
        end
        
        oper = Symbolics.operation(arg)
        monomial = arguments(arg)
        if oper === (^)
            if any(arg -> isequal(arg, var), monomial) 
                continue
            end
            # filter(args[1]), filter[args[2]] and then merge
            subs1, monomial[1] = filter_poly(monomial[1], var)
            subs2, monomial[2] = filter_poly(monomial[2], var)

            new_subs, monomial[2] = merge_filtered_exprs(subs1, monomial[1], subs2, monomial[2])
            subs, args[i] = merge_filtered_exprs(subs, expr, new_subs, arg)
            continue
        end

        if oper === (*)
            subs_of_monom = Dict{Any, Any}()
            for (j, x) in enumerate(monomial)
                type_x = typeof(x)
                vars = Symbolics.get_variables(x)
                if (!isequal(vars, []) && isequal(vars[1], var))  || isequal(type_x, Int64) || isequal(type_x, Rational{Int64})
                    continue
                end
                # filter each arg and then merge
                new_subs, monomial[j] = filter_poly(monomial[j], var)
                subs_of_monom, monomial[j] = merge_filtered_exprs(subs_of_monom, arg, new_subs, monomial[j])
            end
            subs, args[i] = merge_filtered_exprs(subs, expr, subs_of_monom, arg)
        end

        if oper === (/) || oper === (+)
            for (j, x) in enumerate(monomial)
                new_subs, new_filtered = filter_poly(monomial[j], var)
                subs, monomial[j] = merge_filtered_exprs(subs, expr, new_subs, new_filtered)
                if typeof(monomial[j]) == Num
                    monomial[j] = monomial[j].val
                end
            end
        end
    end

    # reassemble expr to avoid variables remembering original values issue
    args = arguments(expr)
    oper = operation(expr)
    new_expr = clean_f(Symbolics.term(oper, args...), var, subs)

    return (subs, new_expr)
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

@variables m x a
# filter_poly(a^10 + x, x)
expr = (BigInt(1)//1) - (BigInt(12)//1)*m + (BigInt(4)//1)*(m^2)
filter_poly(expr, x)