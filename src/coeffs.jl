using Symbolics 

function sub(subs, place_to_sub)
    sub_var = gensym()
    sub_var = (@variables $sub_var)[1]

    subs[sub_var] = deepcopy(place_to_sub)
    place_to_sub = sub_var.val

    return place_to_sub
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
        # TODO:  
        expr = isequal(expr, true) ? 1 : expr
        subs = Dict{Any, Any}()

        expr = sub(subs, expr)
        return subs, expr 
    end
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

        subs = merge(subs1, subs2)
        i_var = gensym()
        i_var = (@variables $i_var)[1]

        subs[i_var] = im
        expr = expr1 + i_var*expr2

        return subs, expr
    end

    subs = Dict{Any, Any}()
    for (i, arg) in enumerate(args)
        # handle constants
        vars = Symbolics.get_variables(arg)
        type_arg = typeof(arg)
        if isequal(vars, [])
            if type_arg == Int64 || type_arg == Rational{Int64}
                continue
            end
            args[i] = sub(subs, args[i])
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

            merge!(subs, merge(subs1, subs2))
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
                merge!(subs_of_monom, new_subs)
            end
            merge!(subs, subs_of_monom)
        end

        if oper === (/) || oper === (+)
            for (j, x) in enumerate(monomial)
                new_subs, new_filtered = filter_poly(monomial[j], var)
                merge!(subs, new_subs)
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
