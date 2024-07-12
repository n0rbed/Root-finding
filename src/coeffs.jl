using Symbolics 

# TODO(Alex): it is not clear to me at the moment if we should use
# Symbolics.term or TermInterface.maketerm

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

    if expr isa Integer || expr isa Rational
        return Dict(), expr

    else
        expr = isequal(expr, true) ? 1 : expr
        subs = Dict{Any, Any}()

        expr = sub(subs, expr)
        return subs, expr 
    end
end

function filter_subexpr(expr, var)
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

        if !isequal(expr.re, 0)
            subs1, expr1 = filter_subexpr(expr.re, var)
        end
        if !isequal(expr.im, 0)
            subs2, expr2 = filter_subexpr(expr.im, var)
        end

        subs = merge(subs1, subs2)
        i_var = gensym()
        i_var = (@variables $i_var)[1]

        subs[i_var] = im
        expr = Symbolics.unwrap(expr1 + i_var*expr2)

        args = arguments(expr)
        oper = operation(expr)
        return subs, Symbolics.term(oper, args...)
    end

    subs = Dict{Any, Any}()
    for (i, arg) in enumerate(args)
        # handle constants
        vars = Symbolics.get_variables(arg)
        if isequal(vars, [])
            if arg isa Integer || arg isa Rational
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
            subs1, monomial[1] = filter_subexpr(monomial[1], var)
            subs2, monomial[2] = filter_subexpr(monomial[2], var)

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
                new_subs, monomial[j] = filter_subexpr(monomial[j], var)
                merge!(subs_of_monom, new_subs)
            end
            merge!(subs, subs_of_monom)
        end

        if oper === (/) || oper === (+)
            for (j, x) in enumerate(monomial)
                new_subs, new_filtered = filter_subexpr(monomial[j], var)
                merge!(subs, new_subs)
            end
        end
    end

    args = arguments(expr)
    oper = operation(expr)
    expr = Symbolics.term(oper, args...)
    return subs, expr
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

    subs, expr = filter_subexpr(expr, var)

    # reassemble expr to avoid variables remembering original values issue
    # CANT DO THIS FOR EVERY RECURSION ALONE, HAVE TO BE AT FINAL STAGE OF FILTER_POLY ONLY
    args = arguments(expr)
    oper = operation(expr)
    new_expr = clean_f(Symbolics.term(oper, args...), var, subs)

    return subs, new_expr
end


function sdegree(coeffs, var)
    degree = 0
    vars = collect(keys(coeffs))
    for n in vars
        isequal(n, 1) && continue
        isequal(n, var) && degree > 1 && continue

        if isequal(n, var) && degree < 1 
            degree = 1
            continue
        end

        args = arguments(n)
        if args[2] > degree 
            degree = args[2]
        end
    end
    return degree
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

