using Symbolics 

# TODO(Alex): it is not clear to me at the moment if we should use
# Symbolics.term or TermInterface.maketerm

"""
    sub(subs, place_to_sub)
Helper function for filter_poly. Generates a symbolics variable and adds it
to subs (changes the original array passed to the sub function). returns the
subbed value so the filter_poly function could change it in its scope.

# Arguments
- subs: Vector of dicts which consist of var_subbed => value_subbed in filter_poly
- place_to_sub: The place which should be substituted by a new variable.

"""
function sub(subs, place_to_sub)
    sub_var = gensym()
    sub_var = (@variables $sub_var)[1]

    subs[sub_var] = deepcopy(place_to_sub)
    place_to_sub = sub_var.val

    return place_to_sub
end

"""
    clean_f(filtered_expr, var, subs)

Helper function for `filter_poly` which is called directly before returning
the `filtered_expression`. This function aims to get the filtered expressions
resulting from `filter_poly` ready to get used by solve, `get_roots`, and any other
function in the library. An important feature of it is that it simplifies fractions
to make the output "valid" in the eyes of `factor_use_nemo` and other Nemo functions.

# Arguments
- `filtered_expr`: The output from `_filter_poly`.
- var: The variable which is filtered for.
- subs: Vector of dicts which consist of `var_subbed` => `value_subbed` in `_filter_poly`.

"""
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

"""
    filter_stuff(expr)
Helper function for _filter_poly and filter_poly which aims to filter
constants and things that dont contain variables like Symbolics.term(sqrt, 2).
Filters only if the passed expr is scary (i.e. not a Rational or Int).

# Arguments
- expr: The detected constant in _filter_poly

# Examples
```jldoctest
julia> filter_stuff(Symbolics.term(sqrt, 2))
(Dict{Any, Any}(var"##278" => sqrt(2)), var"##278")

julia> filter_stuff(123)
(Dict{Any, Any}(), 123)
```
"""
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

"""
    _filter_poly(expr, var)

Main mechanism for filter_poly. Traverses arguments as deep 
as needed/makes sense by recalling itself as many times as necessary.
Output is then returned to filter_poly. This function should not be used
(although functional) instead of filter_poly since it does not contain
the final touch of clean_f, which is needed if the output is going to get passed
to other functions.

# Arguments
- expr: Expression to be filtered passed from filter_poly
- var: Var that the expression should be filtered with respect to.

# Examples
```jldoctest
julia> _filter_poly(x + sqrt(2), x)
(Dict{Any, Any}(var"##239" => 1.4142135623730951), var"##239" + x)

julia> RootFinding._filter_poly(x*sqrt(2), x)
(Dict{Any, Any}(var"##240" => 1.4142135623730951), var"##240"*x)
```
"""
function _filter_poly(expr, var)
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
            subs1, expr1 = _filter_poly(expr.re, var)
        end
        if !isequal(expr.im, 0)
            subs2, expr2 = _filter_poly(expr.im, var)
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
        arg = Symbolics.unwrap(arg)
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
            if any(arg -> isequal(arg, var), monomial) continue
            end
            # filter(args[1]), filter[args[2]] and then merge
            subs1, monomial[1] = _filter_poly(monomial[1], var)
            subs2, monomial[2] = _filter_poly(monomial[2], var)

            merge!(subs, merge(subs1, subs2))
            continue
        end

        if oper === (*)
            subs_of_monom = Dict{Any, Any}()
            for (j, x) in enumerate(monomial)
                type_x = typeof(x)
                vars = Symbolics.get_variables(x)
                if (!isequal(vars, []) && isequal(vars[1], x))  || isequal(type_x, Int64) || isequal(type_x, Rational{Int64})
                    continue
                end
                # filter each arg and then merge
                new_subs, monomial[j] = _filter_poly(monomial[j], var)
                merge!(subs_of_monom, new_subs)
            end
            merge!(subs, subs_of_monom)
            continue
        end

        if oper === (/) || oper === (+)
            for (j, x) in enumerate(monomial)
                new_subs, new_filtered = _filter_poly(monomial[j], var)
                merge!(subs, new_subs)
            end
            continue
        end
    end

    args = arguments(expr)
    oper = operation(expr)
    expr = Symbolics.term(oper, args...)
    return subs, expr
end


"""
    filter_poly(og_expr, var)

Filters polynomial expressions from weird irrationals that mess up
the functionality of other helper functions for our solvers such as
polynomial_coeffs and factor_use_nemo.

# Arguments
- og_expr: Original passed expression. This is deepcopied so that the
user's expression is not altered unwillingly.
- var: Var that the expression should be filtered with respect to.

# Examples
```jldoctest
julia> filter_poly(x + 2im, x)
(Dict{Any, Any}(var"##244" => im), 2var"##244" + x)

julia> filter_poly((1/im)*x + 3*y*z, x)
(Dict{Any, Any}(var"##245" => -1.0, var"##246" => im), 3y*z + var"##245"*var"##246"*x)

julia> filter_poly((x+1)*Symbolics.term(log, 3), x)
(Dict{Any, Any}(var"##247" => log(3)), var"##247"*(1 + x))
```
"""
function filter_poly(og_expr, var)
    expr = deepcopy(og_expr)
    expr = Symbolics.unwrap(expr)
    vars = Symbolics.get_variables(expr)

    # handle edge cases
    if !isequal(vars, []) && isequal(vars[1], expr)
        return (Dict{Any, Any}(), expr)
    elseif isequal(vars, [])
        return filter_stuff(expr)
    end

    # core filter
    subs, expr = _filter_poly(expr, var)

    # reassemble expr to avoid variables remembering original values issue and clean
    args = arguments(expr)
    oper = operation(expr)
    new_expr = clean_f(Symbolics.term(oper, args...), var, subs)

    return subs, new_expr
end


"""
    sdegree(coeffs, var)
Gets the degree of a polynomial by traversing the `coeffs`
output from polynomial_coeffs.

# Arguments
- coeffs: output from polynomial_coeffs
- var: var present in polynomial_coeffs

# Examples
```jldoctest
julia> coeffs, constant = polynomial_coeffs(x^2 + x + 1, [x])
(Dict{Any, Any}(x^2 => 1, x => 1, 1 => 1), 0)

julia> sdegree(coeffs, x)
2


julia> coeffs, constant = polynomial_coeffs(x^12 + 3x^2, [x])
(Dict{Any, Any}(x^2 => 3, x^12 => 1), 0)

julia> sdegree(coeffs, x)
12
```
"""
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



########## NOT USED ########## 

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

