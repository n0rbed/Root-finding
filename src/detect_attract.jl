# ADD LOGS #

function check_poly_inunivar(poly, var)
    subs, filtered = filter_poly(poly, var)
    coeffs, constant = polynomial_coeffs(filtered, [var])
    return isequal(constant, 0)
end

function arg_contains_log(arg, var)
    oper = Symbolics.operation(arg)
    isequal(oper, log) && return true

    if oper == (*)
        return find_logandexpon(arg, var, log, 1)
    end

    return false
end

function detect_addlogs(lhs, var)
    args = Symbolics.arguments(Symbolics.unwrap(lhs))
    oper = Symbolics.operation(Symbolics.unwrap(lhs))
    !isequal(oper, (+)) && return false
    
    found = [false, false]
    c = 1
    for arg in args
        !iscall(arg) && continue
        isequal(n_occurrences(arg, var), 0) && continue

        # arg has target var and iscall
        if arg_contains_log(arg, var) && c < 3
            found[c] = true
            c += 1
        end
    end

    return all(found)
end

function find_logandexpon(arg, var, oper, poly_index)
    args_arg = Symbolics.arguments(arg)

    oper_term, constant_term = 0, 0 

    for a in args_arg
        if n_occurrences(a, var) != 0 && Symbolics.iscall(a) && Symbolics.operation(a) == (oper) &&
                check_poly_inunivar(arguments(a)[poly_index], var)
            oper_term = a
        elseif n_occurrences(a, var) == 0
            constant_term = a
        end
    end

    !isequal(oper_term, 0) && !isequal(constant_term, 0) && return true
    return false
end

function detect_exponential(lhs, var)
    args = Symbolics.arguments(Symbolics.unwrap(lhs))
    oper = Symbolics.operation(Symbolics.unwrap(lhs))
    !isequal(oper, (+)) && return false
    length(args) != 2 && return false

    found = [false, false]
    c = 1
    for arg in args
        n_occurrences(arg, var) == 0 && continue
        !iscall(arg) && continue

        oper = Symbolics.operation(arg)
        if isequal(oper, ^)
            args_arg = Symbolics.arguments(arg)
            n_occurrences(args_arg[1], var) != 0 && continue
            !check_poly_inunivar(args_arg[2], var) && continue

            found[c] = true
            c += 1

        elseif isequal(oper, *) && find_logandexpon(arg, var, ^, 2)
            found[c] = true
            c += 1
        end
    end

    return all(found) 
end
