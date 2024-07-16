using Symbolics


# tries to find polynomial substituions for transcendental functions
# e.g. 1 + log(x) + log(x)^2 => 1 + X + X^2
function turn_to_poly(expr, var)
    expr = Symbolics.unwrap(expr)
    !Symbolics.iscall(expr) && return expr

    args = Symbolics.arguments(expr)

    sub = 0
    broken = Ref(false)

    for (i, arg) in enumerate(args)
        !Symbolics.iscall(arg) && continue
        arg_oper = Symbolics.operation(arg)


        if arg_oper === (^)
            tp = trav_pow(args, i, var, broken, sub)
            sub = isequal(tp, false) ? sub : tp
            continue
        end
        if arg_oper === (*)
            sub = trav_mult(arg, var, broken, sub)
            continue
        end
        isequal(add_sub(sub, arg, var, broken), false) && continue
        sub = arg
    end

    if broken[] || isequal(sub, 0) 
        return (expr, Dict{Any, Any}())
    end

    new_var = gensym()
    new_var = (@variables $new_var)[1]
    return ssubs(expr, Dict(sub=>new_var)), Dict{Any,Any}(new_var=>sub)
    
end

function trav_pow(args, index, var, broken, sub)
    args_arg = Symbolics.arguments(args[index])
    base = args_arg[1]
    power = args_arg[2]
    
    # case 1: log(x)^2 .... 9^x = 3^2^x = 3^2x = (3^x)^2
    !isequal(add_sub(sub, base, var, broken), false) && power isa Integer && return base

    # case 2: int^f(x)
    # n_func_occ may not be strictly 1, we could attempt attracting it after solving
    if base isa Integer && n_func_occ(power, var) == 1
        factors = prime_factors(base)
        length(factors) != 1 && return false
        b, p = factors[1] 
        new_b = b^power
        sub = isequal(sub, 0) ? new_b : sub
        if !isequal(sub, new_b)
            broken[] = true
            return false
        end
        new_b = Symbolics.term(^, new_b, p)
        args[index] = new_b
        return sub
    end

    return false
end

function trav_mult(arg, var, broken, sub)
    args_arg = Symbolics.arguments(arg)
    for (i, arg2) in enumerate(args_arg)
        !Symbolics.iscall(arg2) && continue

        oper = Symbolics.operation(arg2)
        if oper === (^) 
            sub = isequal(trav_pow(args_arg, i, var, broken, sub), false) ? sub : trav_pow(args_arg, i, var, broken, sub)
            continue
        end

        isequal(add_sub(sub, arg2, var, broken), false) && continue
        sub = arg2
    end
    return sub
end

function add_sub(sub, arg, var, broken::Ref{Bool})
    if contains_transcendental(arg, var)
        if isequal(sub, 0) 
            return true
        elseif !isequal(sub, arg)
            broken[] = true
            return false
        end
    end

    cond1 = n_occurrences(arg, var) > 0
    cond2 = isequal(sub, arg)
    return (cond1 && cond2)
end

function contains_transcendental(arg, var, n_occ=1)
    !Symbolics.iscall(arg) && return false
    arg_oper = Symbolics.operation(arg)

    friends = [sin, log, log2, log10, cos, tan, asin, acos, atan, exp]
    if any(isequal(arg_oper, oper) for oper in friends) && n_func_occ(arg, var) == n_occ
        return true
    end
    return false
end

# prime_factors(9) = [(3, 2)] i.e. 3^2 
# prime_factors(36) = [(2, 2), (3,2)] i.e. 2^2 + 3^2
# a.i. generated
function prime_factors(n::Integer)
    factors = []
    d = 2
    while n > 1
        count = 0
        while n % d == 0
            n ÷= d
            count += 1
        end
        if count > 0
            push!(factors, (d, count))
        end
        d += 1
        if d * d > n
            if n > 1
                push!(factors, (n, 1))
                break
            end
        end
    end
    return factors
end
