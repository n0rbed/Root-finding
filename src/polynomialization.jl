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
            sub = isequal(trav_pow(arg, var, broken, sub), false) ? sub : trav_pow(arg, var, broken, sub)
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

function trav_pow(arg, var, broken, sub)
    args_arg = Symbolics.arguments(arg)
    isequal(add_sub(sub, args_arg[1], var, broken), false) && return false
    return args_arg[1]
end

function trav_mult(arg, var, broken, sub)
    args_arg = Symbolics.arguments(arg)
    for arg2 in args_arg 
        !Symbolics.iscall(arg2) && continue

        oper = Symbolics.operation(arg2)
        if oper === (^) 
            sub = isequal(trav_pow(arg2, var, broken, sub), false) ? sub : trav_pow(arg2, var, broken, sub)
            continue
        end

        isequal(add_sub(sub, arg2, var, broken), false) && continue
        sub = arg2
    end
    return sub
end

function add_sub(sub, arg, var, broken::Ref{Bool})
    !Symbolics.iscall(arg) && return false
    arg_oper = Symbolics.operation(arg)

    friends = [sin, log, log2, log10, cos, tan, asin, acos, atan, exp]
    if any(isequal(arg_oper, oper) for oper in friends) && n_occurrences(arg, var) > 0
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



@variables x
expr = log(x)^2 + 3log(x) + 1
turn_to_poly(expr, x)