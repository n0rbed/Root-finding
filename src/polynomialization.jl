using Symbolics
include("coeffs.jl")
include("nemo_stuff.jl")
include("solve_helpers.jl")
include("postprocess.jl")
include("univar.jl")
include("isoa_helpers.jl")


function turn_to_poly(expr, var)
    expr = Symbolics.unwrap(expr)
    !Symbolics.iscall(expr) && return expr

    args = Symbolics.arguments(expr)

    sub = 0
    for (i, arg) in enumerate(args)
        !Symbolics.iscall(arg) && continue
        arg_oper = Symbolics.operation(arg)


        args_arg = Symbolics.arguments(arg)
        if arg_oper === (^)
            isequal(add_sub(sub, args_arg[1], var), false) && continue 
            sub = args_arg[1]
            continue
        end
        if arg_oper === (*)
            for arg2 in args_arg 
                isequal(add_sub(sub, arg2, var), false) && continue
                sub = arg2
                continue
            end
        end
        isequal(add_sub(sub, arg, var), false) && continue
        sub = arg
    end

    new_var = gensym()
    new_var = (@variables $new_var)[1]
    return ssubs(expr, Dict(sub=>new_var)), Dict{Any,Any}(new_var=>sub)
    
end

function add_sub(sub, arg, var)
    !Symbolics.iscall(arg) && return sub
    arg_oper = Symbolics.operation(arg)

    friends = [sin, log, log2, log10, cos, tan, asin, acos, atan, exp]
    if any(isequal(arg_oper, oper) for oper in friends) && n_occurrences(arg, var) > 0
        if isequal(sub, 0) 
            return true
        elseif !isequal(sub, arg)
            return false
        end
    end

    cond1 = n_occurrences(arg, var) > 0
    cond2 = !isequal(sub, arg)
    return (cond1 && cond2)
end



@variables x
expr = log(x)^2 + log(x) + 1
turn_to_poly(expr, x)