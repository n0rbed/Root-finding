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

function attract_logs(lhs, var)
    contains_var(arg) = occursin(string(var), string(arg))

    r_addlogs = Vector{Any}() 
    push!(r_addlogs, @acrule log(~x::(contains_var)) + log(~y::(contains_var)) => log(~x * ~y))
    push!(r_addlogs, @acrule ~z*log(~x::(contains_var)) + log(~y::(contains_var)) => log((~x)^(~z) * ~y))
    push!(r_addlogs, @acrule ~z*log(~x::(contains_var)) + ~h*log(~y::(contains_var)) => log((~x)^(~z) * (~y)^(~h)))


    for r in r_addlogs
        lhs = SymbolicUtils.Fixpoint(r)(Symbolics.unwrap(lhs))
    end

    return lhs
end


function attract_exponential(lhs, var)
    contains_var(arg) = occursin(string(var), string(arg))

    r_addexpon = Vector{Any}()
    push!(r_addexpon, @acrule (~b)^(~f::(contains_var)) + (~d)^(~g::(contains_var)) => ~f*Symbolics.term(log, ~b) - ~g*Symbolics.term(log, ~d) + Symbolics.term(log, Symbolics.term(complex, -1)))
    push!(r_addexpon, @acrule (~a)*(~b)^(~f::(contains_var)) + (~d)^(~g::(contains_var)) => Symbolics.term(log, -~a) + ~f*Symbolics.term(log, ~b) - ~g*Symbolics.term(log, ~d))
    push!(r_addexpon, @acrule (~a)*(~b)^(~f::(contains_var)) + (~c)*(~d)^(~g::(contains_var)) => Symbolics.term(log, -(~a)/(~c)) + ~f*Symbolics.term(log, ~b) - ~g*Symbolics.term(log, ~d))

    for r in r_addexpon
        lhs = SymbolicUtils.Fixpoint(r)(Symbolics.unwrap(lhs))
    end

    return expand(lhs)
end


