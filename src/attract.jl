using Symbolics, Test
function detect_doubleangle1(arg, var)
    # detect 2sin(x)cos(x) = sin(2x)
    !Symbolics.iscall(arg) && return false
    oper_arg = Symbolics.operation(arg)

    found = [false, false, false]

    c = 1
    sub_args = Symbolics.arguments(arg)
    oper_arg != (*) && return false
    length(sub_args) != 3 && return false

    sin_content = NaN
    cos_content = NaN
    if any(isequal(2, n) for n in sub_args)
        found[c] = true
        c += 1
    end
    for n in sub_args
        !Symbolics.iscall(n) && continue
        cur_content = Symbolics.arguments(n)[1]
        if (isequal(Symbolics.operation(n), sin) || isequal(Symbolics.operation(n), cos)) && n_occurrences(n, var) > 0
            Symbolics.operation(n) == sin ? sin_content = cur_content : cos_content = cur_content
            found[c] = true
            c += 1
        end
    end

    !isequal(sin_content, cos_content) && return false
    return all(found)
end

function detect_trig(lhs, var)
    u_lhs = Symbolics.unwrap(lhs)
    args = Symbolics.arguments(u_lhs)
    oper = Symbolics.operation(u_lhs)

    # should traverse entire expresssion in the future, so
    # u_lhs as a whole, then u_lhs arg by arg recursively,
    # say we have 2sin(x)cos(x) / y
    # or do we do this traversing outside? in the mother func attract?
    b = false
    if oper === (+)
        for arg in args
            b |= detect_doubleangle1(arg, var)
        end
    elseif oper === (*)
        b |= detect_doubleangle1(u_lhs, var)
    end
    return b
end


function detect_addlogs(lhs, var)
    u_lhs = Symbolics.unwrap(lhs)
    args = Symbolics.arguments(u_lhs)
    oper = Symbolics.operation(u_lhs)
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
    push!(r_addlogs, @acrule slog(~x::(contains_var)) + slog(~y::(contains_var)) => slog(~x * ~y))
    push!(r_addslogs, @acrule ~z*slog(~x::(contains_var)) + slog(~y::(contains_var)) => slog((~x)^(~z) * ~y))
    push!(r_addslogs, @acrule ~z*slog(~x::(contains_var)) + ~h*slog(~y::(contains_var)) => slog((~x)^(~z) * (~y)^(~h)))


    lhs = expand(simplify(lhs, SymbolicUtils.Postwalk(SymbolicUtils.Chain(r_addlogs))))

    return lhs
end


function attract_exponential(lhs, var)
    lhs = Symbolics.unwrap(lhs)
    contains_var(arg) = occursin(string(var), string(arg))

    r_addexpon = Vector{Any}()
    # we might wanna create our own slog function which handles weird input too
    push!(r_addexpon, @acrule (~b)^(~f::(contains_var)) + (~d)^(~g::(contains_var)) => ~f*Symbolics.term(slog, ~b) - ~g*Symbolics.term(slog, ~d) + Symbolics.term(log, Symbolics.term(complex, -1)))
    push!(r_addexpon, @acrule (~a)*(~b)^(~f::(contains_var)) + (~d)^(~g::(contains_var)) => ~f*Symbolics.term(slog, ~b) - ~g*Symbolics.term(slog, ~d) + Symbolics.term(slog, -~a))
    push!(r_addexpon, @acrule (~a)*(~b)^(~f::(contains_var)) + (~c)*(~d)^(~g::(contains_var)) => ~f*Symbolics.term(slog, ~b) - ~g*Symbolics.term(slog, ~d) + Symbolics.term(slog, -(~a)//(~c)))

    lhs = expand(simplify(lhs, SymbolicUtils.Postwalk(SymbolicUtils.Chain(r_addexpon))))

    return expand(lhs)
end

function attract_trig(lhs, var)
    lhs = Symbolics.unwrap(lhs)
    contains_var(arg) = occursin(string(var), string(arg))

    # r_doubleangle1 = @acrule 2*sin(~x::(contains_var))*cos(~x::(contains_var)) => sin(2*~x)
    r_trig = [
        @acrule(sin(~x::(contains_var))^2 + cos(~x::(contains_var))^2 => one(~x::(contains_var)))
        @acrule(sin(~x::(contains_var))^2 + -1        => -1*cos(~x::(contains_var))^2)
        @acrule(cos(~x::(contains_var))^2 + -1        => -1*sin(~x::(contains_var))^2)

        @acrule(cos(~x::(contains_var))^2 + -1*sin(~x::(contains_var))^2 => cos(2 * ~x))
        @acrule(sin(~x::(contains_var))^2 + -1*cos(~x::(contains_var))^2 => -cos(2 * ~x))
        @acrule(cos(~x::(contains_var)) * sin(~x::(contains_var)) => sin(2 * ~x)/2)

        @acrule(tan(~x::(contains_var))^2 + -1*sec(~x::(contains_var))^2 => one(~x))
        @acrule(-1*tan(~x::(contains_var))^2 + sec(~x::(contains_var))^2 => one(~x))
        @acrule(tan(~x::(contains_var))^2 +  1 => sec(~x)^2)
        @acrule(sec(~x::(contains_var))^2 + -1 => tan(~x)^2)

        @acrule(cot(~x::(contains_var))^2 + -1*csc(~x)^2 => one(~x))
        @acrule(cot(~x::(contains_var))^2 +  1 => csc(~x)^2)
        @acrule(csc(~x::(contains_var))^2 + -1 => cot(~x)^2)

        @acrule(cosh(~x::(contains_var))^2 + -1*sinh(~x)^2 => one(~x))
        @acrule(cosh(~x::(contains_var))^2 + -1  => sinh(~x)^2)
        @acrule(sinh(~x::(contains_var))^2 +  1  => cosh(~x)^2)

        @acrule(cosh(~x::(contains_var))^2 + sinh(~x::(contains_var))^2 => cosh(2 * ~x))
        @acrule(cosh(~x::(contains_var)) * sinh(~x::(contains_var)) => sinh(2 * ~x)/2)
    ]

    lhs = expand(simplify(lhs, SymbolicUtils.Postwalk(SymbolicUtils.Chain(r_trig))))

    return lhs
end

@variables x
lhs = 2sin(x+1)cos(x+1)
attract_trig(lhs, x) 
