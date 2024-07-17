using Symbolics, Test

"""
    detect_doubleangle1(arg, variable)

Detects if 2sin(x)cos(x) is present in any form in an expression
However, this function is currently useless since we can just 
apply the rewrite rule to the expression and see if it works or not.
For future reference, see if there are any positives of detecting the presence
of the identity first in the expression. Returns true if found,
false if not.

# Arguments
- arg: Should be SymbolicUtils.BasicSymbolic{Real} only.
- variable: Num and should satisfy Symbolics.is_singleton

# Examples
```jldoctest
julia> detect_doubleangle1(Symbolics.unwrap(2sin(x)cos(x)), x)
true
```
"""
function detect_doubleangle1(arg, var)
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

"""
    detect_trig(lhs, var)

Traverses lhs arg by arg (if oper === (+)) and
calls detect_(any trig identity detector) on it and
if any detection function returns true, the whole
function returns true. P.S. currently only
detect_doubleangle1 is implemented as the detect function.
Returns true if found, false if not.

# Arguments
- lhs: A Symbolics Num expression.
- variable: Num and should satisfy Symbolics.is_singleton

# Examples
```jldoctest
julia> detect_trig(2sin(x)cos(x), x)
true
```
"""
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


"""
    detect_addlog(lhs, var)

Traverses expression arg by arg and if two logs of
degree 1 are found it returns true if both of them
have n_func_occ of 1, so they're solvable after 
getting attracted. Returns true if found, 
false if not.

# Arguments
- lhs: A Symbolics Num expression.
- variable: Num and should satisfy Symbolics.is_singleton

# Examples
```jldoctest
julia> detect_addlogs(log(x+1) + log(x-1), x)
true
```
"""
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

"""
    detect_exponential(lhs, var)

Traverses an expression and attempts to match the
expression to ab^f(x) + cd^g(x). Hence why it
expects only 2 arguments in the input expression.
Returns true if found, false if not.

# Arguments
- lhs: A Symbolics Num expression.
- variable: Num and should satisfy Symbolics.is_singleton

# Examples
```jldoctest
julia> detect_exponential(2*3^(x+1) + 5^(x^2 + 3), x)
true
```
"""
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

"""
    attract_logs(lhs, var)

Rewrites any 2 logs in lhs where the log contains the
input variable var. so log(x) + log(y) is not attracted.
When attracted, it returns an slog.

# Arguments
- lhs: A Symbolics Num expression.
- variable: Num and should satisfy Symbolics.is_singleton

# Examples
```jldoctest
julia> RootFinding.attract_logs(log(x+1)+log(x-1), x)
Main.RootFinding.slog(-1 + x^2)
```
"""
function attract_logs(lhs, var)
    contains_var(arg) = occursin(string(var), string(arg))

    r_addlogs = Vector{Any}() 
    push!(r_addlogs, @acrule log(~x::(contains_var)) + log(~y::(contains_var)) => slog(~x * ~y))
    push!(r_addlogs, @acrule ~z*log(~x::(contains_var)) + log(~y::(contains_var)) => slog((~x)^(~z) * ~y))
    push!(r_addlogs, @acrule ~z*log(~x::(contains_var)) + ~h*log(~y::(contains_var)) => slog((~x)^(~z) * (~y)^(~h)))


    lhs = expand(simplify(lhs, rewriter=SymbolicUtils.Postwalk(SymbolicUtils.Chain(r_addlogs))))

    return lhs
end

"""
    attract_exponential(lhs, var)

Rewrites ab^f(x) + cd^g(x) into 
f(x) * log(b) - g(x) * log(d) + log(-a/c)
Which is a solvable polynomial, hence reducing the smart
number of occurrences of x in the expression.

# Arguments
- lhs: A Symbolics Num expression.
- variable: Num and should satisfy Symbolics.is_singleton

# Examples
```jldoctest
julia> attract_exponential(2^(x+1) + 5^(x+3), x)
Main.RootFinding.slog(2) + log(complex(-1)) - 3Main.RootFinding.slog(5) + x*Main.RootFinding.slog(2) - x*Main.RootFinding.slog(5)
"""
function attract_exponential(lhs, var)
    lhs = Symbolics.unwrap(lhs)
    contains_var(arg) = occursin(string(var), string(arg))

    r_addexpon = Vector{Any}()
    push!(r_addexpon, @acrule (~b)^(~f::(contains_var)) + (~d)^(~g::(contains_var)) => ~f*Symbolics.term(slog, ~b) - ~g*Symbolics.term(slog, ~d) + Symbolics.term(log, Symbolics.term(complex, -1)))
    push!(r_addexpon, @acrule (~a)*(~b)^(~f::(contains_var)) + (~d)^(~g::(contains_var)) => ~f*Symbolics.term(slog, ~b) - ~g*Symbolics.term(slog, ~d) + Symbolics.term(slog, -~a))
    push!(r_addexpon, @acrule (~a)*(~b)^(~f::(contains_var)) + (~c)*(~d)^(~g::(contains_var)) => ~f*Symbolics.term(slog, ~b) - ~g*Symbolics.term(slog, ~d) + Symbolics.term(slog, -(~a)//(~c)))

    lhs = expand(simplify(lhs, rewriter=SymbolicUtils.Postwalk(SymbolicUtils.Chain(r_addexpon))))

    return expand(lhs)
end

"""
    attract_trig(lhs, var)

This function brute forces a bunch of trig identities
on the input lhs as an attempt to reduce the occurrences
of x in the lhs.

# Arguments
- lhs: A Symbolics Num expression.
- variable: Num and should satisfy Symbolics.is_singleton

# Examples
```jldoctest
julia> attract_trig(2sin(x)cos(x), x)
sin(2x)

julia> attract_trig(sin(x)^2 + cos(x)^2, x)
1

julia> RootFinding.attract_trig(cosh(x)^2 + sinh(x)^2, x)
cosh(2x)
"""
function attract_trig(lhs, var)
    lhs = Symbolics.unwrap(lhs)
    contains_var(arg) = occursin(string(var), string(arg))

    # r_doubleangle1 = @acrule 2*sin(~x::(contains_var))*cos(~x::(contains_var)) => sin(2*~x)
    r_trig = [
        @acrule(sin(~x::(contains_var))^2 + cos(~x::(contains_var))^2 => one(~x))
        @acrule(sin(~x::(contains_var))^2 + -1 => -1*cos(~x)^2)
        @acrule(cos(~x::(contains_var))^2 + -1 => -1*sin(~x)^2)

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

    lhs = expand(simplify(lhs, rewriter=SymbolicUtils.Postwalk(SymbolicUtils.Chain(r_trig))))

    return lhs
end

@variables x
lhs = 2sin(x+1)cos(x+1)
attract_trig(lhs, x) 
