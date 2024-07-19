using Symbolics


"""
    turn_to_poly(expr, var)

Traverses an expression and keeps track of transcendental functions,
if multiple types or versions of transcendental functions are found,
it aborts subbing into the expression since this expression becomes
and unsolvable by polynomial substitution. Say ``log(x)^2 + sin(x) + 1``
What should we sub here? no sub makes it a solvable polynomial.
Moreover, if the transcendental functions are not exactly the same,
it also aborts. E.g. ``log(x+1)^2 + log(x) + 1``. here we have two
"versions" of logs, so we can not sub here too.

# Arguments
- expr A Symbolics Num expression.
- var: Num and should satisfy Symbolics.is_singleton

# Examples
```jldoctest
julia> new_expr, sub = turn_to_poly(log(x)^2 + log(x) + 1, x)
(1 + var"##231" + var"##231"^2, Dict{Any, Any}(var"##231" => log(x)))

julia> new_expr
1 + var"##231" + var"##231"^2

julia> sub
Dict{Any, Any} with 1 entry:
  var"##231" => log(x)



julia> new_expr, sub = turn_to_poly(9^x + 3^x + 2, x)
(2 + var"##237" + var"##237"^2, Dict{Any, Any}(var"##237" => 3^x))

julia> new_expr
2 + var"##237" + var"##237"^2

julia> sub
Dict{Any, Any} with 1 entry:
  var"##237" => 3^x
```
"""
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

"""
    trav_pow(args, index, var, broken, sub)

Traverses an argument passed from ``turn_to_poly`` if it
satisfies ``oper === (^)``. Returns sub if changed from 0
to a new transcendental function or its value is 
kept the same, and false if these 2 cases do not occur.

# Arguments
- args: The original arguments array of the expression passed to ``turn_to_poly``
- index: The index of the argument in the args array that needs to be traversed
- var: Num and should satisfy Symbolics.is_singleton
- broken: A reference bool which the value of is altered if a different transc function is found than sub
- sub: The sub stored in memory in ``turn_to_poly``. 

# Examples
```jldoctest
julia> trav_pow([Symbolics.unwrap(9^x)], 1, x, Ref(false), 3^x)
3^x

julia> trav_pow([Symbolics.unwrap(x^2)], 1, x, Ref(false), 3^x)
false
```
"""
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


"""
    trav_mult(arg, var, broken, sub)

Traverses an argument passed from ``turn_to_poly`` if it
satisfies ``oper === (*)``. Returns sub whether its changed from 0
to a new transcendental function or its value is 
kept the same, but changes broken if these 2 cases do not occur. It
traverses the * argument by sub_arg and compares it to sub using
the function ``add_sub``

# Arguments
- arg: The arg to be traversed.
- var: Num and should satisfy Symbolics.is_singleton
- broken: A reference bool which the value of is altered if a different transc function is found than sub
- sub: The sub stored in memory in ``turn_to_poly``. 

# Examples
```jldoctest
julia> trav_mult(Symbolics.unwrap(9*log(x)), x, Ref(false), log(x))
log(x)

julia> trav_mult(Symbolics.unwrap(9*log(x)^2), x, Ref(false), log(x))
log(x)

# value of broken is changed here to true
julia> trav_mult(Symbolics.unwrap(9*log(x+1)), x, Ref(false), log(x))
log(x)
```
"""
function trav_mult(arg, var, broken, sub)
    args_arg = Symbolics.arguments(arg)
    for (i, arg2) in enumerate(args_arg)
        !Symbolics.iscall(arg2) && continue

        oper = Symbolics.operation(arg2)
        if oper === (^) 
            tp = trav_pow(args_arg, i, var, broken, sub)
            sub = isequal(tp, false) ? sub : tp 
            continue
        end

        isequal(add_sub(sub, arg2, var, broken), false) && continue
        sub = arg2
    end
    return sub
end


"""
    add_sub(sub, arg, var, broken::Ref{Bool})

Checks if arg contains var, and compares it directly to sub,
if they're the same, returns true, if sub is 0 and arg contains
a transcendental function in var, returns true. If both cases are
not met, false is returned indicating that the previous function
in the stack should not sub. ``add_sub`` also changes ``broken``
to true if arg contains a transcendental function thats different from
sub in any way.

# Arguments
- sub: The sub stored in memory in ``turn_to_poly``. 
- arg: The arg to be compared to sub.
- var: Num and should satisfy Symbolics.is_singleton
- broken: A reference bool which the value of is altered if a different transc function is found than sub

# Examples
```jldoctest
julia> add_sub(3^x, Symbolics.unwrap(3^x), x, Ref(false))
true

julia> add_sub(0, Symbolics.unwrap(log(x)), x, Ref(false))
true

# broken here changed to true
julia> add_sub(log(x+1), Symbolics.unwrap(log(x)), x, Ref(false))
false
```
"""
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


"""
    contains_transcendental(arg, var, n_occ=1)

Simple helper function to check if arg contains a julia
transcendental function and has number of occurrences of var = n_occ.

# Arguments
- arg: The arg to be checked. Should be of type ``SymbolicUtils.BasicSymbolic{Real}``.
- var: Num and should satisfy Symbolics.is_singleton
- n_occ: The number of occurrences that arg should contain. Default value is 1.

# Examples
```jldoctest
julia> contains_transcendental(Symbolics.unwrap(log(x)), x)
true

julia> contains_transcendental(Symbolics.unwrap(x), x)false

julia> contains_transcendental(Symbolics.unwrap(x^2), x)
false

julia> contains_transcendental(Symbolics.unwrap(1), x)
false

julia> contains_transcendental(Symbolics.unwrap(sin(x)), x)
true
```
"""
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
