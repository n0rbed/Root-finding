using Symbolics, Groebner, SymbolicUtils

"""
    solve(expr, x, multiplicities=false)
The base solver chooses the appropriate solver after ensuring that the input is valid (i.e., it does some assertions).
All the examples in the REPL in the documentation of the other solvers can be repeated using RootFinding.solve too.
This makes it more convenient For future reference, we intend to add nl solve in this solver too. We can do this by
first checking if the input is a valid polynomial, if not, attempt nl solve.

# Available solvers
- `solve_univar`
- `solve_multivar`
- `solve_multipoly`
- TODO: `nl_solve`

# Arguments
- expr: Could be a single univar Symbolics expression in the form of a poly 
or multiple univar expressions or multiple multivar polys or a transcendental nonlinear function
which is solved by nl_solve.

- x: Could be a single variable or an array of variables which should be solved

- multiplicities: Should the output be printed `n` times where `n` is the number of occurrence of the root?
Say we have `(x+1)^2`, we then have 2 roots `x = -1`, by default the output is `[-1]`, If multiplicites is inputed as true,
then the output is `[-1, -1]`.

# Examples

## Example 1
```jldoctest
julia> solve(x^7 - 1, x)
2-element Vector{Any}:
roots_of((1//1) + x + x^2 + x^3 + x^4 + x^5 + x^6)
 1//1
```

## Example 2
```jldoctest
julia> expr = expand((x + 3)*(x^2 + 2x + 1)*(x + 2))
6 + 17x + 17(x^2) + 7(x^3) + x^4

julia> solve(expr, x)
3-element Vector{Any}:
 -2//1
 -1//1
 -3//1

julia> solve(expr, x, true)
4-element Vector{Any}:
 -2//1
 -1//1
 -1//1
 -3//1
```

## Example 3
```jldoctest
julia> eqs = [x+y^2+z, z*x*y, z+3x+y]
3-element Vector{Num}:
 x + z + y^2
       x*y*z
  3x + y + z

julia> solve(eqs, [x,y,z])
3-element Vector{Any}:
 Dict{Num, Any}(z => 0//1, y => 0//1, x => 0//1)
 Dict{Num, Any}(z => 0//1, y => 1//3, x => -1//9)
 Dict{Num, Any}(z => -1//1, y => 1//1, x => 0//1)
```


## Example 4
```jldoctest
julia> solve([x-1, x^3 - 1, x^2 - 1, (x-1)^20], x)
1-element Vector{Rational{BigInt}}:
 1
```
"""
function solve(expr, x, multiplicities=false)
    type_x = typeof(x)
    expr_univar = false
    x_univar = false


    if (type_x == Num || type_x == SymbolicUtils.BasicSymbolic{Real})
        x_univar = true
        @assert Symbolics.is_singleton(Symbolics.unwrap(x)) "Expected a variable, got $x"
    else
        for var in x
            @assert Symbolics.is_singleton(Symbolics.unwrap(var)) "Expected a variable, got $x"
        end
    end

    if !(expr isa Vector) 
        expr_univar = true
        check_expr_validity(expr)
    else
        for e in expr
            check_expr_validity(e)
        end
    end


    if x_univar

        sols = []
        if expr_univar
            sols = check_poly_inunivar(expr, x) ? solve_univar(expr, x, multiplicities) : nl_solve(expr, x)
        else
            exprs_ispoly = []
            for e in expr
                push!(e, check_poly_inunivar(e, x))
            end

            sols = all(exprs_ispoly) ? solve_multipoly(expr, x, multiplicities) : throw("Solve can not solve this input currently")
        end

        sols = map(postprocess_root, sols)
        return sols
    end

    if !expr_univar && !x_univar
        for e in expr
            for var in x
                @assert check_poly_inunivar(e, var) "This system can not be currently solved by solve."
            end
        end

        sols = solve_multivar(expr, x, multiplicities)
        for sol in sols
            for var in x
                sol[var] = postprocess_root(sol[var])
            end
        end

        return sols 
    end
end

function solve_univar(expression, x, mult=false)
    args = []
    mult_n = 1
    expression = Symbolics.unwrap(expression)
    expression = expression isa PolyForm ? SymbolicUtils.toterm(expression) : expression

    # handle multiplicities, i.e. (x+1)^20
    if Symbolics.iscall(expression)
        expr = Symbolics.simplify(deepcopy(expression))
        args = arguments(expr)
        operation = SymbolicUtils.operation(expr)
        if isequal(operation, ^) && args[2] isa Int64
            expression = Symbolics.wrap(args[1])
            mult_n = args[2]
        end
    end

    subs, filtered_expr = filter_poly(expression, x)
    coeffs, constant = polynomial_coeffs(filtered_expr, [x])

    degree = sdegree(coeffs, x)
    u, factors = factor_use_nemo(Symbolics.wrap(filtered_expr))
    factors = convert(Vector{Any}, factors)

    factors_subbed = map(factor -> ssubs(factor, subs), factors)
    arr_roots = []

    if degree < 5 && length(factors) == 1
        arr_roots = get_roots(expression, x)

        # multiplicities
        if mult
            og_arr_roots = deepcopy(arr_roots)
            for i = 1:(mult_n-1)
                append!(arr_roots, og_arr_roots)    
            end
        end

        return arr_roots
    end

    if length(factors) != 1
        @assert isequal(expand(filtered_expr - u*expand(prod(factors))), 0)

        for factor in factors_subbed
            roots = solve_univar(factor, x, mult)
            if isequal(typeof(roots), RootsOf)
                push!(arr_roots, roots)
            else
                append!(arr_roots, roots)
            end
        end
    end


    if isequal(arr_roots, [])
        return RootsOf(Symbolics.wrap(expression), Symbolics.wrap(x))
    end

    return arr_roots
end

# You can compute the GCD between a system of polynomials by doing the following:
# Get the GCD between the first two polys,
# and get the GCD between this result and the following index,
# say: solve([x^2 - 1, x - 1, (x-1)^20], x)
# the GCD between the first two terms is obviously x-1,
# now we call gcd_use_nemo() on this term, and the following,
# gcd_use_nemo(x - 1, (x-1)^20), which is again x-1.
# now we just need to solve(x-1, x) to get the common root in this
# system of equations.
function solve_multipoly(polys::Vector, x::Num, mult=false)
    polys = unique(polys)

    if length(polys) < 1
        throw("No expressions entered")
    end
    if length(polys) == 1
        return solve_univar(polys[1], x, mult)
    end

    gcd = gcd_use_nemo(polys[1], polys[2])

    for i = 3:length(polys)
        gcd = gcd_use_nemo(gcd, polys[i])
    end
    
    if isequal(gcd, 1)
        @info "Nemo gcd is 1."
        return []
    end

    return solve_univar(gcd, x, mult)
end


function solve_multivar(eqs::Vector{Num}, vars::Vector{Num}, mult=false)
    # do the trick
    @variables HAT
    old_len = length(vars)
    push!(vars, HAT)
    new_eqs = []
    generating = true
    while generating
        new_eqs = deepcopy(eqs)
        eq = HAT
        for i = 1:(old_len)
            eq -= rand(1:10)*vars[i]
        end
        push!(new_eqs, eq)
        new_eqs = convert(Vector{Any}, Symbolics.groebner_basis(new_eqs, ordering=Lex(vars)))

        if length(new_eqs) <= length(vars) 
            generating &= false
        end

        for i = 2:length(new_eqs)
            generating |= all(Symbolics.degree(var) > 1 for var in Symbolics.get_variables(new_eqs[i]))
        end
    end

    solutions = []

    # handle "unsolvable" cases
    if isequal(1, new_eqs[1])
        return solutions
    end
    if length(new_eqs) < length(vars)
        throw("Infinite number of solutions")
    end


    # first, solve any single variable equations
    i = 1
    while !(i > length(new_eqs))
            present_vars = Symbolics.get_variables(new_eqs[i])
        for var in vars
            if size(present_vars, 1) == 1 && isequal(var, present_vars[1])
                new_sols = solve(Symbolics.wrap(new_eqs[i]), var, mult)

                if length(solutions) == 0
                    append!(solutions, [Dict{Num, Any}(var => sol) for sol in new_sols])
                else
                    solutions = add_sol_to_all(solutions, new_sols, var)
                end

                deleteat!(new_eqs, i)
                i = i - 1
                break
            end
        end
        i = i + 1
    end


    # second, iterate over eqs and sub each found solution
    # then add the roots of the remaining unknown variables 
    for eq in new_eqs
        solved = false
        present_vars = Symbolics.get_variables(eq)
        size_of_sub = length(solutions[1])

        if size(present_vars, 1) <= (size_of_sub + 1)
            while !solved 
                subbed_eq = eq
                for (var, root) in solutions[1]
                    subbed_eq = Symbolics.substitute(subbed_eq, Dict([var => root]), fold=false)
                end

                var_tosolve = Symbolics.get_variables(subbed_eq)[1]
                new_var_sols = solve_univar(subbed_eq, var_tosolve, mult)
                add_sol!(solutions, new_var_sols, var_tosolve, 1)

                solved = all(x -> length(x) == size_of_sub+1, solutions)
            end
        end
    end

    pop!(vars)
    for roots in solutions
        delete!(roots, HAT)
    end

    return solutions
end

@variables a b c d e x
f = -19e + 93d*x + 7c*(x^2) + 45b*(x^3) - 10a*(x^4)
#get_roots_deg4(f, x)
