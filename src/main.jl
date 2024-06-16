using Symbolics, Groebner, SymbolicUtils
using Nemo
include("univar.jl")
include("coeffs.jl")


function sub_roots(arr_roots, subs)
    for i = 1:length(arr_roots)
        vars = Symbolics.get_variables(arr_roots[i])
        for var in vars
#            try
                arr_roots[i] = substitute(arr_roots[i], Dict([var => subs[var]]))
#            catch e
#            end
        end
    end
end

function solve(expression, x)
    try
        if isequal(SymbolicUtils.operation(expression.val), ^) && SymbolicUtils.arguments(expression.val)[2] isa Int64
            expression = Symbolics.wrap(SymbolicUtils.arguments(expression.val)[1])
        end
    catch e
    end
    

    expression = expand(expression)
    expression = simplify.(expression)
    degree = Symbolics.degree(expression, x)

    subs, filtered_expression = filter_poly(expression)
    u, factors = factor_use_nemo(filtered_expression)

    arr_roots = []

    if degree < 5 && length(factors) == 1
        append!(arr_roots, get_roots(filtered_expression, x, subs))
        sub_roots(arr_roots, subs)
        return arr_roots
    end

    if length(factors) != 1
        @assert isequal(expand(filtered_expression - u*expand(prod(factors))), 0)

        for factor in factors
            append!(arr_roots, solve(factor, x))
        end
    end


    if isequal(arr_roots, [])
        throw("This expression does not have an exact solution, use a numerical method instead.")
    end

    sub_roots(arr_roots, subs)
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
function solve(polys::Vector, x::Num)
    polys = unique(polys)

    if length(polys) < 1
        throw("No expressions entered")
    end
    if length(polys) == 1
        return solve(polys[1], x)
    end

    gcd = gcd_use_nemo(polys[1], polys[2])

    for i = 3:length(polys)
        gcd = gcd_use_nemo(gcd, polys[i])
    end
    
    if isequal(gcd, 1)
        @info "Nemo gcd is 1."
        return []
    end
    return solve(gcd, x)
end

function contains(var, vars)
    for variable in vars
        if isequal(var, variable)
            return true
        end
    end
    return false
end

function add_sol(solutions, new_sols, var, index)
    sol_used = solutions[index]
    deleteat!(solutions, index)
    for new_sol in new_sols
        sol_used[var] = new_sol
        push!(solutions, deepcopy(sol_used))
    end
    return solutions
end

function add_sol_to_all(solutions, new_sols, var)
    existing_solutions = deepcopy(solutions)
    solutions = []
    for new_sol in new_sols
        copy_sol = deepcopy(existing_solutions)
        for i = 1:length(copy_sol)
            copy_sol[i][var] = new_sol
        end
        append!(solutions, copy_sol)
    end
    return solutions
end

function solve(eqs::Vector{Num}, vars::Vector{Num})
    eqs = convert(Vector{Any}, Symbolics.groebner_basis(eqs, ordering=Lex(vars)))

    solutions = []

    # handle "unsolvable" cases
    if isequal(1, eqs[1])
        return solutions
    end

    if length(eqs) < length(vars)
        throw("Infinite number of solutions")
    end

    # first, solve any single variable equations
    i = 1
    while !(i > length(eqs))
            present_vars = Symbolics.get_variables(eqs[i])
        for var in vars
            if size(present_vars, 1) == 1 && isequal(var, present_vars[1])
                new_sols = solve(Symbolics.wrap(eqs[i]), var)

                if length(solutions) == 0
                    append!(solutions, [Dict(var => sol) for sol in new_sols])
                else
                    solutions = add_sol_to_all(solutions, new_sols, var)
                end

                deleteat!(eqs, i)
                i = i - 1
                break
            end
        end
        i = i + 1
    end

    # second, iterate over eqs and sub each found solution
    # then add the roots of the remaining unknown variables 
    j = 1
    for eq in eqs
        solved = false
        present_vars = Symbolics.get_variables(eq)
        size_of_sub = length(solutions[1])

        if size(present_vars, 1) == (size_of_sub + 1)
            while !solved 
                subbed_eq = eq
                for (var, root) in solutions[1]
                    subbed_eq = substitute(subbed_eq, Dict([var => root]))
                end
                subbed_eq = Symbolics.wrap(subbed_eq)

                if isequal(subbed_eq, 0)
                    # this is my assumption for which an equation is redundant (most likely wrong)
                    # i think its redundant if the equation is divisible by something
                    # i.e. simplifiable, not sure how to describe this as code though. 
                    break
                end

                var_tosolve = Symbolics.get_variables(subbed_eq)[1]
                new_var_sols = solve(subbed_eq, var_tosolve)
                solutions = add_sol(solutions, new_var_sols, var_tosolve, 1)

                solved = all(x -> length(x) == j+1, solutions)
            end
        end
        if solved
            j = j + 1
        end
    end

    return solutions
end
    

# Map each variable of the given poly.
# Can be used to transform Nemo polynomial to Symbolics expression.
function nemo_crude_evaluate(poly::Nemo.MPolyRingElem, varmap)
    @assert Nemo.coefficient_ring(poly) in (Nemo.ZZ, Nemo.QQ)
    new_poly = 0
    for (i, term) in enumerate(Nemo.terms(poly))
        new_term = Rational(Nemo.coeff(poly, i))
        for var in Nemo.vars(term)
            exp = Nemo.degree(term, var)
            exp == 0 && continue
            new_var = varmap[var]
            new_term *= new_var^exp
        end
        new_poly += new_term
    end
    new_poly
end

# Checks that the expression is a polynomial with integer or rational
# coefficients
function check_polynomial(poly::Num)
    vars = Symbolics.get_variables(poly)
    distr, rem = Symbolics.polynomial_coeffs(poly, vars)
    @assert isequal(rem, 0) "Not a polynomial"
    @assert all(c -> c isa Integer || c isa Rational, collect(values(distr))) "Coefficients must be integer or rational"
end

# factor(x^2*y + b*x*y - a*x - a*b)  ->  (x*y - a)*(x + b)
function factor_use_nemo(poly::Num)
    check_polynomial(poly)
    Symbolics.degree(poly) == 0 && return poly, Num[]
    vars = Symbolics.get_variables(poly)
    nemo_ring, nemo_vars = Nemo.polynomial_ring(Nemo.QQ, map(string, vars))
    sym_to_nemo = Dict(vars .=> nemo_vars)
    nemo_to_sym = Dict(v => k for (k, v) in sym_to_nemo)
    nemo_poly = Symbolics.substitute(poly, sym_to_nemo)
    nemo_fac = Nemo.factor(nemo_poly)
    nemo_unit = Nemo.unit(nemo_fac)
    nemo_factors = collect(keys(nemo_fac.fac)) # TODO: do not forget multiplicities
    sym_unit = Rational(Nemo.coeff(nemo_unit, 1))
    sym_factors = map(f -> Symbolics.wrap(nemo_crude_evaluate(f, nemo_to_sym)), nemo_factors)

    for (i, fac) in enumerate(sym_factors)
        sym_factors[i] = fac^(collect(values(nemo_fac.fac))[i])
    end

    return sym_unit, sym_factors
end

# gcd(x^2 - y^2, x^3 - y^3) -> x - y
function gcd_use_nemo(poly1::Num, poly2::Num)
    check_polynomial(poly1)
    check_polynomial(poly2)
    vars1 = Symbolics.get_variables(poly1)
    vars2 = Symbolics.get_variables(poly2)
    vars = vcat(vars1, vars2)
    nemo_ring, nemo_vars = Nemo.polynomial_ring(Nemo.QQ, map(string, vars))
    sym_to_nemo = Dict(vars .=> nemo_vars)
    nemo_to_sym = Dict(v => k for (k, v) in sym_to_nemo)
    nemo_poly1 = Symbolics.substitute(poly1, sym_to_nemo)
    nemo_poly2 = Symbolics.substitute(poly2, sym_to_nemo)
    nemo_gcd = Nemo.gcd(nemo_poly1, nemo_poly2)
    sym_gcd = Symbolics.wrap(nemo_crude_evaluate(nemo_gcd, nemo_to_sym))
    return sym_gcd
end

# NOTE:
# We can use Gcd to solve systems of polynomial equations in 1 variable:
#   f_1(x) = ... = f_m(x) = 0.
#
# The algorithm goes as follows:
#   1. Compute g = gcd(f_1, ..., f_m).
#   2. Solve   g = 0.
#
# Example. Consider the system
#   f_1(x) = x^2 - 1,
#   f_2(x) = x^3 - 1.
#
#   1. Compute g = x - 1 = gcd(x^2 - 1, x^3 - 1).
#   2. Solve   g = 0  => x = 1.
#
# Therefore, 1 is the only solution. Indeed,
# - The roots of f_1(x) = 0 are 1, -1.
# - The roots of f_2(x) = 0 are 1, (-1 +- sqrt(3)*i)/2.
# - The solution of f_1 = f_2 = 0 is their common root: 1.
#@variables x y z
#eqs = [x-y-z, x+y-z^2, x^2 + y^2 - 1]
#solve(eqs, [x,y,z])
@variables x
exp = x^4 + 1
solve(exp, x)