using Symbolics, Groebner, SymbolicUtils
using Nemo

coeff = Symbolics.coeff

function get_roots_deg2(expression, x)
    # ax^2 + bx + c = 0
    coeffs, constant = polynomial_coeffs(expression, [x])

    a = rationalize(coeffs[x^2])
    b = rationalize(get(coeffs, x, 0))
    c = rationalize(get(coeffs, x^0, 0))

    root1 = simplify((-b + Symbolics.term(sqrt, Symbolics.term(complex, (b^2 - 4(a*c))))) / 2a)
    root2 = simplify((-b - Symbolics.term(sqrt, Symbolics.term(complex, (b^2 - 4(a*c))))) / 2a)

    return [root1, root2]
end

function get_roots_deg3(expression, x)
    coeffs, constant = polynomial_coeffs(expression, [x])

    a = rationalize(coeffs[x^3])
    b = rationalize(get(coeffs, x^2, 0))
    c = rationalize(get(coeffs, x, 0))
    d = rationalize(get(coeffs, x^0, 0))

    Q = ((3*a*c) - b^2)//(9a^2)
    R = (9*a*b*c - ((27*(a^2)*d)+2b^3))//(54a^3)

    S = 0
    T = 0

    # cbrt(negative real numbers) evaluates normally with symbolics
    # while (im)^(1//3) also evaluates normally,
    # the other cases cbrt(im) and (-real)^(1/3) do not evaluate
    if eval(Symbolics.toexpr((Q^3+R^2))) > 0 
        S = Symbolics.term(cbrt, (R + Symbolics.term(sqrt, (Q^3+R^2))))
        T = Symbolics.term(cbrt, (R - Symbolics.term(sqrt, (Q^3+R^2))))
    else
        S = (Symbolics.term(complex, (R + im*Symbolics.term(sqrt, -(Q^3+R^2)))))^(1//3)
        T = (Symbolics.term(complex, (R - im*Symbolics.term(sqrt, -(Q^3+R^2)))))^(1//3)
    end


    root1 = S + T - (b//(3*a))
    root2 = -((S+T)//2) - (b//(3*a)) + (im*(Symbolics.term(sqrt, 3))/2)*(S-T)
    root3 = -((S+T)//2) - (b//(3*a)) - (im*(Symbolics.term(sqrt, 3))/2)*(S-T)

    return [root1, root2, root3]
end


function get_roots_deg4(expression, x)
    coeffs, constant = polynomial_coeffs(expression, [x])

    a = rationalize(coeffs[x^4])
    b = rationalize(get(coeffs, x^3, 0))
    c = rationalize(get(coeffs, x^2, 0))
    d = rationalize(get(coeffs, x, 0))
    e = rationalize(get(coeffs, x^0, 0))

    p = (8(a*c)-3(b^2))//(8(a^2))

    q = (b^3 - 4(a*b*c) + 8(d*a^2))//(8*a^3)

    r = (-3(b^4) + 256(e*a^3) - 64(d*b*a^2) + 16(c*(b^2)*a))//(256*a^4)

    @variables m y
    eq_m = 8m^3 + 8(p)*m^2 + (2(p^2) - 8r)m - q^2
    roots_m = solve(eq_m, m)
    m = 0
    for root in roots_m
        if !isequal(root, 0)
            m = root
            break
        end
    end

    arr = get_yroots(m, p, q)
    for (i, root) in enumerate(arr)
        arr[i] = root - (b//(4a))
    end

    return arr
end

function get_yroots(m, p, q)
    a = 1
    b1 = Symbolics.term(sqrt, Symbolics.term(complex, 2m))
    c1 = (p//2) + m - (q//(2*Symbolics.term(sqrt, Symbolics.term(complex, 2m))))
    b2 = -Symbolics.term(sqrt, Symbolics.term(complex, 2m))
    c2 = (p//2) + m + (q//(2*Symbolics.term(sqrt, Symbolics.term(complex, 2m))))

    root1 = simplify((-b1 + Symbolics.term(sqrt, Symbolics.term(complex, (b1^2 - 4(a*c1))))) / 2a)
    root2 = simplify((-b1 - Symbolics.term(sqrt, Symbolics.term(complex, (b1^2 - 4(a*c1))))) / 2a)
    root3 = simplify((-b2 + Symbolics.term(sqrt, Symbolics.term(complex, (b2^2 - 4(a*c2))))) / 2a)
    root4 = simplify((-b2 - Symbolics.term(sqrt, Symbolics.term(complex, (b2^2 - 4(a*c2))))) / 2a)

    return [root1, root2, root3, root4]
end

function get_roots(expression, x)
    degree = Symbolics.degree(expression, x)

    if degree == 0 && expression == 0
        return 0
    elseif degree == 0 && expression != 0
        throw("Not a valid statement")
    end

    if degree == 1
        # mx + c = 0
        coeffs, constant = polynomial_coeffs(expression, [x])
        m = rationalize(get(coeffs, x, 0))
        c = rationalize(get(coeffs, x^0, 0))
        root = -c // m
        return root
    end

    if degree == 2
        return get_roots_deg2(expression, x)
    end

    if degree == 3
        return get_roots_deg3(expression, x)
    end

    if degree == 4
        return get_roots_deg4(expression, x)
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

    u, factors = factor_use_nemo(expression)

    arr_roots = []

    if degree < 5 && length(factors) == 1
        return get_roots(expression, x)
    end

    if length(factors) != 1
        @assert isequal(expand(expression - u*prod(factors)), 0)

        for factor in factors
            append!(arr_roots, solve(factor, x))
        end
    end


    if isequal(arr_roots, [])
        throw("This expression does not have an exact solution, use a numerical method instead.")
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
        throw("GCD found between the expressions is 1, solve() can not solve this system of equations exactly")
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


function solve(eqs::Vector{Num}, vars::Vector{Num})
    eqs = convert(Vector{Any}, Symbolics.groebner_basis(eqs))
    all_roots = Dict()

    # initialize the place of each var
    for var in vars
        all_roots[var] = [] 
    end


    # get roots for first var (z in this case)
    solved = false
    while !solved
        # first, solve any single variable equations
        for eq in eqs
            for var in vars
                present_vars = Symbolics.get_variables(eq)
                if contains(var, present_vars) && size(present_vars, 1) == 1 
                    var = Symbolics.get_variables(eq)[1]
                    append!(all_roots[var], solve(Symbolics.wrap(eq), var))
                end
            end
        end

        # second, Substitute the roots of the variables where found
        for (i, eq) in pairs(eqs)
            for var in vars
                present_vars = Symbolics.get_variables(eq)

                if contains(var, present_vars) && !isequal(all_roots[var], [])
                    deleteat!(eqs, i)

                    for root in all_roots[var]
                        insert!(eqs, i, substitute(eq, Dict([var => root])))
                    end

                end
            end
        end


        solved = true
        for (var, value) in all_roots
            if isequal(value, [])
                solved = false
            end
        end
    end
    return all_roots
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
    sym_unit, sym_factors
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
    sym_gcd
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



#@variables x
#solve(x^4 + 1, x)