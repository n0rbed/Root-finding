using Symbolics
import Nemo


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
    return true
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
    nemo_factors = collect(keys(nemo_fac.fac)) 
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

