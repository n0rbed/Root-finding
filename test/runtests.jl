using RootFinding, Symbolics
using Test

function sort_roots(roots)
    return sort(roots, lt = (x,y) -> real(x)==real(y) ? imag(x)<imag(y) : real(x)<real(y))
end

@variables x y z

@test isequal(solve(x+1, x), -1)

@test isequal(solve(2x+1, x), -1/2)

@test isequal(solve(x, x), 0) 

@test isequal(solve((x+1)^20, x), -1)

exp = x^2 + 1
arr_calcd_roots = sort_roots(eval.(Symbolics.toexpr.(solve(exp, x))))
arr_known_roots = sort_roots([-im, im])
@test isequal(arr_calcd_roots .≈ arr_known_roots, [1,1])

exp = x^2 + 2x + 10
arr_calcd_roots = sort_roots(eval.(Symbolics.toexpr.(solve(exp, x))))
arr_known_roots = sort_roots([-1 + 3im, -1 - 3im])
@test isequal(arr_calcd_roots .≈ arr_known_roots, [1,1])

exp = x^2 - 10x + 25
arr_calcd_roots = sort_roots(eval.(Symbolics.toexpr.(solve(exp, x))))
arr_known_roots = [5,5]
@test isequal(arr_calcd_roots .≈ arr_known_roots, [1,1])


exp = x^3 - 2x^2 + x - 2 
arr_calcd_roots = sort_roots(eval.(Symbolics.toexpr.(solve(exp, x))))
arr_known_roots = sort_roots([2, -im, im])
@test isequal(arr_calcd_roots .≈ arr_known_roots, [1,1,1])

exp = x^4 + 1
arr_calcd_roots = sort_roots(eval.(Symbolics.toexpr.(solve(exp, x))))
arr_known_roots = sort_roots(eval.([-(complex(-1))^(1/4),(complex(-1))^(1/4), (complex(-1))^(3/4), -(complex(-1))^(3/4)]))
@test isequal(arr_calcd_roots .≈ arr_known_roots, [1,1,1,1])


# Alex:

# solve returns evaluated sqrt(2)
# expr = x - Symbolics.term(sqrt, 2)
# @test solve(expr, x) == Symbolics.term(sqrt, 2)

# solve errors
# expr = x + im
# @test solve(expr, x) == -im

# Factorisation #

f = 10x
u, factors = RootFinding.factor_use_nemo(f)
@test isequal(u, 10) && isequal(factors, [x])

f = Symbolics.wrap(10)
u, factors = RootFinding.factor_use_nemo(f)
@test isequal(u, 10) && isempty(factors)

f = x^2 - 1
u, factors = RootFinding.factor_use_nemo(f)
@test isequal(u, 1) && isequal(expand(u*prod(factors) - f), 0)

f = expand((x + 1//3) * ((x*y)^2 + 2x*y + y^2) * (x - z))
u, factors = RootFinding.factor_use_nemo(f)
@test isequal(expand(u*prod(factors) - f), 0)



# Gcd #

f1, f2 = x^2 - y^2, x^3 - y^3
@test isequal(x - y, RootFinding.gcd_use_nemo(f1, f2))
