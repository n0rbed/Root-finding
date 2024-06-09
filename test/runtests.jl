using RootFinding, Symbolics
using Test

@variables x y z

@test isequal(solve(x+1, x), -1)

@test isequal(solve(2x+1, x), -1/2)

@test isequal(solve(x, x), 0) 

@test isequal(solve(x^2 + 1, x), [-1, -1]) 

@test isequal(solve(x^2 + 2x + 10) ,[-1 + 3im, -1 - 3im])

@test isequal(solve(x^2 - 10x + 25) ,[5, 5])

@test isequal(solve((x+1)^20, x), -1)

@test isequal(solve(x^3 - 2x^2 + x - 2, x), [2, -im, im])

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
