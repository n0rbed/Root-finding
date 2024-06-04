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
