using Test, RootFinding, Symbolics

@variables x

@test isequal(solvepoly(x - 1, x), 1)

