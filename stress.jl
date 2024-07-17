using Symbolics, RootFinding
using Test

######## POLYNOMIALS OF DEGREE <=4 ########

#=
Bad examples found:

f = 1
CRASH

f = 386314 - 412163x - 357800(x^2) + 1029179(x^3) - 111927(x^4)
OVERFLOW
=#

@variables x

boot = 10

for i in 1:boot
    deg = rand(1:4)
    f = rand(-2^20:2^20)
    for j in 1:deg
        f = f * x
        f += rand(-2^20:2^20)
    end
    f = expand(f)
    @info "" f
    roots = RootFinding.solve(f, x)
    @info "" roots
    roots_eval = eval.(Symbolics.toexpr.(roots))
    f_eval = map(root -> substitute(f, x => root), roots_eval)
    @info "" f_eval
    @test all(isapprox.(f_eval, 0, atol=1e-10))
end

######## PARAMETRIC POLYNOMIALS OF DEGREE <=4 ########

# TODO
# @variables x a b c d e
# params = [a,b,c,d,e]

# boot = 10

# for i in 1:boot
#     deg = rand(1:4)
#     f = rand(-100:100) * params[1]
#     for j in 1:deg
#         f = f * x
#         f += rand(-100:100) * params[j + 2]
#     end
#     f = expand(f)
#     @info "" f
#     roots = RootFinding.solve(f, x)
#     @info "" roots
# end

