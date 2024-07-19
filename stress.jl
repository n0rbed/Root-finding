using Symbolics, RootFinding
using Test

######## POLYNOMIALS OF DEGREE <=4 ########

@variables x

if false # switch to true to run
    boot = 100

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
        @test all(isapprox.(f_eval, 0, atol=1e-5))
    end
end

######## PARAMETRIC POLYNOMIALS OF DEGREE <=4 ########

#=
Bad example found:

RootFinding.solve(3b, x)
CRASH (solved)

RootFinding.solve(-19e + 93d*x + 7c*(x^2) + 45b*(x^3) - 37a*(x^4), x)
OVERFLOW
=#

if false # switch to true to run
    @variables x a b c d e
    params = [a,b,c,d,e]

    boot = 1000

    for i in 1:boot
        deg = rand(1:4)
        f = rand(-100:100) * params[1]
        for j in 1:deg
            f = f * x
            f += rand(-100:100) * params[j + 1]
        end
        f = expand(f)
        @info "" f
        roots = RootFinding.solve(f, x)
        @info "" roots
        isempty(roots) && (@assert length(f) == 1; true) && continue
        # getting 0 by chance is very unlucky
        point = Dict(params .=> rand(length(params)))
        @info "" point
        roots_eval = map(root -> Symbolics.substitute(root, point), roots)
        f_eval = map(root -> substitute(substitute(f, x => root), point), roots_eval)
        @info "" f_eval
        @test all(isapprox.(f_eval, 0, atol=1e-5))
    end
end

######## PARAMETRIC POLYNOMIALS ########

#=

=#

if true # switch to true to run
    @variables x p[1:100]

    boot = 1000

    for i in 1:boot
        # make a large random polynomial as a product of many tiny polynomials
        f = 1
        nfactors = rand(2:4)
        for j in 1:nfactors
            deg = rand(1:4)
            g = rand(-10:10) * p[rand(1:10)]
            for k in 1:deg
                g = g * x
                g += rand(-10:10) * p[rand(1:10)]
            end
            f = f*expand(g)
        end
        @info "" f expand(f)
        f = expand(f)
        roots = RootFinding.solve(f, x)
        @info "" roots
        isempty(roots) && (@assert length(f) == 1; true) && continue
        # getting 0 by chance is very unlucky
        point = Dict(collect(p) .=> rand(length(p)))
        @info "" point
        roots_eval = map(root -> Symbolics.substitute(root, point), roots)
        f_eval = map(root -> substitute(substitute(f, x => root), point), roots_eval)
        @info "" f_eval
        @test all(isapprox.(f_eval, 0, atol=1e-2))
    end
end
