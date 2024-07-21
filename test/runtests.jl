using RootFinding, Symbolics
using Test

E = Base.MathConstants.e

function sort_roots(roots)
    return sort(roots, lt = (x,y) -> real(x)==real(y) ? imag(x)<imag(y) : real(x)<real(y))
end

function sort_arr(sols, vars)
    for i = 1:length(sols)
        sols[i] = convert(Dict{Num, Any}, sols[i])
        for (var, root) in sols[i]
            sols[i][var] = eval(Symbolics.toexpr(root))
        end
    end

    extract_root(c) = (real(c), imag(c))
    sort_root(dict) = Tuple(extract_root(dict[var]) for var in vars)

    # Sorting the answers lexicographically
    sols = sort(sols, by=sort_root)
end

function check_equal(arr1, arr2)
    l1 = length(arr1)
    if l1 != length(arr2)
        return false
    end
    for i = 1:l1
        if !isequal(keys(arr1[i]), keys(arr2[i]))
            return false
        end
        if !all(values(arr1[i]) .≈ values(arr2[i]))
            return false
        end
    end
    return true
end

@variables x y z

@testset "Invalid input" begin
    @test_throws AssertionError get_roots(x, x^2)
    @test_throws AssertionError get_roots(x^3 + sin(x), x)
    @test_throws AssertionError get_roots(1/x, x)
end

@testset "Deg 1 univar" begin
    @test isequal(solve(x+1, x), [-1])

    @test isequal(solve(2x+1, x), [-1/2])

    @test isequal(solve(x, x), [0]) 

    @test isequal(solve((x+1)^20, x), [-1])

    @test isequal(get_roots_deg1(x + y^3, x), [-y^3])

    expr = x - Symbolics.term(sqrt, 2)
    @test isequal(solve(expr, x)[1], Symbolics.term(sqrt, 2))

    expr = x + im
    @test solve(expr, x)[1] == -im
end

@testset "Deg 2 univar" begin
    expr = x^2 + 1
    arr_get_roots = sort_roots(eval.(Symbolics.toexpr.(get_roots_deg2(expr, x))))
    arr_solve_roots = sort_roots(eval.(Symbolics.toexpr.(solve(expr, x))))
    arr_known_roots = sort_roots([-im, im])
    @test all(arr_get_roots .≈ arr_known_roots)
    @test all(arr_solve_roots .≈ arr_known_roots)

    expr = x^2 + 2x + 10
    arr_get_roots = sort_roots(eval.(Symbolics.toexpr.(get_roots_deg2(expr, x))))
    arr_solve_roots = sort_roots(eval.(Symbolics.toexpr.(solve(expr, x))))
    arr_known_roots = sort_roots([-1 + 3im, -1 - 3im])
    @test all(arr_get_roots .≈ arr_known_roots)
    @test all(arr_solve_roots .≈ arr_known_roots)

    expr = x^2 - 10x + 25
    arr_get_roots = sort_roots(eval.(Symbolics.toexpr.(get_roots_deg2(expr, x))))
    arr_solve_roots = sort_roots(eval.(Symbolics.toexpr.(solve(expr, x))))
    arr_known_roots = [5,5]
    @test all(arr_get_roots .≈ arr_known_roots)
    @test all(arr_solve_roots .≈ arr_known_roots)
end

@testset "Deg 3 univar" begin
    expr = x^3 - 2x^2 + x - 2 
    arr_get_roots = sort_roots(eval.(Symbolics.toexpr.(get_roots_deg3(expr, x))))
    arr_solve_roots = sort_roots(eval.(Symbolics.toexpr.(solve(expr, x))))
    arr_known_roots = sort_roots([2, -im, im])
    @test all(isapprox.(arr_get_roots, arr_known_roots, atol=0.0000001))
    @test all(arr_solve_roots .≈ arr_known_roots)

    expr = x^3 + x^2 + x + 1
    arr_get_roots = sort_roots(eval.(Symbolics.toexpr.(get_roots_deg3(expr, x))))
    arr_solve_roots = sort_roots(eval.(Symbolics.toexpr.(solve(expr, x))))
    arr_known_roots = sort_roots([-1, -im, im])
    @test all(isapprox.(arr_get_roots, arr_known_roots, atol=0.0000001))
    @test all(arr_solve_roots .≈ arr_known_roots)

    expr = x^3 + 10x
    arr_get_roots = sort_roots(eval.(Symbolics.toexpr.(get_roots_deg3(expr, x))))
    arr_solve_roots = sort_roots(eval.(Symbolics.toexpr.(solve(expr, x))))
    arr_known_roots = sort_roots([0, -sqrt(10)*im, sqrt(10)*im])
    @test all(arr_get_roots .≈ arr_known_roots)
    @test all(arr_solve_roots .≈ arr_known_roots)
end

@testset "Deg 4 univar" begin
    expr = x^4 + 1
    arr_get_roots = sort_roots(eval.(Symbolics.toexpr.(get_roots_deg4(expr, x))))
    arr_solve_roots = sort_roots(eval.(Symbolics.toexpr.(solve(expr, x))))
    arr_known_roots = sort_roots(eval.([-(complex(-1))^(1/4),(complex(-1))^(1/4), (complex(-1))^(3/4), -(complex(-1))^(3/4)]))
    @test all(arr_get_roots .≈ arr_known_roots)
    @test all(arr_solve_roots .≈ arr_known_roots)

    expr = x^4 - 3x^2 + 2
    arr_get_roots = sort_roots(eval.(Symbolics.toexpr.(get_roots_deg4(expr, x))))
    arr_solve_roots = sort_roots(eval.(Symbolics.toexpr.(solve(expr, x))))
    arr_known_roots = sort_roots(eval.([-1, 1, sqrt(2), -sqrt(2)]))
    @test all(arr_get_roots .≈ arr_known_roots)
    @test all(arr_solve_roots .≈ arr_known_roots)

    expr = x^4 - x^3 - 2x^2 + 6x - 4
    arr_get_roots = sort_roots(eval.(Symbolics.toexpr.(get_roots_deg4(expr, x))))
    arr_solve_roots = sort_roots(eval.(Symbolics.toexpr.(solve(expr, x))))
    arr_known_roots = sort_roots(eval.([-2, 1, 1-im, 1+im]))
    @test all(arr_get_roots .≈ arr_known_roots)
    @test all(arr_solve_roots .≈ arr_known_roots)

    expr = 386314 - 412163x - 357800(x^2) + 1029179(x^3) - 111927(x^4)
    arr_get_roots = sort_roots(eval.(Symbolics.toexpr.(get_roots_deg4(expr, x))))
    arr_solve_roots = sort_roots(eval.(Symbolics.toexpr.(solve(expr, x))))
    arr_known_roots = sort_roots(eval.([0.5840484 + 0.4176250im, 0.5840484 - 0.4176250im,
    8.788773679354421, -0.76177906049]))
    @test all(isapprox.(arr_known_roots, arr_get_roots, atol=0.00001))
    @test all(isapprox.(arr_known_roots, arr_solve_roots, atol=0.00001))
end

@testset "Complex coeffs univar" begin
    expr = x^4 + sqrt(complex(-2//1))
    arr_get_roots = sort_roots(eval.(Symbolics.toexpr.(solve(expr, x))))
    arr_known_roots = sort_roots(eval.([-Complex(-1)^(3/8)*2^(1/8), Complex(-1)^(3/8)*2^(1/8), 
        Complex(-1)^(7/8)*2^(1/8), -Complex(-1)^(7/8)*2^(1/8)]))
    @test all(arr_get_roots .≈ arr_known_roots)

    # standby
    # expr = x^3 + sqrt(complex(-2//1))*x + 2
end




@testset "Multivar solver" begin
    eqs = [x*y + 2x^2, y^2 -1]
    arr_calcd_roots = sort_arr(solve(eqs, [x,y]), [x,y])
    arr_known_roots = [Dict(x=>-1/2, y=>1), Dict(x=>0, y=>-1), Dict(x=>0, y=>1), Dict(x=>1/2, y=>-1)]
    arr_known_roots = sort_arr(arr_known_roots, [x,y])
    @test check_equal(arr_calcd_roots, arr_known_roots)   

    eqs = [x-y-z, x+y-z^2, x^2 + y^2 - 1]
    arr_calcd_roots = sort_arr(solve(eqs, [x,y,z]), [x,y,z])
    arr_known_roots = sort_arr([Dict(x => 0, y=>1, z=>-1), Dict(x=>1, y=>0, z=>1),
        Dict(x=>(1/2)*(-2-sqrt(2)*im), y=>(1/2)*(-2+sqrt(2)*im), z=>-sqrt(2)*im),
        Dict(x=>(1/2)*(-2+sqrt(2)*im), y=>(1/2)*(-2-sqrt(2)*im), z=>sqrt(2)*im)], [x,y,z])
    @test check_equal(arr_calcd_roots, arr_known_roots)   

    eqs = [x^2, y, z]
    arr_calcd_roots = sort_arr(solve(eqs, [x,y,z], true), [x,y,z])
    arr_known_roots = sort_arr([Dict(x=>0, y=>0, z=>0), Dict(x=>0, y=>0, z=>0)], [x,y,z])
    @test check_equal(arr_calcd_roots, arr_known_roots)   

    eqs = [y^2 - 1, x]
    arr_calcd_roots = sort_arr(solve(eqs, [x,y]), [x,y])
    arr_known_roots = sort_arr([Dict(y=>1//1, x=>0//1), Dict(y=>-1//1, x=>0//1)], [x,y])
    @test check_equal(arr_calcd_roots, arr_known_roots)   

    eqs = [x^5 + x, y]
    arr_calcd_roots = sort_arr(solve(eqs, [x,y]), [x,y])
    arr_known_roots = sort_arr([Dict(x=>0, y=>0), Dict(x=>-(complex(-1))^(1/4), y=>0),
    Dict(x=>(complex(-1))^(1/4), y=>0), Dict(x=>-(complex(-1))^(3/4), y=>0),
    Dict(x=>(complex(-1))^(3/4), y=>0)], [x,y])
    @test check_equal(arr_calcd_roots, arr_known_roots)   

    @test isequal(solve([x*y - 1, y], [x,y]), [])
    @test isequal(solve([x+y+1, x+y+2], [x,y]), [])

end


@testset "Factorisation" begin
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
end

@testset "GCD" begin
    f1, f2 = x^2 - y^2, x^3 - y^3
    @test isequal(x - y, RootFinding.gcd_use_nemo(f1, f2))
end


# Post Process roots #


# Filter Poly # For future reference, i think @check_polynomial is not enough of a test.
@testset "Filter poly" begin
    @variables x y z c1 c2 
    poly = x*sqrt(complex(-2)) + 2.23324234
    subs, filtered_poly = filter_poly(poly, x)
    @test check_polynomial(filtered_poly)

    poly = x + 2im
    subs, filtered_poly = filter_poly(poly, x)
    @test check_polynomial(filtered_poly)

    poly = im*x + Symbolics.wrap(Symbolics.term(sqrt, 2))
    subs, filtered_poly = filter_poly(poly, x)
    @test check_polynomial(filtered_poly)

    poly = (1/im)*x + 3*y*z
    subs, filtered_poly = filter_poly(poly, x)
    @test check_polynomial(filtered_poly)

    poly = (x+1)*Symbolics.term(log, 3)
    subs, filtered_poly = filter_poly(poly, x)
    @test check_polynomial(filtered_poly)
end


@testset "n func occ (smart occurrence counter)" begin
    @test n_func_occ(x, x) == 1
    @test n_func_occ(log(x), x) == 1
    @test n_func_occ(log(x) + x, x) == 2

    # standby
    # @test n_func_occ(log(y) + x , x) == 1
    
    @test n_func_occ(log(x + sin((x^2 + x)/log(x))), x) == 3
    @test n_func_occ(x^2 + x + x^3, x) == 1
    @test n_func_occ(log(x)^2 - 17, x) == 1
    @test n_func_occ(2^(x^2 + x) + 5^(x+3), x) == 2

    expr = log( log(x) + log(x) ) + log( log(x) + log(x) ) - 11
    @test n_func_occ(expr, x) == 1

    # log(2) - 3log(5) + x*log(2) - x*log(5)
    expr = expand((1 + x)*Symbolics.term(log, 2) - (3 + x)*Symbolics.term(log, 5))
    @test n_func_occ(expr, x) == 1
end



@testset "NL solve" begin
    @variables a b c d e x
    lhs = ia_solve(a*x^b + c, x)[1]
    lhs2 = solve(a*x^b + c, x)[1]
    rhs = Symbolics.term(^, -c.val/a.val, 1/b.val) 
    #@test isequal(lhs, rhs)
    
    expr = x + 2
    lhs = eval.(Symbolics.toexpr.(ia_solve(expr, x)))
    lhs_solve = eval.(Symbolics.toexpr.(solve(expr, x)))
    rhs = [-2]
    @test lhs[1] ≈ rhs[1]
    @test lhs_solve[1] ≈ rhs[1]

    expr = sqrt(log(cbrt(x^2)))
    lhs = sort_roots(eval.(Symbolics.toexpr.(ia_solve(expr, x))))
    lhs_solve = sort_roots(eval.(Symbolics.toexpr.(solve(expr, x))))
    rhs = sort_roots([1, -1])
    @test isequal(lhs, rhs)

    expr = 2^(x+1) + 5^(x+3)
    lhs = eval.(Symbolics.toexpr.(ia_solve(expr, x)))
    lhs_solve = eval.(Symbolics.toexpr.(solve(expr, x)))
    rhs = [(-im*Base.MathConstants.pi - log(2) + 3log(5))/(log(2) - log(5))]
    @test lhs[1] ≈ rhs[1]
    @test lhs_solve[1] ≈ rhs[1]

    expr = 3*2^(x+3) + 2*5^(x+1)
    lhs = eval.(Symbolics.toexpr.(ia_solve(expr, x)))
    lhs_solve = eval.(Symbolics.toexpr.(solve(expr, x)))
    rhs = [(-im*Base.MathConstants.pi + log(5) - log(12))/(log(2) - log(5))]
    @test lhs[1] ≈ rhs[1]
    @test lhs_solve[1] ≈ rhs[1]

    expr = exp(2x)*exp(x^2 + 3) + 3
    lhs = sort_roots(eval.(Symbolics.toexpr.(ia_solve(expr, x))))
    lhs_solve = sort_roots(eval.(Symbolics.toexpr.(solve(expr, x))))
    rhs = sort_roots([-1 + sqrt(-2 + im*Base.MathConstants.pi + log(3)),
        -1 - sqrt(-2 + im*Base.MathConstants.pi + log(3))])
    @test all(lhs .≈ rhs)
    @test all(lhs_solve .≈ rhs)

    expr = x/5 + 3x^2
    lhs = sort_roots(eval.(Symbolics.toexpr.(ia_solve(expr, x))))
    lhs_solve = sort_roots(eval.(Symbolics.toexpr.(solve(expr, x))))
    rhs = sort_roots([0, -1//15])
    @test all(lhs .≈ rhs)
    @test all(lhs_solve .≈ rhs)


    expr = log(x)^2 + log(x) + 1
    lhs = sort_roots(eval.(Symbolics.toexpr.(ia_solve(expr, x))))
    lhs_solve = sort_roots(eval.(Symbolics.toexpr.(solve(expr, x))))
    rhs = sort_roots([E^(-(complex(-1))^(1/3)), E^((complex(-1))^(2/3))])
    @test all(lhs .≈ rhs)
    @test all(lhs_solve .≈ rhs)

    expr = 2log(x^2 + 1)^2 + 3log(x^2 + 1) + 1
    lhs = sort_roots(eval.(Symbolics.toexpr.(ia_solve(expr, x))))
    lhs_solve = sort_roots(eval.(Symbolics.toexpr.(solve(expr, x))))
    rhs = sort_roots([-im*ssqrt(1 - 1/E), im*ssqrt(1 - 1/E),
        -im*ssqrt(1 - 1/ssqrt(E)), im*ssqrt(1 - 1/ssqrt(E))])
    @test all(lhs .≈ rhs)
    @test all(lhs_solve .≈ rhs)


    # this has c1 subbed, ia_solve still needs to incl c1
    # much like sin and cos solutions (i.e. infinite solutions)
    expr = 9^(x^2 + 1) + 3^(x^2 + 1) + 2
    lhs = sort_roots(eval.(Symbolics.toexpr.(ia_solve(expr, x))))
    lhs_solve = sort_roots(eval.(Symbolics.toexpr.(solve(expr, x))))
    rhs = sort_roots([-ssqrt(-1 + slog(1/2*(-1 - im*ssqrt(7)))/slog(3)),
        ssqrt(-1 + slog(1/2*(-1 - im*ssqrt(7)))/slog(3)),
        ssqrt(-1 + slog(1/2*(-1 + im*ssqrt(7)))/slog(3)),
        -ssqrt(-1 + slog(1/2*(-1 + im*ssqrt(7)))/slog(3))])

    @test all(lhs .≈ rhs)
    @test all(lhs_solve .≈ rhs)
end

@testset "Turn to poly" begin
    @variables x
    # does not sub because these can not be solved as polys
    expr = log(x+1)^2 + log(x) + 1
    expr, sub = turn_to_poly(expr, x)
    @test isequal(sub, Dict())

    expr = 2log(x)^2 + sin(x) + 1
    expr, sub = turn_to_poly(expr, x)
    @test isequal(sub, Dict())

    
    # subs and turns to polys
    expr = log(x)^2 + log(x) + 1
    expr, sub = turn_to_poly(expr, x)
    @test check_polynomial(expr)

    expr = 2log(x^2 + 1)^2 + 3log(x^2 + 1) + 1
    expr, sub = turn_to_poly(expr, x)
    @test check_polynomial(expr)

    expr = sin(x^2)^2 + 3sin(x^2) + 1
    expr, sub = turn_to_poly(expr, x)
    @test check_polynomial(expr)

    expr = 9^(x^2 + 1) + 3^(x^2 + 1) + 2
    expr, sub = turn_to_poly(expr, x)
    @test check_polynomial(expr)
end


