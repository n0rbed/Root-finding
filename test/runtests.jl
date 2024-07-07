using RootFinding, Symbolics
using Test

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

# univar 
@test isequal(solve(x+1, x), [-1])

@test isequal(solve(2x+1, x), [-1/2])

@test isequal(solve(x, x), [0]) 

@test isequal(solve((x+1)^20, x), [-1])

@test isequal(get_roots_deg1(x + y^3, x), [-y^3])

exp = x^2 + 1
arr_get_roots = sort_roots(eval.(Symbolics.toexpr.(get_roots_deg2(exp, x))))
arr_solve_roots = sort_roots(eval.(Symbolics.toexpr.(solve(exp, x))))
arr_known_roots = sort_roots([-im, im])
@test all(arr_get_roots .≈ arr_known_roots)
@test all(arr_solve_roots .≈ arr_known_roots)

exp = x^2 + 2x + 10
arr_get_roots = sort_roots(eval.(Symbolics.toexpr.(get_roots_deg2(exp, x))))
arr_solve_roots = sort_roots(eval.(Symbolics.toexpr.(solve(exp, x))))
arr_known_roots = sort_roots([-1 + 3im, -1 - 3im])
@test all(arr_get_roots .≈ arr_known_roots)
@test all(arr_solve_roots .≈ arr_known_roots)

exp = x^2 - 10x + 25
arr_get_roots = sort_roots(eval.(Symbolics.toexpr.(get_roots_deg2(exp, x))))
arr_solve_roots = sort_roots(eval.(Symbolics.toexpr.(solve(exp, x))))
arr_known_roots = [5,5]
@test all(arr_get_roots .≈ arr_known_roots)
@test all(arr_solve_roots .≈ arr_known_roots)


exp = x^3 - 2x^2 + x - 2 
arr_get_roots = sort_roots(eval.(Symbolics.toexpr.(get_roots_deg3(exp, x))))
arr_solve_roots = sort_roots(eval.(Symbolics.toexpr.(solve(exp, x))))
arr_known_roots = sort_roots([2, -im, im])
@test all(arr_get_roots .≈ arr_known_roots)
@test all(arr_solve_roots .≈ arr_known_roots)

exp = x^3 + x^2 + x + 1
arr_get_roots = sort_roots(eval.(Symbolics.toexpr.(get_roots_deg3(exp, x))))
arr_solve_roots = sort_roots(eval.(Symbolics.toexpr.(solve(exp, x))))
arr_known_roots = sort_roots([-1, -im, im])
@test all(arr_get_roots .≈ arr_known_roots)
@test all(arr_solve_roots .≈ arr_known_roots)

exp = x^3 + 10x
arr_get_roots = sort_roots(eval.(Symbolics.toexpr.(get_roots_deg3(exp, x))))
arr_solve_roots = sort_roots(eval.(Symbolics.toexpr.(solve(exp, x))))
arr_known_roots = sort_roots([0, -sqrt(10)*im, sqrt(10)*im])
@test all(arr_get_roots .≈ arr_known_roots)
@test all(arr_solve_roots .≈ arr_known_roots)

exp = x^4 + 1
arr_get_roots = sort_roots(eval.(Symbolics.toexpr.(get_roots_deg4(exp, x))))
arr_solve_roots = sort_roots(eval.(Symbolics.toexpr.(solve(exp, x))))
arr_known_roots = sort_roots(eval.([-(complex(-1))^(1/4),(complex(-1))^(1/4), (complex(-1))^(3/4), -(complex(-1))^(3/4)]))
@test all(arr_get_roots .≈ arr_known_roots)
@test all(arr_solve_roots .≈ arr_known_roots)

exp = x^4 - 3x^2 + 2
arr_get_roots = sort_roots(eval.(Symbolics.toexpr.(get_roots_deg4(exp, x))))
arr_solve_roots = sort_roots(eval.(Symbolics.toexpr.(solve(exp, x))))
arr_known_roots = sort_roots(eval.([-1, 1, sqrt(2), -sqrt(2)]))
@test all(arr_get_roots .≈ arr_known_roots)
@test all(arr_solve_roots .≈ arr_known_roots)

exp = x^4 - x^3 - 2x^2 + 6x - 4
arr_get_roots = sort_roots(eval.(Symbolics.toexpr.(get_roots_deg4(exp, x))))
arr_solve_roots = sort_roots(eval.(Symbolics.toexpr.(solve(exp, x))))
arr_known_roots = sort_roots(eval.([-2, 1, 1-im, 1+im]))
@test all(arr_get_roots .≈ arr_known_roots)
@test all(arr_solve_roots .≈ arr_known_roots)

#complex univar
exp = x^4 + sqrt(complex(-2//1))
arr_get_roots = sort_roots(eval.(Symbolics.toexpr.(solve(exp, x))))
arr_known_roots = sort_roots(eval.([-Complex(-1)^(3/8)*2^(1/8), Complex(-1)^(3/8)*2^(1/8), 
    Complex(-1)^(7/8)*2^(1/8), -Complex(-1)^(7/8)*2^(1/8)]))
@test all(arr_get_roots .≈ arr_known_roots)


# standby
# exp = x^3 + sqrt(complex(-2//1))*x + 2

# multivar
eqs = [x*y + 2x^2, y^2 -1]
arr_calcd_roots = sort_arr(solve(eqs, [x,y]), [x,y])
arr_known_roots = [Dict(x=>-1/2, y=>1), Dict(x=>0, y=>-1), Dict(x=>0, y=>1), Dict(x=>1/2, y=>-1)]
arr_known_roots = sort_arr(arr_known_roots, [x,y])
@test check_equal(arr_calcd_roots, arr_known_roots)   

@test_throws AssertionError get_roots(x, x^2)
@test_throws AssertionError get_roots(x^3 + sin(x), x)
@test_throws AssertionError get_roots(1/x, x)


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

# solve returns un-evaluated sqrt(2)
expr = x - Symbolics.term(sqrt, 2)
@test isequal(solve(expr, x)[1], Symbolics.term(sqrt, 2))

# solve errors
 expr = x + im
 #@test solve(expr, x) == -im

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


# Post Process roots #


# Filter Poly #
@variables x y z c1 c2 
poly = x*sqrt(complex(-2)) + 2.23324234
subs, filtered_poly = filter_poly(poly, x)
@test check_polynomial(filtered_poly)

poly = x + 2im
subs, filtered_poly = filter_poly(poly, x)
@test check_polynomial(filtered_poly)

# note: poly = im*x + Symbolics.term(sqrt, 2) fails
poly = (1/im)*x + 3*y*z
subs, filtered_poly = filter_poly(poly, x)
@test check_polynomial(filtered_poly)

poly = (x+1)*Symbolics.term(log, 3)
subs, filtered_poly = filter_poly(poly, x)
@test check_polynomial(filtered_poly)


# n_func_occ #
@test n_func_occ(x, x) == 1
@test n_func_occ(log(x), x) == 1
@test n_func_occ(log(x) + x, x) == 2
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



# NL SOLVE #
@variables a b c d e x
lhs = nl_solve(a*x^b + c, x)
rhs = Symbolics.term(^, -c/a, 1/b) 
@test isequal(lhs, rhs)

expr = sqrt(log(cbrt(x^2)))
lhs = sort_roots(eval.(Symbolics.toexpr.(nl_solve(expr, x))))
rhs = sort_roots([1, -1])
@test isequal(lhs, rhs)



