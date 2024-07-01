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

@test isequal(get_roots_deg1(x + y^3, x), [-y^3])
@test_throws DomainError get_roots_deg1(x, x^2)
@test_throws DomainError get_roots_deg1(x + sin(x), x)
@test_throws DomainError get_roots_deg1(x^2, x)

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
