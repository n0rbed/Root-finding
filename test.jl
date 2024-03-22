
using Symbolics

coeff = Symbolics.coeff
function solve_poly(expression, x)
    simple_expression = expand(expression)
    simple_expression = simplify.(simple_expression)
    degree = Symbolics.degree(simple_expression, x)


    expression = simple_expression




    if degree == 0 && expression == 0
        return 0
    elseif degree == 0 && expression != 0
        return "Not a valid statement"
    end

    if degree == 1
        m = Symbolics.coeff(expression, x)
        c = coeff(expression)
        root = -c / m
        return root
    end

    if degree == 2
        coeffs, constant = polynomial_coeffs(expression, [x])

        a = coeffs[x^2]
        b = coeffs[x]
        c = coeffs[x^0]

        root1 = simplify(expand((-b + Symbolics.Term(sqrt, [(b^2 - 4(a*c))])) / 2a))
        root2 = simplify(expand((-b - Symbolics.Term(sqrt, [(b^2 - 4(a*c))])) / 2a))

        return [root1, root2]
    end
end

function solve_poly(polys::Vector, x)
    polys = unique(polys)
    expression = polys[1]
    for i = 2:(length(polys))
        expression = expression - polys[i]
    end
    expression = expand(expression)
    expression = simplify.(expression)
    return solve_poly(expression, x)
end

@variables x a
ex = x^2 + 2a*x + a^2
println("x --> ", solve_poly(ex, x))
  # Output --> 2.0
