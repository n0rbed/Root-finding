using Symbolics

coeff = Symbolics.coeff
function solve_poly(expression, x)
    if isequal(SymbolicUtils.operation(expression.val), ^) && SymbolicUtils.arguments(expression.val)[2] isa Int64
        expression = SymbolicUtils.arguments(expression.val)[1]
    end

    expression = expand(expression)
    expression = simplify.(expression)
    degree = Symbolics.degree(expression, x)


    if degree == 0 && expression == 0
        return 0
    elseif degree == 0 && expression != 0
        return "Not a valid statement"
    end

    if degree == 1
        # mx + c = 0
        coeffs, constant = polynomial_coeffs(expression, [x])
        m = coeffs[x]
        c = coeffs[x^0]
        root = -c / m
        return root
    end

    if degree == 2
        # ax^2 + bx + c
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
ex = (x^2 + 2*x + 1)^20
println("x --> ", solve_poly(ex, x))
  # Output --> 2.0
