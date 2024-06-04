using Symbolics, Groebner, SymbolicUtils

coeff = Symbolics.coeff
function solve(expression, x)
    if isequal(imag.(expression), 0)
        expression = real.(expression)
        if isequal(SymbolicUtils.operation(expression.val), ^) && SymbolicUtils.arguments(expression.val)[2] isa Int64
            expression = SymbolicUtils.arguments(expression.val)[1]
        end
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
        m = get(coeffs, x, 0)
        c = get(coeffs, x^0, 0)
        root = -c / m
        return root
    end

    if degree == 2
        # ax^2 + bx + c
        coeffs, constant = polynomial_coeffs(expression, [x])

        a = coeffs[x^2]
        b = get(coeffs, x, 0)
        c = get(coeffs, x^0, 0)

        root1 = simplify(expand((-b + Symbolics.Term(sqrt, [(b^2 - 4(a*c))])) / 2a))
        root2 = simplify(expand((-b - Symbolics.Term(sqrt, [(b^2 - 4(a*c))])) / 2a))

        return [root1, root2]
    end

    if degree == 3
        coeffs, constant = polynomial_coeffs(expression, [x])

        a = coeffs[x^3]
        b = get(coeffs, x^2, 0)
        c = get(coeffs, x, 0)
        d = get(coeffs, x^0, 0)

        Q = ((3*a*c) - b^2)//(9a^2)
        R = (9*a*b*c - ((27*(a^2)*d)+2b^3))//(54a^3)

        S = Symbolics.term(cbrt, (R + Symbolics.term(sqrt, (Q^3+R^2))))
        T = Symbolics.term(complex, Symbolics.term(cbrt, (R - Symbolics.term(sqrt, (Q^3+R^2)))))

        root1 = S + T - (b//(3*a))
        root2 = -((S+T)//2) - (b//(3*a)) + (im*(Symbolics.term(sqrt, 3))/2)*(S-T)
        root3 = -((S+T)//2) - (b//(3*a)) - (im*(Symbolics.term(sqrt, 3))/2)*(S-T)
        
        if imag.(eval(Symbolics.toexpr(root2))) == 0
            root2 = real.(root2)
        end
        if imag.(eval(Symbolics.toexpr(root3))) == 0
            root3 = real.(root3)
        end
        return [real.(root1), root2, root3]
    end

    if degree == 4

        coeffs, constant = polynomial_coeffs(expression, [x])

        a = coeffs[x^4]
        b = get(coeffs, x^3, 0)
        c = get(coeffs, x^2, 0)
        d = get(coeffs, x, 0)
        e = get(coeffs, x^0, 0)

        p = (8(a*c)-3(b^2))/(8(a^2))

        q = (b^3 - 4(a*b*c) + 8(d*a^2))/(8*a^3)

        r = (-3(b^4) + 256(e*a^3) - 64(d*b*a^2) + 16(c*(b^2)*a))/(256*a^4)

        @variables m y
        eq_m = 8m^3 + 8(p)*m^2 + (2(p^2) - 8r)m - q^2
        roots_m = solve(eq_m, m)
        m = 0
        for root in roots_m
            if root != 0
                m = root
                break
            end
        end

        root1, root2 = solve(y^2 + (p/2) + m + ((2m)^(1/2))*y - (q/(2*((2m)^(1/2)))), y)
        root3, root4 = solve(y^2 + (p/2) + m - ((2m)^(1/2))*y + (q/(2*((2m)^(1/2)))), y)

        arr = [root1, root2, root3, root4]
        for (i, root) in enumerate(arr)
            arr[i] = root - (b/(4a))
        end

        return arr
    end

end

function solve(polys::Vector, x::Num)
    polys = unique(polys)
    expression = polys[1]
    for i = 2:(length(polys))
        expression = expression - polys[i]
    end
    expression = expand(expression)
    expression = simplify.(expression)
    return solve(expression, x)
end




function contains(var, vars)
    for variable in vars
        if isequal(var, variable)
            return true
        end
    end
    return false
end


function solve(eqs::Vector{Num}, vars::Vector{Num})
    eqs = convert(Vector{Any}, Symbolics.groebner_basis(eqs))
    all_roots = Dict()

    # initialize the place of each var
    for var in vars
        all_roots[var] = [] 
    end


    #get roots for first var (z in this case)
    solved = false
    while !solved
        # first, solve any single variable equations
        for eq in eqs
            for var in vars
                present_vars = Symbolics.get_variables(eq)
                if contains(var, present_vars) && size(present_vars, 1) == 1 
                    var = Symbolics.get_variables(eq)[1]
                    append!(all_roots[var], solve(Symbolics.wrap(eq), var))
                end
            end
        end

        # second, Substitute the roots of the variables where found
        for (i, eq) in pairs(eqs)
            for var in vars
                present_vars = Symbolics.get_variables(eq)

                if contains(var, present_vars) && !isequal(all_roots[var], [])
                    deleteat!(eqs, i)

                    for root in all_roots[var]
                        insert!(eqs, i, substitute(eq, Dict([var => root])))
                    end

                end
            end
        end


        solved = true
        for (var, value) in all_roots
            if isequal(value, [])
                solved = false
            end
        end
    end
    return all_roots
end

@variables x
print(solve(x^3 -2x^2 + x - 2, x))