using Symbolics
coeff = Symbolics.coeff

function get_roots_deg2(expression, x)
    # ax^2 + bx + c = 0
    coeffs, constant = polynomial_coeffs(expression, [x])

    a = rationalize(coeffs[x^2])
    b = rationalize(get(coeffs, x, 0))
    c = rationalize(get(coeffs, x^0, 0))


    root1 = simplify((-b + Symbolics.term(sqrt, (b^2 - 4(a*c)))) / 2a)
    root2 = simplify((-b - Symbolics.term(sqrt, (b^2 - 4(a*c)))) / 2a)
    if eval(Symbolics.toexpr(b^2 - 4(a*c))) < 0
        root1 = simplify((-b + Symbolics.term(sqrt, Symbolics.term(complex, (b^2 - 4(a*c))))) / 2a)
        root2 = simplify((-b - Symbolics.term(sqrt, Symbolics.term(complex, (b^2 - 4(a*c))))) / 2a)
    end

    return [root1, root2]
end

function get_roots_deg3(expression, x)
    coeffs, constant = polynomial_coeffs(expression, [x])

    a = rationalize(coeffs[x^3])
    b = rationalize(get(coeffs, x^2, 0))
    c = rationalize(get(coeffs, x, 0))
    d = rationalize(get(coeffs, x^0, 0))

    Q = ((3*a*c) - b^2)//(9a^2)
    R = (9*a*b*c - ((27*(a^2)*d)+2b^3))//(54a^3)

    S = 0
    T = 0

    # cbrt(negative real numbers) evaluates normally with symbolics
    # while (im)^(1//3) also evaluates normally,
    # the other cases cbrt(im) and (-real)^(1/3) do not evaluate
    if eval(Symbolics.toexpr((Q^3+R^2))) > 0 
        S = Symbolics.term(cbrt, (R + Symbolics.term(sqrt, (Q^3+R^2))))
        T = Symbolics.term(cbrt, (R - Symbolics.term(sqrt, (Q^3+R^2))))
    else
        S = (Symbolics.term(complex, (R + im*Symbolics.term(sqrt, -(Q^3+R^2)))))^(1//3)
        T = (Symbolics.term(complex, (R - im*Symbolics.term(sqrt, -(Q^3+R^2)))))^(1//3)
    end


    root1 = S + T - (b//(3*a))
    root2 = -((S+T)//2) - (b//(3*a)) + (im*(Symbolics.term(sqrt, 3))/2)*(S-T)
    root3 = -((S+T)//2) - (b//(3*a)) - (im*(Symbolics.term(sqrt, 3))/2)*(S-T)

    return [root1, root2, root3]
end


function get_roots_deg4(expression, x)
    coeffs, constant = polynomial_coeffs(expression, [x])

    a = rationalize(coeffs[x^4])
    b = rationalize(get(coeffs, x^3, 0))
    c = rationalize(get(coeffs, x^2, 0))
    d = rationalize(get(coeffs, x, 0))
    e = rationalize(get(coeffs, x^0, 0))

    p = (8(a*c)-3(b^2))//(8(a^2))

    q = (b^3 - 4(a*b*c) + 8(d*a^2))//(8*a^3)

    r = (-3(b^4) + 256(e*a^3) - 64(d*b*a^2) + 16(c*(b^2)*a))//(256*a^4)

    @variables m y
    eq_m = 8m^3 + 8(p)*m^2 + (2(p^2) - 8r)m - q^2
    roots_m = solve(eq_m, m)
    m = 0
    for root in roots_m
        if !isequal(root, 0)
            m = root
            break
        end
    end

    arr = get_yroots(m, p, q)
    for (i, root) in enumerate(arr)
        arr[i] = root - (b//(4a))
    end

    return arr
end

function get_yroots(m, p, q)
    a = 1
    b1 = Symbolics.term(sqrt, Symbolics.term(complex, 2m))
    c1 = (p//2) + m - (q//(2*Symbolics.term(sqrt, Symbolics.term(complex, 2m))))
    b2 = -Symbolics.term(sqrt, Symbolics.term(complex, 2m))
    c2 = (p//2) + m + (q//(2*Symbolics.term(sqrt, Symbolics.term(complex, 2m))))

    root1 = simplify((-b1 + Symbolics.term(sqrt, (b1^2 - 4(a*c1)))) / 2a)
    root2 = simplify((-b1 - Symbolics.term(sqrt, (b1^2 - 4(a*c1)))) / 2a)
    if eval(Symbolics.toexpr(b1^2 - 4(a*c1))) < 0
        root1 = simplify((-b1 + Symbolics.term(sqrt, Symbolics.term(complex, (b1^2 - 4(a*c1))))) / 2a)
        root2 = simplify((-b1 - Symbolics.term(sqrt, Symbolics.term(complex, (b1^2 - 4(a*c1))))) / 2a)
    end


    root3 = simplify((-b2 + Symbolics.term(sqrt, (b2^2 - 4(a*c2)))) / 2a)
    root4 = simplify((-b2 - Symbolics.term(sqrt, (b2^2 - 4(a*c2)))) / 2a)
    if eval(Symbolics.toexpr(b2^2 - 4(a*c2))) < 0
        root3 = simplify((-b2 + Symbolics.term(sqrt, Symbolics.term(complex, (b2^2 - 4(a*c2))))) / 2a)
        root4 = simplify((-b2 - Symbolics.term(sqrt, Symbolics.term(complex, (b2^2 - 4(a*c2))))) / 2a)
    end

    return [root1, root2, root3, root4]
end


function get_roots(expression, x)
    degree = Symbolics.degree(expression, x)

    if degree == 0 && expression == 0
        return 0
    elseif degree == 0 && expression != 0
        throw("Not a valid statement")
    end

    if degree == 1
        # mx + c = 0
        coeffs, constant = polynomial_coeffs(expression, [x])
        m = rationalize(get(coeffs, x, 0))
        c = rationalize(get(coeffs, x^0, 0))
        root = -c // m
        return root
    end

    if degree == 2
        return get_roots_deg2(expression, x)
    end

    if degree == 3
        return get_roots_deg3(expression, x)
    end

    if degree == 4
        return get_roots_deg4(expression, x)
    end

end

