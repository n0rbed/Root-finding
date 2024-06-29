using Symbolics

function comp_rational(x,y)
    try
        r = x//y
        return r
    catch e
        r = nothing
        if x isa Float64 || x isa ComplexF64
            real_p = real(x)
            imag_p = imag(x)
            r = Rational{BigInt}(real_p)//y
            r += (Rational{BigInt}(imag_p)//y)*im
        end
        return r
    end
end

function get_roots_deg1(expression, x)
    subs, expr = filter_poly(expression, x)
    coeffs, constant = polynomial_coeffs(simplify(real(expr)), [x])
    m = get(coeffs, x, 0)
    c = get(coeffs, x^0, 0)
    root = -comp_rational(c,m)
    for (var, sub) in subs
        root = Symbolics.substitute(root, Dict([var => sub]), fold=false)
    end
    return [root]
end
function get_roots_deg2(expression, x)
    # ax^2 + bx + c = 0
    subs, filtered_expr = filter_poly(expression, x)
    coeffs, constant = polynomial_coeffs(filtered_expr, [x])

    results = (substitute(get(coeffs, x^i, 0), subs, fold=false) for i in 2:-1:0)
    a, b, c = results

    root1 = comp_rational(-b + Symbolics.term(sqrt, comp_rational((b^2 - 4(a*c)), 1)), 2a)
    root2 = comp_rational(-b - Symbolics.term(sqrt, comp_rational((b^2 - 4(a*c)), 1)), 2a)
    try
        if eval(Symbolics.toexpr(b^2 - 4(a*c))) < 0
            root1 = simplify((-b + Symbolics.term(sqrt, Symbolics.term(complex, (b^2 - 4(a*c))))) // 2a)
            root2 = simplify((-b - Symbolics.term(sqrt, Symbolics.term(complex, (b^2 - 4(a*c))))) // 2a)
        end
    catch e
        if typeof(e) == UndefVarError
            println("Roots: ", root1, " and ", root2, " assume that ", b^2 - 4(a*c), " > 0")
        end
    end

    return [root1, root2]
end

function get_roots_deg3(expression, x)
    subs, filtered_expr = filter_poly(expression, x)
    coeffs, constant = polynomial_coeffs(filtered_expr, [x])

    results = (substitute(get(coeffs, x^i, 0), subs, fold=false) for i in 3:-1:0)
    a, b, c, d = results

    
    Q = comp_rational((((3*a*c) - b^2)), (9a^2))
    R = comp_rational(((9*a*b*c - ((27*(a^2)*d)+2b^3))), (54a^3))


    # cbrt(negative real numbers) evaluates normally with symbolics
    # while (im)^(1//3) also evaluates normally,
    # the other cases cbrt(im) and (-real)^(1/3) do not evaluate
    S = Symbolics.term(cbrt, (R + Symbolics.term(sqrt, (Q^3+R^2))))
    T = Symbolics.term(cbrt, (R - Symbolics.term(sqrt, (Q^3+R^2))))
    try 
        eval(Symbolics.toexpr(S))
        eval(Symbolics.toexpr((Q^3+R^2))) < 0 
    catch e
        if typeof(e) == UndefVarError
            println("Condition: ", Q^3+R^2, " > 0")
        else
            S = (Symbolics.term(complex, (R + im*Symbolics.term(sqrt, -(Q^3+R^2)))))^(1//3)
            T = (Symbolics.term(complex, (R - im*Symbolics.term(sqrt, -(Q^3+R^2)))))^(1//3)
        end
    end

    root1 = S + T - (b//(3*a))
    root2 = -((S+T)//2) - (b//(3*a)) + (im*(Symbolics.term(sqrt, 3))/2)*(S-T)
    root3 = -((S+T)//2) - (b//(3*a)) - (im*(Symbolics.term(sqrt, 3))/2)*(S-T)

    return [root1, root2, root3]
end



function get_roots_deg4(expression, x)
    subs, filtered_expr = filter_poly(expression, x)
    coeffs, constant = polynomial_coeffs(filtered_expr, [x])

    results = (substitute(get(coeffs, x^i, 0), subs, fold=false) for i in 4:-1:0)
    a, b, c, d, e = results

    p = comp_rational((8(a*c)-3(b^2)), (8(a^2)))

    q = comp_rational(b^3 - 4(a*b*c) + 8(d*a^2), (8*a^3))

    r = comp_rational((-3(b^4) + 256(e*a^3) - 64(d*b*a^2) + 16(c*(b^2)*a)), (256*a^4))

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
        r = comp_rational(b, 4a)
        arr[i] = root - (r)
    end

    return arr
end

function get_yroots(m, p, q)
    a = 1
    if m isa Real && (eval(Symbolics.toexpr(m))) < 0
        b1 = im*Symbolics.term(sqrt,  -2m)
        c1 = (p//2) + m - (q//(2*im*Symbolics.term(sqrt, -2m)))
        b2 = -im*Symbolics.term(sqrt, -2m)
        c2 = (p//2) + m + (q//(2*im*Symbolics.term(sqrt, -2m)))
    else
        b1 = Symbolics.term(sqrt,  2m)
        c1 = (p//2) + m - (q//(2*Symbolics.term(sqrt, 2m)))
        b2 = -Symbolics.term(sqrt, 2m)
        c2 = (p//2) + m + (q//(2*Symbolics.term(sqrt, 2m)))
    end

    root1 = simplify((-b1 + Symbolics.term(sqrt, (b1^2 - 4(a*c1)))) / 2a)
    root2 = simplify((-b1 - Symbolics.term(sqrt, (b1^2 - 4(a*c1)))) / 2a)
    if real(eval(Symbolics.toexpr(b1^2 - 4(a*c1)))) < 0
        root1 = simplify((-b1 + im*Symbolics.term(sqrt, -(b1^2 - 4(a*c1)))) / 2a)
        root2 = simplify((-b1 - im*Symbolics.term(sqrt, -(b1^2 - 4(a*c1)))) / 2a)
    end


    root3 = simplify((-b2 + Symbolics.term(sqrt, (b2^2 - 4(a*c2)))) / 2a)
    root4 = simplify((-b2 - Symbolics.term(sqrt, (b2^2 - 4(a*c2)))) / 2a)
    if real(eval(Symbolics.toexpr(b2^2 - 4(a*c2)))) < 0
        root3 = simplify((-b2 + im*Symbolics.term(sqrt, -(b2^2 - 4(a*c2)))) / 2a)
        root4 = simplify((-b2 - im*Symbolics.term(sqrt, -(b2^2 - 4(a*c2)))) / 2a)
    end

    return [root1, root2, root3, root4]
end


function get_roots(expression, x)
    subs, filtered_expr = filter_poly(expression, x)
    degree = Symbolics.degree(filtered_expr, x)

    if degree == 0 && expression == 0
        return []
    elseif degree == 0 && expression != 0
        throw("Not a valid statement")
    end



    if degree == 1
        # mx + c = 0
        return get_roots_deg1(expression, x)
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
