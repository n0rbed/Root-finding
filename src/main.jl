using Symbolics

coeff = Symbolics.coeff
function solve(expression, x)
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


function solve(eqs::Vector{Num}, vars::Vector{Num})
    eqs = Symbolics.groebner_basis(eqs)
    all_roots = Dict()

    # initialize the place of each var
    for var in vars
        all_roots[var] = nothing
    end


    #get roots for first var (z in this case)
    for eq in eqs
        for var in vars
            present_vars = Symbolics.get_variables(eq)
            if isequal(var, present_vars[1]) && size(present_vars, 1) == 1 
                var = Symbolics.get_variables(eq)[1]
                all_roots[var] = solve(Symbolics.wrap(eq), var)
            end
        end
    end

    
    # sub z, this is really messy and not systematic
    # as changing the number of vars will force me to
    # write more code, so i should think
    # of a better way to design this
    for (var, roots) in all_roots
        if !isnothing(roots)
            for eq in eqs
                if var in Symbolics.get_variables(eq)
                    substitute(eq, Dict([var => roots]))
                end
            end
        end
    end
end

@variables x y z
equations = [x*y + z - 40, x*z + y - 51, x + y + z - 19]

expr = [x^2 - 3*x - 10, x^2 + 20*x + 36]
println("x --> ", roots(expr, x))
