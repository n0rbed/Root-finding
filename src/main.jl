using Symbolics, Groebner, SymbolicUtils

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
                        push!(eqs, substitute(eq, Dict([var => root])))
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

@variables x y z
equations = [x*y + z - 40, x*z + y - 51, x + y + z - 19]

expr = [x^2 - 3*x - 10, x^2 + 20*x + 36]
println("equations to solve: ", equations)
println(solve(equations, [x, y , z]))
