import SymbolicUtils
import Symbolics
import TermInterface: maketerm

using Test

# Alex: can be used in the following way.
#
# sol = RootFinding.solve(x^12 - 1, x)
# map(RootFinding.postprocess_root, sol)

# Alex: make sure `Num`s are not processed here as they'd break it.
_is_const_number(x::Number) = true
function _is_const_number(x::SymbolicUtils.BasicSymbolic)
    !istree(x) && return false
    all(_is_const_number, arguments(x))
end

SymbolicUtils.@syms __x
@test !_is_const_number(__x) && !_is_const_number(sqrt(__x))
@test _is_const_number(1) && _is_const_number(2 // 3) && _is_const_number(3 + 4im)
@test _is_const_number(SymbolicUtils.term(sqrt, 2) + 21)
@test _is_const_number((SymbolicUtils.term(exp, 2) * SymbolicUtils.term(exp, 2)) // 99)

function _postprocess_root(x::Number)
    # N // 1 => N
    if x isa Rational
        if isone(denominator(x))
            @info x => big(numerator(x))
            return big(numerator(x))
        end
    end

    # A + im*0 => A
    if x isa Complex
        if iszero(imag(x))
            # Alex: maybe return _postprocess_root(real(x))
            @info x => real(x)
            return real(x)
        end
    end
    
    return x
end

function _postprocess_root(x::SymbolicUtils.BasicSymbolic)
    !istree(x) && return x

    x = maketerm(
        typeof(x), 
        operation(x),
        map(_postprocess_root, arguments(x)),
        nothing
    )

    # sqrt(0), cbrt(0) => 0
    # sqrt(1), cbrt(1) => 1
    if istree(x) && (operation(x) === sqrt || operation(x) === cbrt)
        arg = arguments(x)[1]
        if isequal(arg, 0) || isequal(arg, 1)
            @info x => arg
            return arg
        end
    end

    # sqrt(N^2) => N
    if istree(x) && operation(x) === sqrt
        arg = arguments(x)[1]
        if arg isa Integer && arg == (isqrt(arg))^2
            @info x => isqrt(arg)
            return isqrt(arg)
        end
    end

    # (sqrt(N))^2 => N
    if istree(x) && operation(x) === (^) && isequal(arguments(x)[2], 2)
        arg1 = arguments(x)[1]
        if istree(arg1) && operation(arg1) === sqrt
            if _is_const_number(arguments(arg1)[1])
                @info x => arguments(arg1)[1]
                return arguments(arg1)[1]
            end
        end
    end

    return x
end

function postprocess_root(x)
    _postprocess_root(x)
end

SymbolicUtils.@syms __x
__symsqrt(x) = SymbolicUtils.term(sqrt, x)
@test postprocess_root(2 // 1) == 2 && postprocess_root(2 + 0*im) == 2
@test postprocess_root(__symsqrt(__symsqrt(0)) - 11) == -11
@test postprocess_root(3*__symsqrt(2)^2) == 6
@test postprocess_root(__symsqrt(4)) == 2
@test isequal(postprocess_root(__symsqrt(__x)^2), __symsqrt(__x)^2)