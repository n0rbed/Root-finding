module RootFinding

include("coeffs.jl")
include("nemo_stuff.jl")
include("solve_helpers.jl")
include("postprocess.jl")
include("univar.jl")
include("isoa_helpers.jl")
include("polynomialization.jl")
include("attract.jl")
include("isoa_main.jl")
include("main.jl")

export solve

# helper solvers
export solve_univar
export solve_multipoly
export solve_multivar
export ia_solve

end
