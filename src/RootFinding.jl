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
export get_roots
export get_roots_deg1
export get_roots_deg2
export get_roots_deg3
export get_roots_deg4
export filter_poly
export lead_term
export lead_coeff
export gcd_use_nemo
export check_polynomial
export nl_solve
export n_occurrences
export n_func_occ
export split_by_variable
export postprocess_root
export attract
export isolate


end
