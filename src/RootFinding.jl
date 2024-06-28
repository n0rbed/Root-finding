module RootFinding

include("nemo_stuff.jl")
include("nonlinear.jl")
include("main.jl")
export solve
export get_roots
export filter_poly
export lead_term
export lead_coeff
export gcd_use_nemo
export nl_solve
export n_occurrences
export split_by_variable


end
