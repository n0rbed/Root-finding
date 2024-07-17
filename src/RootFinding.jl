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

# polynomial solvers
export solve
export solve_univar
export solve_multipoly
export solve_multivar

# univar solvers (analytic solutions)
export get_roots
export get_roots_deg1
export get_roots_deg2
export get_roots_deg3
export get_roots_deg4

# polynomial solver helpers
export filter_poly
export lead_term
export lead_coeff
export gcd_use_nemo
export check_polynomial
export postprocess_root

# transcendental solver
export nl_solve

# trans solver helpers
export n_occurrences
export n_func_occ
export attract
export isolate
export turn_to_poly

# tools
export ssqrt
export slog
export scbrt
export ssubs


end
