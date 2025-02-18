# Running all scripts for creating synthetic time series for all model systems

include("generate_data_fw5.jl")
include("generate_data_switching.jl")
include("generate_data_logistic.jl")
include("generate_jmat_discretised_fw5.jl")
include("generate_jmat_discretised_switching.jl")