module Bretthorst

export schuster, student, brett1d, brett2d

using LinearAlgebra

BLAS.set_num_threads(Threads.nthreads())

include("schuster.jl")
include("student.jl")
include("bretthorst1d.jl")
include("bretthorst2d.jl")
end
