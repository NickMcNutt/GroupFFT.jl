module GroupFFT

using Combinatorics, Groups

export
    fft

include("fft_O3.jl")
include("fft_symmetric_group.jl")

end
