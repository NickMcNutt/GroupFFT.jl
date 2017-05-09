# This code generates the 2l+1 dimensional matrix Fourier coefficients
# for a function defined over the rotation and rotation/reflection groups
#
# The routines in this file are based upon:
#
#   FFTs on the Rotation Group
#   Peter J. Kostelec and Daniel N. Rockmore
#   J. Fourier. Anal. Appl., Vol. 14, Issue. 2, p. 145-179, 2008

import Base.fft

function w(B::Int, k::Int)
    s = 0.0
    
    for j in 0:B-1
        s += sin((2j + 1) * (2k + 1) * (π / (4B))) / (2j + 1)
    end
    
    s *= (2/B) * sin(π * (2k + 1) / (4B))
end

function S₁(f::Function, B::Int, k::Int, j₂::Int)
    β = π * (2k + 1) / (4B)
    γ = (2π * j₂) / (2B)
    
    A = [exp(-im * ((2π * j₁) / (2B)) * (B - 1)) * f(O3((2π * j₁) / (2B), β, γ)) for j₁ in 0:2B-1]
    A′ = ifft(A)
end

function S₂(f::Function, B::Int, k::Int)
    S = [S₁(f, B, k, j₂) for j₂ in 0:2B-1]
    A = Vector{Vector{Complex128}}(2B)
    A′ = Vector{Vector{Complex128}}(2B)
    
    for M′ in -(B-1):B
        A[M′ + B] = [exp(-im * ((2π * j₂) / (2B)) * (B - 1)) * S[j₂ + 1][M′ + B] for j₂ in 0:2B-1]
        A′[M′ + B] = ifft(A[M′ + B])
    end
    
    return hcat(A′...)
end

function fast_fourier_coeffs(f::Function, B::Int)
    S = [S₂(f, B, k) for k in 0:2B-1]
    f̂ = [zeros(Complex128, 2l + 1, 2l + 1) for l in 0:B-1]
    
    for l in 0:B-1
        for M in -l:l, M′ in -l:l
            for k in 0:2B-1
                β = π * (2k + 1) / (4B)
                f̂[l + 1][M′ + l + 1, M + l + 1] += wigner_d̃(l, M, M′, β) * w(B, k) * S[k + 1][M′ + B, M + B] / sqrt(2 / (2l + 1))
            end
        end
    end
    
    return f̂
end

fft(f::Function, ::Type{O3}, L::Int) = fast_fourier_coeffs(f, L)
