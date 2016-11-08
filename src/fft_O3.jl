# This code generates the 2l+1 dimensional matrix Fourier coefficients
# for a function defined over the rotation and rotation/reflection groups
#
# The routines in this file are based upon:
#
#   FFTs on the Rotation Group
#   Peter J. Kostelec and Daniel N. Rockmore
#   J. Fourier. Anal. Appl., Vol. 14, Issue. 2, p. 145-179, 2008

import Base.fft

# Slow, reference version of the Fourier transform

function w(B::Int, k::Int)
    s = 0.0
    
    for j in 0:B-1
        s += sin((2j + 1) * (2k + 1) * (π / (4B))) / (2j + 1)
    end
    
    s *= (2/B) * sin(π * (2k + 1) / (4B))
end

function f̂(f::Function, B::Int, l::Int, M::Int, M′::Int)
    s = 0.0
    
    for j₁ in 0:2B-1, j₂ in 0:2B-1, k in 0:2B-1
        α = (2π * j₁) / (2B)
        β = π * (2k + 1) / (4B)
        γ = (2π * j₂) / (2B)
        
        s += w(B, k) * f(O3(α, β, γ)) * wigner_D̃(l, M, M′, α, β, γ)'
    end
    
    s /= (2B)^2
end

f̂_matrix(f::Function, B::Int, l::Int) = [f̂(f, B, l, M′, M) for M in -l:l, M′ in -l:l]

fourier_coeffs(f::Function, B::Int) = [f̂_matrix(f, B, l) for l in 0:B-1]

function T₂(f::Function, B::Int, k::Int, M::Int, M′::Int)
    s = 0.0
    
    for j₂ in 0:2B-1
        γ = (2π * j₂) / (2B)
        s += exp(im * M′ * γ) * exp(-im * ((2π * j₂) / (2B)) * (B)) * T₁(f, B, k, M, j₂)
    end
    
    s /= 2B
end

function fft_f̂(f::Function, B::Int, l::Int, M::Int, M′::Int)
    s = 0.0
    
    for k in 0:2B-1
        β = π * (2k + 1) / (4B)
        s += w(B, k) * wigner_d̃(l, M, M′, β) * T₂(f, B, k, M, M′)
    end
    
    return s
end

fft_f̂_matrix(f::Function, B::Int, l::Int) = [fft_f̂(f, B, l, M′, M) for M in -l:l, M′ in -l:l]

fft_fourier_coeffs(f::Function, B::Int) = [fft_f̂_matrix(f, B, l) for l in 0:B-1]

function T₁(f::Function, B::Int, k::Int, M::Int, j₂::Int)
    β = π * (2k + 1) / (4B)
    γ = (2π * j₂) / (2B)
    
    s = 0.0
    for j₁ in 0:2B-1
        α = (2π * j₁) / (2B)
        s += exp(im * M * α) * exp(-im * ((2π * j₁) / (2B)) * (B)) * f(O3(α, β, γ))
    end
    
    s /= 2B
end

# Fast Fourier transform

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
                f̂[l + 1][M′ + l + 1, M + l + 1] += wigner_d̃(l, M, M′, β) * w(B, k) * S[k + 1][M′ + B, M + B]
            end
        end
    end
    
    return f̂
end

fft(f::Function, ::Type{O3}, L::Int) = fast_fourier_coeffs(f, L)
