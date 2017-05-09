function ft_symmetric_group(f::Function, n::Int)
    table = standard_tableaux(n)
    d_λ = Dict(λ => length(table[n][λ]) for λ in integer_partitions(n))
    f̂ = Dict(λ => zeros(Float64, d_λ[λ], d_λ[λ]) for λ in integer_partitions(n))
            
    for σ in permutations(1:n, n)
        g = SymmetricGroup(σ)
        for λ in integer_partitions(n)
            f̂[λ] .+= f(g) * irrep(g, λ)'
        end
    end

    for λ in integer_partitions(n)
        f̂[λ] *= d_λ[λ] / factorial(n)
    end
            
    return f̂
end

fft{N}(f::Function, ::Type{SymmetricGroup{N}}, L::Int) = ft_symmetric_group(f, N)
