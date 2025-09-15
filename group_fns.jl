
function solution_to_matrix(v, n)
    t=1
    M = zeros(n, n)
    for i in 1:n
        for j in i+1:n
            M[i,j] = v[t]
            M[j,i] = v[t]
            t += 1
        end
    end
    M
end

function solution_to_matrix2(v, n)
    t=1
    M = zeros(n, n)
    for i in 1:n
        for j in 1:i-1
            M[i,j] = v[t]
            M[j,i] = v[t]
            t += 1
        end
    end
    M
end


function matrix_to_solution(M)
    [M[i,j] for i in 1:n for j in i+1:n]
end



function ij_representation(σ, n)
    M = zeros(binomial(n, 2), binomial(n, 2))
    for (t, (i,j)) in zip(1:binomial(n,2), Combinatorics.combinations(1:n, 2))
        for (r, (k,l)) in zip(1:binomial(n,2), Combinatorics.combinations(1:n, 2))
            M[t,r] = issetequal([σ[i], σ[j]], [k,l])
        end
    end
    M
end

function defining_representation(σ, n)
    [σ[i] == j for i in 1:n, j in 1:n]
end

function find_orbits(S, actors, action)
    D = Dict()
    for s in S
        found=false
        s_orbit = map(a -> action(s, a), actors)
        for key in keys(D)
            for s_rep in s_orbit
                if all(isapprox.(key, s_rep))
                    found=true
                    push!(D[key], s)
                    break
                end                
            end
        end
        if !found
            D[s] = []
        end
    end
    D
end

function volume_of_metric(v)
    prod(v.^2)
end

function only_nonzero(solution)
    return all(abs.(solution) .> 1e-10)
end

function latex_metric(s)
    s = round.(s, digits=5)
    print(s[1])
    for i in s[2:end]
        print(" & ")
        print(i)
    end
    print(" \\\\")
    println()
    print("\\hline")
    println()
end