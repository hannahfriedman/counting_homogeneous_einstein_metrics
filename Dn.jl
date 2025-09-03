using HomotopyContinuation, Combinatorics, MixedSubdivisions, LinearAlgebra;
include("group_fns.jl");

n = 4;

# Set up the system
@var x[1:n, 1:n], y[1:n, 1:n];
x = Symmetric(x);
y = Symmetric(y);
xvars = [x[i,j] for (i,j) in combinations(1:n, 2)];
yvars = [y[i,j] for (i,j) in combinations(1:n, 2)];
@var b[1:n, 1:n], d[1:n, 1:n], L[1:n, 1:n, 1:n], K[1:n, 1:n, 1:n];

eqns1 = [b[i,j]/(2*x[i,j]) - 1 + 1/(4*d[i,j]) * (sum([L[i,j,k]*(x[i,j]/(x[i,k]*x[k,j]) - x[i,k]/(x[i,j]*x[k,j]) - x[j,k]/(x[i,j]*x[i,k])) + K[i,j,k]*(x[i,j]/(y[i,k]*y[k,j]) - y[i,k]/(x[i,j]*y[k,j]) - y[j,k]/(x[i,j]*y[i,k])) for k in 1:n if (k != j && k != i)], init=0)) for (i,j) in combinations(1:n, 2)];
eqns2 = [b[i,j]/(2*y[i,j]) - 1 + 1/(4*d[i,j]) * (sum([L[k,j,i]*(y[i,j]/(x[i,k]*y[k,j]) - x[i,k]/(y[i,j]*y[k,j]) - y[j,k]/(y[i,j]*x[i,k])) + K[k,j,i]*(y[i,j]/(y[i,k]*x[k,j]) - y[i,k]/(y[i,j]*x[k,j]) - x[j,k]/(y[i,j]*y[i,k])) for k in 1:n if (k != j && k != i)], init=0)) for (i,j) in combinations(1:n, 2)];
eqns = [eqns1; eqns2];
F = System(eqns, variables=[xvars; yvars], parameters=[vec(L); vec(K); vec(b); vec(d)]);

# Compute mixed volume
F_eval = expand.(evaluate(F, [xvars; yvars], [(1/(2*n-2))*ones(2*length(vec(L))); ones(length(vec(b))); ones(length(vec(d)))]) * prod(xvars) * prod(yvars));
mv = MixedSubdivisions.mixed_volume(System(F_eval));

# Solve the system + certify solutions
mon_res = monodromy_solve(F, target_solutions_count=mv);
certify(F, solutions(mon_res), parameters(mon_res));
S = solve(F, solutions(mon_res); start_parameters = parameters(mon_res), target_parameters = [(1/(2*n-2))*ones(2*length(vec(L))); ones(length(vec(b))); ones(length(vec(d)))]);
certify(F, solutions(S), [(1/(2*n-2))*ones(2*length(vec(L))); ones(length(vec(b))); ones(length(vec(d)))]);

# Collect singular solutions without zero coordinates
pos_sols = [sol for sol in real.(solutions(S, only_real=true)) if all(sol .> 0)];
pos_sing_sols = [sol for sol in real.(solution.(singular(S, only_real=true))) if all(sol .> 0)];
pos_sols = [pos_sols; pos_sing_sols];
pos_sols = sort(pos_sols, by = prod);

# Analyze Weyl group action on solutions
function action((x, y), (σ, τ))
    new_x = []
    new_y = []
    for ((i,j), t) in zip(combinations(1:n, 2), 1:6)
        if τ[i]*τ[j] == 1
            push!(new_x, x[t])
            push!(new_y, y[t])
        else
            push!(new_x, y[t])
            push!(new_y, x[t])
        end
    end
    σ_mat = ij_representation(σ, n)

    return (σ_mat * new_x, σ_mat * new_y)
end;

# Compute solution orbits under Sn action
Ω = [(sol[1:binomial(n,2)], sol[binomial(n,2)+1:2*binomial(n,2)]) for sol in pos_sols];
actors = vec(collect(Base.product(SymmetricGroup(n), [s for s in Combinatorics.multiset_permutations([ones(n); -ones(n)], n) if prod(s) == 1])));
orbits = find_orbits(Ω, actors, action);

# Check that we found all elements of each orbit
for key in keys(orbits)
    println(length(unique(map(a -> action((round.(key[1], digits=10), round.(key[2], digits=10)), a), actors))) == length(orbits[key]) + 1);
end

# Print representatives for each orbit and their volumes
display.(keys(orbits));
println(volume_of_metric.([key[1]; key[2]] for key in keys(orbits)));