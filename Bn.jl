using HomotopyContinuation, Combinatorics, MixedSubdivisions, LinearAlgebra;
include("group_fns.jl");

n = 3;

# Set up the system
@var x[1:n, 1:n], y[1:n, 1:n], z[1:n];
x = Symmetric(x);
y = Symmetric(y);
xvars = [x[i,j] for (i,j) in combinations(1:n, 2)];
yvars = [y[i,j] for (i,j) in combinations(1:n, 2)];
@var b[1:n, 1:n], d[1:n, 1:n], L[1:1000];

eqns1 = [b[i,j]/(2*x[i,j]) - 1 + 1/(4*d[i,j]) * (rand(L)*(x[i,j]/(z[i]*z[j]) - z[i]/(x[i,j]*z[j]) - z[j]/(x[i,j]*z[i])) + sum([rand(L)*x[i,j]/(x[i,k]*x[k,j]) - rand(L)*x[i,k]/(x[i,j]*x[k,j]) - rand(L)*x[j,k]/(x[i,j]*x[i,k]) + rand(L)*x[i,j]/(y[i,k]*y[k,j]) - rand(L)*y[i,k]/(x[i,j]*y[k,j]) - rand(L)*y[j,k]/(x[i,j]*y[i,k]) for k in 1:n if (k != j && k != i)], init=0)) for (i,j) in combinations(1:n, 2)];
eqns2 = [b[i,j]/(2*y[i,j]) - 1 + 1/(4*d[i,j]) * (rand(L)*(y[i,j]/(z[i]*z[j]) - z[i]/(y[i,j]*z[j]) - z[j]/(y[i,j]*z[i])) + sum([rand(L)*y[i,j]/(x[i,k]*y[k,j]) - rand(L)*x[i,k]/(y[i,j]*y[k,j]) - rand(L)*y[j,k]/(y[i,j]*x[i,k]) + rand(L)*y[i,j]/(x[j,k]*y[k,i]) - rand(L)*x[j,k]/(y[i,j]*y[k,i]) - rand(L)*y[i,k]/(y[i,j]*x[j,k]) for k in 1:n if (k != j && k != i)], init=0)) for (i,j) in combinations(1:n, 2)];
eqns3 = [b[i,i]/(2*z[i]) - 1 + 1/(4*d[i,i]) * (sum(rand(L)*(z[i]/(z[k]*x[k,i]) - x[k,i]/(z[i]*z[k]) - z[k]/(z[i]*x[k,i])) + rand(L)*(z[i]/(z[k]*y[k,i]) - y[k,i]/(z[i]*z[k]) - z[k]/(z[i]*y[k,i])) for k in 1:n if k != i)) for i in 1:n];
eqns = [eqns1; eqns2; eqns3];
F = System(eqns, variables=[xvars; yvars; z], parameters=[L; vec(b); vec(d)]);

# Compute mixed volume
F_eval = System(expand.(evaluate(F, [xvars; yvars; z], [(1/(2*n-1))*ones(length(vec(L))); ones(length(vec(b))); ones(length(vec(d)))]) * prod(xvars) * prod(yvars) * prod(z)));
mv = mixed_volume(F_eval);

# Solve the system + certify solutions
mon_res = monodromy_solve(F, target_solutions_count=mv);
certify(F, solutions(mon_res), parameters(mon_res));
S = solve(F, solutions(mon_res); start_parameters = parameters(mon_res), target_parameters = [(1/(2*n-1))*ones(length(vec(L))); ones(length(vec(b))); ones(length(vec(d)))]);
certify(F, solutions(S),  [(1/(2*n-1))*ones(length(vec(L))); ones(length(vec(b))); ones(length(vec(d)))]);

# Collect singular solutions without zero coordinates
println(setdiff(filter(only_nonzero, solutions(S, only_nonsingular=false)), solutions(S)));

# Analyze Weyl group action on solutions
pos_sols = [sol for sol in real.(solutions(S, only_real=true)) if all(sol .> 0)];
pos_sing_sols = [sol for sol in real.(solution.(singular(S, only_real=true))) if all(sol .> 0)];
pos_sols = [pos_sols; pos_sing_sols];
pos_sols = sort(pos_sols, by = volume_of_metric);

# Define symmetric group action
function action((x, y, z), (σ, τ))
    n = 3
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

    return (σ_mat * new_x, σ_mat * new_y, defining_representation(σ, n) * z)
end;

# Compute solution orbits under Sn action
Ω = [(sol[1:3], sol[4:6], sol[7:end]) for sol in pos_sols];
actors = vec(collect(Base.product(SymmetricGroup(n), Combinatorics.multiset_permutations([1,1,1,-1,-1,-1], n))));
orbits = find_orbits(Ω, actors, action);

# Check that we found all elements of each orbit
for key in keys(orbits)
    println(length(unique(map(a -> action((round.(key[1], digits=11), round.(key[2], digits=11), round.(key[3], digits=11)), a), actors))) == length(orbits[key]) + 1)
end

# Print representatives for each orbit and their volumes
display.(keys(orbits));
println(volume_of_metric.([key[1]; key[2]; key[3]] for key in keys(orbits)));