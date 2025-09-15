using HomotopyContinuation, Combinatorics, MixedSubdivisions, LinearAlgebra;
include("group_fns.jl");

# Here, n differs from the n in the paper by 1
n = 5;

# Set up the system 
@var L[1:n, 1:n, 1:n], x[1:n, 1:n], b[1:n, 1:n], d[1:n, 1:n];
x = Symmetric(x);
xvars = [x[i,j] for i in 1:n for j in i+1:n];

Ls = zeros(Expression, n, n, n);
Lvars = Vector{Variable}();
for (i,j,k) in Combinatorics.permutations(1:n, 3)
    Ls[i,j,k] = L[i,j,k]
    push!(Lvars, L[i,j,k])
end

eqns = [b[i,j]/(2*x[i,j]) - 1 - 1/(4*d[i,j]) * sum(Ls[i,j,k]*(2*x[j,k]^2 - x[i,j]^2)/(x[i,j]*x[j,k]*x[i,k]) + Ls[i,j,k]*(2*x[i,k]^2 - x[i,j]^2)/(x[i,j]*x[j,k]*x[i,k]) for k in 1:n) for (i,j) in Combinatorics.combinations(1:n, 2)];
F = System(eqns, parameters = [vec(b); vec(d); Lvars], variables = xvars);

# Compute mixed volume
target_b = ones(n^2); 
target_d = 2*ones(n^2); 
target_L = (1/n)*ones(length(Lvars));
F_eval = System(expand.(evaluate(F, xvars, [target_b; target_d; target_L])*prod(xvars.^2)), variables = xvars);
mv = mixed_volume(F_eval);

# Solve the system + certify solutions
mon_res = monodromy_solve(F, target_solutions_count=mv, threading=true);
certify(F, solutions(mon_res),  parameters(mon_res));
S = solve(F, solutions(mon_res); start_parameters = parameters(mon_res), target_parameters = [target_b; target_d; target_L]);
certify(F, solutions(S),  [target_b; target_d; target_L])

# Collect singular solutions without zero coordinates
pos_sols = [sol for sol in real.(solutions(S, only_real=true)) if all(sol .> 0)];
pos_sing_sols = [sol for sol in real.(solution.(singular(S, only_real=true))) if all(sol .> 0)];
pos_sols = [pos_sols; pos_sing_sols];
pos_sols = sort(pos_sols, by = prod);

# Define symmetric group action
function action(s, a)
    a*s*a'
end;

# Compute solution orbits under Sn action
Sn = SymmetricGroup(n);
matrix_sols = solution_to_matrix.(pos_sols, n); # reshape solutions into symmetric matrices
orbits = find_orbits(matrix_sols, defining_representation.(Sn,n), action);

# Check that we found all elements of each orbit
for key in keys(orbits)
    println(length(unique(map(a -> action(round.(key, digits=10), defining_representation(a, n)), Sn))) == length(orbits[key])+1)
end;

# Print representatives for each orbit and their volumes
display.(keys(orbits));
println(volume_of_metric.(matrix_to_solution.(keys(orbits))));