using Gurobi;
using JuMP;

N=8; # nombre d'individus dans la population
Nm=4; # nombre d'individus males
Nf=4; # nombre d'individus femelles
C=1; # nombre de paires de chromosome par individu (je crois qu'on s'en fiche comme on a supposé C=1 dans notre mod)
G=5; # nombre de locus, = nombre de gènes je crois
A=2; # nombre d'allèles par gène, il s'avère qu'il y a le même nombre d'allèles pour tous les gènes

# hyperparamètres des coupures theta_r :
T=50; 
init=0.001; 

individus = # mon n_agi
[
[[[2,1],   # les deux allèles du gène 1 de l'individu 1
[1,1],     # les deux allèles du gène 2 de l'individu 1
[2,1],     # ...
[2,1],
[1,2],],],

[[[2,2],   # les deux allèles du gène 1 de l'individu 2
[2,1],     # ...
[1,2],
[1,2],
[1,1],],],

[[[2,2],   # individu 3
[2,1],
[2,1],
[1,2],
[1,2],],],


[[[2,2],  # individu 4
[1,1],
[1,2],
[2,2],
[2,2],],],

[[[2,1],  # individus 5
[1,1],
[1,1],
[2,2],
[1,2],],],

[[[1,2], #individus 6
[1,1],
[2,2],
[2,1],
[2,1],],],

[[[1,1], # individu 7
[1,1],
[2,2],
[1,1],
[1,1],],],


[[[2,2], # individu 8
[1,1],
[2,2],
[2,1],
[2,1],],],
];

# compteur[i,g,a] = nombre de copies (0, 1 ou 2) de l'allèle a au locus g chez l'individu i
compteur = zeros(Int, N, G, A)
for i in 1:N, g in 1:G, k in 1:2
    compteur[i, g, individus[i][1][g][k]] += 1   # [1] car C = 1
end

function theta(r, T=T, init=init)
    return init^((T-r)/(T-1))
end

function proba_disparition_allele(a,g,n_enfants)
    resultat = 1
    for i in 1:N
        resultat *= (1 - compteur[i, g, a]/2)^n_enfants[i]
    end
    return resultat
end

x_max = 2
model = Model(Gurobi.Optimizer)

@variable(model, n_enfants[1:N]) # enlever le Int ici pour la relaxation
@variable(model, y[1:G, 1:A]>=0)
@variable(model, t[1:G, 1:A])

@constraint(model, conservation_m, sum(n_enfants[i] for i in 1:Nm) == N)
@constraint(model, conservation_f, sum(n_enfants[i] for i in Nm+1:Nm+Nf) == N)
@constraint(model,[i in 1:N], 0 <= n_enfants[i] <= x_max)
@constraint(model, [g in 1:G, a in 1:A], y[g,a]>=t[g,a]- sum(n_enfants[i] for i in 1:N if compteur[i, g, a] == 2))
@constraint(model, [g in 1:G, a in 1:A, r in 1:T], log(theta(r)) + 1/(theta(r)) * (t[g,a] - theta(r)) >= - log(2)*sum(n_enfants[i] for i in 1:N if compteur[i, g, a] == 1))
@objective(model, Min, sum((sum(y[g,a] for a in 1:A)) for g in 1:G))

optimize!(model);

temps_calcul = solve_time(model);
nb_noeuds = MOI.get(model, MOI.NodeCount());
valeur_P = objective_value(model); # borne inf : relancer avec la variable Float et non plus Int
n_enfants_opti = value.(n_enfants);
probas = [[proba_disparition_allele(a,g,n_enfants_opti) for a in 1:A] for g in 1:G];

println("valeur", valeur_P)
println("temps", temps_calcul)
println("nb noeuds", nb_noeuds)
println("probas", probas)