using Gurobi, JuMP
using Pkg
Pkg.add("Plots")
using Plots

println("Hello")

"""
Les constantes
"""
m = 10
n = 10

B = 3500
Amin = 70
Amax = 75

d = [sqrt((i-k)^2 + (j-l)^2) for i=1:m, j=1:n, k=1:m, l=1:n];
c = [
 7 3 10 10 2 8 6 4 5 5;
 7 7 10 5 2 8 6 3 9 9;
 7 3 4 6 3 2 4 9 7 8;
 6 2 7 6 4 7 5 10 7 8;
 2 4 3 4 9 6 4 9 8 4;
 7 5 2 9 8 9 5 6 10 10;
 5 2 3 7 9 9 4 9 6 3;
 5 2 9 4 2 8 6 9 3 4;
 9 6 5 4 5 6 8 9 6 6;
 8 8 7 7 3 5 8 3 9 9;
].*10 ;

# 1 point = 2 coordonnées

"""
Les fonctions
"""
function f(ppv)
    return sum(d[i,j,k,l]*ppv[i,j,k,l] for i in 1:m, j in 1:n, k in 1:m, l in 1:n)
end

function g(selection)
    return sum(selection[i,j] for i in 1:m, j in 1:n)
end

"""
Algo de Binkelbach
"""

# 1. prendre un lambda quelconque : 
tolerance = 1e-5
lambda = 10
i = 1
temps_calcul = 0
nb_noeuds = 0

while (true)
    global lambda # dans une boucle, Julia considère qu'on modifie les variables qu'en local, et non hors de la boucle
    global valeur_P_lambda
    global selection_opti
    global ppv_opti
    global i
    global temps_calcul
    global nb_noeuds

    println(i)

    # 2. Calculer la valeur du problème P_λ :
    model = Model(Gurobi.Optimizer);
    set_optimizer_attribute(model, "OutputFlag", 0);

    # variables
    @variable(model, selection[1:m, 1:n], Bin);
    @variable(model, ppv[1:m, 1:n, 1:m, 1:n], Bin);

    # contraintes spatiales & budgétaire
    @constraint(model, aire_max, sum(selection[i,j] for i in 1:m, j in 1:n) <= Amax);
    @constraint(model, aire_min, sum(selection[i,j] for i in 1:m, j in 1:n) >= Amin);
    @constraint(model, cout, sum(c[i,j] * selection[i,j] for i in 1:m, j in 1:n) <= B);

    # contraintes inhérentes à la définition du ppv
    @constraint(model, [i in 1:m, j in 1:n], ppv[i,j,i,j] == 0);
    @constraint(model,[i in 1:m, j in 1:n], sum(ppv[i,j,k,l] for k in 1:m, l in 1:n)==selection[i,j]);
    @constraint(model,[i in 1:m, j in 1:n, k in 1:m, l in 1:n], ppv[i,j,k,l]<=selection[k,l]);

    # objectif
    @objective(model, Min, f(ppv)-lambda*g(selection));

    optimize!(model);
    temps_calcul += solve_time(model)
    nb_noeuds += MOI.get(model, MOI.NodeCount())

    # 3. Récupérer le x optimal :
    valeur_P_lambda = objective_value(model);
    selection_opti = value.(selection);
    ppv_opti = value.(ppv);

    if valeur_P_lambda >= - tolerance
        break                           
    else
        i+=1
        lambda = f(ppv_opti)/g(selection_opti);
    end
end

#if valeur_P_lambda>
println("valeur du problème = ", f(ppv_opti)/g(selection_opti))
println("nombre  d'itérations Dinkelbach = ", i)
println("nombre de noeuds explorés : ", nb_noeuds)
println("temps de calcul total =", temps_calcul)
println(selection_opti)
savefig(plot(heatmap(selection_opti)),  "heatmap.png")
#nothing # sinon Julia a tendance à print la dernière ligne

