using JuMP, Gurobi # 2 bibliothèques qui contiennent des solveurs pour pb d'optimisation

model = Model(Gurobi.Optimizer)
set_optimizer_attribute(model, "TimeLimit", 18000)
set_optimizer_attribute(model, "MIPGap", 0.000001)

# --- Données ---
m, n = 10, 10
p_rares, q_communes = 3, 3
alphaR, alphaC = 0.5, 0.9

raw_proba = [
    1 4 3 0.4; 1 5 3 0.3; 1 6 2 0.4; 1 8 6 0.3; 1 8 7 0.2; 1 9 5 0.2; 1 9 6 0.4;
    2 2 7 0.2; 2 3 7 0.4; 2 4 7 0.2; 2 5 9 0.4; 2 5 10 0.3; 2 7 2 0.5; 2 9 6 0.2; 2 9 7 0.2;
    3 2 4 0.2; 3 2 5 0.3; 3 2 7 0.4; 3 3 7 0.4; 3 4 7 0.5; 3 5 9 0.2; 3 5 10 0.4; 3 9 2 0.3;
    4 1 2 0.3; 4 1 4 0.3; 4 3 2 0.4; 4 4 3 0.4; 4 5 1 0.3; 4 5 3 0.2; 4 5 5 0.2; 4 6 2 0.2; 
    4 6 3 0.4; 4 7 5 0.4; 4 7 9 0.3; 4 8 7 0.2; 4 8 9 0.5; 4 9 7 0.4;
    5 1 10 0.4; 5 2 1 0.3; 5 2 10 0.3; 5 3 5 0.5; 5 6 3 0.4; 5 6 6 0.2; 5 6 7 0.2; 
    5 7 5 0.4; 5 7 9 0.5; 5 8 9 0.4; 5 9 2 0.5; 5 9 3 0.2; 5 9 4 0.4; 5 9 5 0.4;
    6 1 3 0.4; 6 1 4 0.4; 6 1 6 0.5; 6 1 7 0.3; 6 1 8 0.3; 6 2 9 0.2; 6 3 5 0.4; 
    6 3 9 0.4; 6 3 10 0.2; 6 4 9 0.3; 6 5 5 0.4; 6 8 1 0.5; 6 9 1 0.2; 6 10 3 0.2
]

# Initialisation du tenseur 3D avec des zéros (le comportement "default 0" de AMPL)
proba = zeros(Float64, p_rares + q_communes, m, n)

# Remplissage du tenseur proba[k, i, j]
for row in 1:size(raw_proba, 1)
    k = Int(raw_proba[row, 1])
    i = Int(raw_proba[row, 2])
    j = Int(raw_proba[row, 3])
    val = raw_proba[row, 4]
    
    proba[k, i, j] = val
end


ER = 1:p_rares
EC = (p_rares+1):(p_rares+q_communes)

# Initialisation des coûts selon la logique AMPL
c_cost = zeros(Int, m, n)
for i in 1:m, j in 1:n
    if i<=3 && j<=3;     c_cost[i,j] = 6
    elseif i<=3 && j<=7; c_cost[i,j] = 4
    elseif i<=3 && j<=10; c_cost[i,j] = 8
    elseif i<=7 && j<=3; c_cost[i,j] = 5
    elseif i<=7 && j<=7; c_cost[i,j] = 3
    elseif i<=7 && j<=10; c_cost[i,j] = 7
    elseif i<=10 && j<=3; c_cost[i,j] = 4
    elseif i<=10 && j<=7; c_cost[i,j] = 6
    else;                 c_cost[i,j] = 5
    end
end

println("Exécution OK")

@variable(model, x[1:m, 1:n], Bin)
@variable(model, y[1:m, 1:n], Bin)

@objective(model, z, sum())

@constraint()
@constraint()


optimize!(model)

println("valeur du pb : ", objective_value(model))
println("valeurs de x" : value.(x))
println("valeurs de y : ", value.(y))
