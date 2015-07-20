using JuMP
using Distributions
Pkg.add("Clp")
using Clp

y = [0.0101110, 0.0043532, 0.0137058] #vector of returns
cov = [0.00324625 0.00022983 0.00420395; 0.00022983 0.00049937 0.00019247; 0.00420395 0.00019247 0.00764097]
R = 0.011 #minimum acceptable portfolio return
k = 20^4
smpl = rand(MvNormal(y, cov), k)
β = 0.99

m = Model(solver = ClpSolver())

@defVar(m, 0 <= x[1:3] <= 1)
@defVar(m, α)
@defVar(m, u[1:k] >= 0)

@setObjective(m, Min, α + 1/(k*(1-β))*sum{u[i], i = 1:k} )
@addConstraint(m, sum{x[i], i = 1:3} == 1)
@addConstraint(m, sum{y[i]*x[i], i = 1:3} >= R)
@defConstrRef myCons[1:k]
for i in 1:k
  myCons[i] = @addConstraint(m, u[i] + α + sum{smpl[j, i]*x[j], j = 1:3} >= 0)
end

print(m)

status = solve(m)

println("x = ", getValue(x))
println("α = ", getValue(α))
println("Objective value: ", getObjectiveValue(m))
