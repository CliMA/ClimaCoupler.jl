import SciMLBase: step!
using OrdinaryDiffEq
using OrdinaryDiffEq: ODEProblem, SSPRK33, solve, savevalues!, Euler
using Printf
function rhs!(dY, Y, p, t)
    dY .= p
    @printf("Inside RHS t %lf dY %lf Y %lf p %lf \n", t, dY[1], Y[1], p[1])
end

# Solve dY = (t+1)dt from Y= 0, t=0
# Using forward Euler
#  dY(0) = 1 -> Y(1) = 1
#  dY(1) = 2 -> Y(2) = 3

p0 = [1.0];
Y0= [0.0];
dt = 1.0;
prob = ODEProblem(rhs!, Y0, (0.0,5.0), p0);

println("Run init")
# Compute dY(t=0), do not update Y, t
integrator = init(prob, Euler(); dt = dt);
@printf("Prior to step t(0) %lf Y(0) %lf  p(0) %lf  \n", integrator.t, integrator.u[1], integrator.p[1])
# update to p[1]
@. integrator.p = [2.0] # p[1] = t[1] + 1
# Compute Y(1), t(1) using dY(0) THEN compute dY(t=1) using these values
step!(integrator, dt, true)
@printf("Post step  t(1) %lf Y(1) %lf  p(1) %lf  \n", integrator.t, integrator.u[1], integrator.p[1])
# update to p[2]
@. integrator.p = [3.0] # p[2] = t[2] + 1
# Compute Y(2), t(2) using dY(1) THEN compute dY(t=1) using these values
step!(integrator, dt, true)
@printf("Post step  t(2) %lf Y(2) %lf  p(2) %lf  \n", integrator.t, integrator.u[1], integrator.p[1])

println("\n")
### What the coupler code does
p0 = [10.0]; # set to nonsense
Y0= [0.0];
dt = 1.0;
prob = ODEProblem(rhs!, Y0, (0.0,5.0), p0);
println("Run init")
integrator2 = init(prob, Euler(); dt = dt); # dY = 0.0 in this case
@. integrator2.p = [1.0]
@printf("Prior to step t(0) %lf Y(0) %lf  p(0) %lf  \n", integrator2.t, integrator2.u[1], integrator2.p[1])
step!(integrator2, dt, true)
@printf("Post step  t(1) %lf Y(1) %lf  p(1) %lf  \n", integrator2.t, integrator2.u[1], integrator2.p[1])
@. integrator2.p = [2.0]
step!(integrator2, dt, true)
@printf("Post step  t(2) %lf Y(2) %lf  p(2) %lf  \n", integrator2.t, integrator2.u[1], integrator2.p[1])
@. integrator2.p = [3.0]
step!(integrator2, dt, true)
@printf("Post step  t(3) %lf Y(3) %lf  p(3) %lf  \n", integrator2.t, integrator2.u[1], integrator2.p[1])


#=
@printf("\n")
prob2 = ODEProblem(rhs!, Y0, (0.0,1.0), [2.0]);
sol = solve(prob2, Euler(); dt = dt)
@printf("Using solve u[1] %lf u[2] %lf u[3] %lf\n", sol.u[1][1], sol.u[2][1],sol.u[3][1])
@printf("Using solve t[1] %lf t[2] %lf t[3] %lf\n", sol.t[1], sol.t[2],sol.t[3])



=#
post_reinit_FA = sum(F_A)
post_reinit_FE = sum(F_E)
Cd_post_renit = sum(p_atmos_0.flux_coefficients.Cd)
Ch_post_reinit= sum(p_atmos_0.flux_coefficients.Ch)
post_reinit_FE = sum(p_atmos_0.dif_flux_energy)
post_reinit_FW = sum(p_atmos_0.dif_flux_ρq_tot)

Cd_pre_step = sum(atmos_sim.integrator.p.flux_coefficients.Cd)
Ch_pre_step = sum(atmos_sim.integrator.p.flux_coefficients.Ch)
prestep_FE = sum(atmos_sim.integrator.p.dif_flux_energy)
prestep_FW = sum(atmos_sim.integrator.p.dif_flux_ρq_tot)


Cd_post_step = sum(atmos_sim.integrator.p.flux_coefficients.Cd)
Ch_post_step = sum(atmos_sim.integrator.p.flux_coefficients.Ch)
poststep_FE = sum(atmos_sim.integrator.p.dif_flux_energy)
poststep_FW = sum(atmos_sim.integrator.p.dif_flux_ρq_tot)

post_push_FE = sum(F_A)
post_push_FW = sum(F_E)
