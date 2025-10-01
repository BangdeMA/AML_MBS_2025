# Author: Dr. Shizhao MA (shizhaoma@sjtu.edu.cn)
# Date: 2025-04-15
# Description: This is the code for the optimal control, coupled with Latin hyper sampling and death mechanism.
# --------------------------------------------------------
# 1. Load packages
# If you don't have the packages installed, you can install them by running the following commands:
# Pkg.add("Pkg")
using Pkg
Pkg.activate(".")

using Random
using Distributions
using DelimitedFiles
using ComponentArrays
using OptimizationOptimisers
using Plots
using StatsBase
using DiffEqFlux
using Optimization
#using OptimizationOptimJL
using Random
using DataFrames
using CSV
using LaTeXStrings

rng = Random.default_rng()
Random.seed!(1);

# ----------------------------------------------------------------------------------------------------
## Run the dynamic of Single patient
function RunTherapy(Num, beta_L0, beta_L2, d_L2, para)
    T = 3 * 365
    dt = 0.001
     N = Int(T / dt)
    H1 = zeros(N + 1)
    H2 = zeros(N + 1)
    H3 = zeros(N + 1)
    B1 = zeros(N + 1)
    B2 = zeros(N + 1)
    B3 = zeros(N + 1)
    L1 = zeros(N + 1)
    L2 = zeros(N + 1)
     R = zeros(N + 1)
     P = zeros(N + 1)
     sigma = zeros(N + 1)
    lambda = zeros(N + 1)
    D1 = zeros(N + 1)
    D2 = zeros(N + 1)
    
    epsilon_1 = 1.0
    epsilon_2 = 1.0
    epsilon_3 = 1.0
    epsilon_4 = 1.0
    
    K_H1B1 = 5e4
    K_H2B3 = 2.6444e+08
    K_B1L2 = 7.5e6
    K_L2B2 = 2.9042e10 #3e10
    
    theta_H1H1 = 6.1566e-07
    theta_H1L1 = 5e-3
    theta_H1H3 = 3.1498e-12
    theta_B1B1 = 3.8095e-05
    theta_B2B2 = 2.8154e-09
    theta_B2L2 = 1e-5
    theta_L1H1 = 1e-10
    theta_L1L1 = 2e-7

    beta_H0 = 0.5
    beta_H2 = 0.2450
    beta_B1 = 0.4800
    beta_B2 = 0.1852

    kappa_H0 = 0.0455
    kappa_H2 = 0.5
    kappa_B1 = 0.1000
    kappa_B2 = 0.1000
    kappa_L1 = 0.2500

    d_H1 = 0.0005 * 0
    d_H2 = 0.0213 * 0
    d_H3 = 0.0080
    d_B1 = 0.0600
    d_B2 = 0.0625
    d_B3 = 0.1000
    d_L1 = 0.0688

    delta_L = 4.5e-3 * 0.005
    epsilon_0 = 7.0
    m_B2 = 5.0
    theta_B2 = 1e8

    #mu_p = 0.50
    #sigma_p = 0.065
    death = 0.0

    B1star = (beta_B1 / (kappa_B1 + d_B1) - 1) / theta_B1B1
    first = B1star * kappa_B1 * theta_B2B2 + beta_B2 - d_B2 - kappa_B2
    B2star = (sqrt(first^2 + 4 * B1star * kappa_B1 * theta_B2B2 * (d_B2 + kappa_B2)) + first) / 2 / theta_B2B2 / (d_B2 + kappa_B2)
    B3star = kappa_B2 * B2star / d_B3
    B_H1 = beta_H0 * (1 + B1star / (B1star + K_H1B1))
    K_H2 = kappa_H2 * (1 + B3star / (B3star + K_H2B3))
    FstH1 = -K_H2 * (B_H1^2 * theta_H1H3 + d_H3 * theta_H1H1 * (B_H1 - 2 * kappa_H0)) + beta_H2 * d_H3 * theta_H1H1 * (B_H1 - 2 * kappa_H0)
    SecH1 = (B_H1 * K_H2 * theta_H1H3 + d_H3 * theta_H1H1 * (beta_H2 - K_H2))^2 - 4 * d_H3 * K_H2 * theta_H1H1 * theta_H1H3 * (B_H1 - kappa_H0) * (beta_H2 - K_H2)
    H1star = (FstH1 - B_H1 * sqrt(SecH1)) / (2 * d_H3 * kappa_H0 * theta_H1H1 * theta_H1H1 * (beta_H2 - K_H2))
    FstH2 = B_H1 * K_H2 * theta_H1H3 + beta_H2 * d_H3 * theta_H1H1 - d_H3 * K_H2 * theta_H1H1
    H2star = (-FstH2 - sqrt(SecH1)) / (2 * K_H2 * theta_H1H1 * theta_H1H3 * (beta_H2 - K_H2))
    H3star = K_H2 * H2star / d_H3
    H1[1] = H1star
    H2[1] = H2star
    H3[1] = H3star
    B1[1] = B1star
    B2[1] = B2star
    B3[1] = B3star
    L1[1] = 1e2
    L2[1] = 1e2
     R[1] = L2[1] / (H3[1] + L2[1])
     P[1] = 0.0
    D1[1] = 0.0
    D2[1] = 0.0

    t = zeros(N)
    k = 1
    OS = 0.0; 
    UST = 0.05+0.15*rand()
    #println("UST=$UST")
    for i in 1:N
        lambda[i] = R[i] * sigma[i] * H1[i]
        sigma[i] = delta_L * (epsilon_0 + (B2[i] / theta_B2)^m_B2) / (1 + (B2[i] / theta_B2)^m_B2)

        if abs(R[i] - UST) < 1e-3
            t[k] = i
            k += 1
            #println("$Num $(R[i]) $i $k $(t[1] * dt)")
        end

        if i > t[1] && t[1] != 0
            D1[i] = para[3]#0.30
            D2[i] = D1[i]
        else
            D1[i] = 0.0
            D2[i] = 0.0
        end

        mutation = rand(Poisson(lambda[i] * dt))
        H1[i+1] = H1[i] + dt * H1[i] * (beta_H0 / (1 + theta_H1H1 * H1[i] + theta_H1L1 * L1[i]) * (1 + epsilon_1 * B1[i] / (B1[i] + K_H1B1))) - dt * H1[i] * kappa_H0 / (1 + theta_H1H3 * H3[i]) - dt * d_H1 * H1[i] - mutation
        H2[i+1] = H2[i] + dt * H1[i] * kappa_H0 / (1 + theta_H1H3 * H3[i]) + dt * H2[i] * beta_H2 - dt * H2[i] * kappa_H2 * (1 + epsilon_2 * B3[i] / (B3[i] + K_H2B3)) - dt * d_H2 * H2[i]
        H3[i+1] = H3[i] + dt * H2[i] * kappa_H2 * (1 + epsilon_2 * B3[i] / (B3[i] + K_H2B3)) - dt * H3[i] * d_H3
        B1[i+1] = B1[i] + dt * B1[i] * beta_B1 / (1 + theta_B1B1 * B1[i]) - dt * B1[i] * kappa_B1 * (1 + epsilon_3 * L2[i] / (K_B1L2 + L2[i])) - dt * d_B1 * B1[i]
        B2[i+1] = B2[i] + dt * B1[i] * kappa_B1 * (1 + epsilon_3 * L2[i] / (K_B1L2 + L2[i])) + dt * B2[i] * beta_B2 / (1 + theta_B2B2 * B2[i]) - dt * B2[i] * kappa_B2 / (1 + theta_B2L2 * L2[i]) - dt * d_B2 * B2[i]
        B3[i+1] = B3[i] + dt * B2[i] * kappa_B2 / (1 + theta_B2L2 * L2[i]) - dt * d_B3 * B3[i]
        L1[i+1] = L1[i] + dt * L1[i] * (beta_L0 / (1 + theta_L1H1 * H1[i] + theta_L1L1 * L1[i]) - kappa_L1 - d_L1 - D1[i]) + mutation
        L2[i+1] = L2[i] + dt * L1[i] * kappa_L1 + dt * L2[i] * beta_L2 * (1 + epsilon_4 * B2[i] / (K_L2B2 + B2[i])) - dt * (d_L2 + D2[i]) * L2[i]
         R[i+1] = L2[i+1] / (H3[i+1] + L2[i+1])
         P[i+1] = rand()
        #death = (1.0 + exp(-(1 - 0.5) / 0.082)) / (1.0 + exp(-(R[i+1] - 0.5) / 0.082));
        death = (1.0 + exp(-(1 - para[1]) / para[2])) / (1.0 + exp(-(R[i+1] - para[1]) / para[2]));
        if P[i+1] < death[1] * dt
            #println("$Num: death probability=$(P[i+1]) $(death * dt) R=$(R[i+1])")
            OS = (i+1)*dt
            break
        end
        if i == N
            #println("$Num: death probability=$(P[i+1]) $(death * dt) R=$(R[i+1])")
            OS = N*dt
            break
        end
    end
    return OS, t[1]
end

# ----------------------------------------------------------------------------------------------------
## generate the parameters file of Virtual Patients (Digtal Twins)
# Pkg.add("LatinHypercubeSampling")
# Pkg.add("QuasiMonteCarlo")
# Pkg.add("Surrogates")
#------------------------------------------------------------------------------------------------------
using QuasiMonteCarlo
using DelimitedFiles

function VPparameters(Num,dsz)
    # baseline
    beta_L0 = 0.5000
    beta_L2 = 0.1240
    d_L2    = 0.1330
    # dsz = 2.0


    lower_bounds = [beta_L0 * dsz[4] / dsz[1], beta_L2 * dsz[5] / dsz[2], d_L2 * dsz[6] / dsz[3]]
    upper_bounds = [beta_L0 * dsz[4] * dsz[1], beta_L2 * dsz[5] * dsz[2], d_L2 * dsz[6] * dsz[3]]

    # LSH
    Y = QuasiMonteCarlo.sample(Num, lower_bounds, upper_bounds, LatinHypercubeSample())

    # save as ASCII 
    writedlm("vpatients.txt", Y')

    return Y'
end

x = 0:1095;

# -----------------------------------------------------------------------------------------
# ## Read the Overall Survival clinical data 
file2 = open("TCGA_OS.csv", "r") 
OSclinical = readdlm(file2)[:,1] 
close(file2) 
# ----------------------------------------------------------------------------------------------------
gcdfclinmcal   = StatsBase.ecdf(OSclinical)
Clinicaldata   = gcdfclinmcal(x);
# loss = sum(abs2,(Clinicaldata .- Simulationdata))

# plot
plot(x, Clinicaldata, label="ECDF of OS", xlabel="Overall Survival (OS)", ylabel="Cumulative Probability", lw=2)
## -----------------------------------------------------------------------------------------------------
using Optimization, OptimizationOptimisers

# Define the loss function
function loss_function(para)
    Num = 1000
    D = VPparameters(Num, para[4:9])
    
    # save OS[i] with T[i] > 0  
    OS_selected = Float64[]

    for i in 1:Num
        os, t = RunTherapy(i, D[i,1], D[i,2], D[i,3], para[1:3])
        if t > 0
            push!(OS_selected, os)
        end
    end

    #  ecdf
    ecdf_func = StatsBase.ecdf(OS_selected)
    sim_data = ecdf_func.(x)
    return sum(abs2, Clinicaldata .- sim_data)
end
# function loss_function(para)
#     D = VPparameters(Num,para[4]);
#     sim_data = StatsBase.ecdf([RunTherapy(i, D[i,1], D[i,2], D[i,3], para[1:3]) for i in 1:Num])(x)
#     return sum(abs2, Clinicaldata .- sim_data)
# end

# intial value
prob = OptimizationProblem(OptimizationFunction((p, _) -> loss_function(p), Optimization.AutoForwardDiff()), [0.5692850210234907,0.0789246091076159,0.3610983400073908,2.550633099502173, 2.5285273106986703, 2.3561948714475895, 0.40096147583088987, 2.213659193380096, 1.450483904772558])


function callback_function(state, loss)
    println("Iteration: Current Loss = $(loss), Parameters = $(state.u)")
    return false
end

using OptimizationOptimJL  # use BFGS need OptimizationOptimJL
sol = solve(prob, Optim.NelderMead(), maxiters=1000, callback=callback_function)
println("the optimal parameter set: ", sol.u)
