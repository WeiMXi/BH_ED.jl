# baiaiyu
# 2021,Jul,31
# Repetition of https://stone-zeng.github.io/2019-10-03-exact-diagonalization/
# 统摄不懈，只为破敌，岂可擅退

module BH_ED

using SparseArrays
using LinearAlgebra
using KrylovKit
using Memoize
# using Distributed

struct BH_Model
    states::Vector{String}
    # states_n::Matrix{Int}
    N::Int
    L::Int
    Lp::Int
    μ::Float64
    U::Float64
    J::Float64
    V::Float64
    ε::Float64
end

include("./bh_hamiltonian.jl")

export give_a_BH_Model, give_H, BH_solver, give_the_num_of_particles_of_eigvector


@doc raw"""
    give_a_BH_Model(; N::Int        = 4
                    , Lp::Int       = 2
                    , L::Int        = N*Lp
                    , μ::Float64    = 1.0
                    , U::Float64    = 1.0
                    , J::Float64    = 1.0
                    , V::Float64    = 1.0
                    , ε::Float64    = 1.0
                    )

定义一个Bose–Hubbard model

```math
\hat{H}=-J \sum_{i=1}^{N}\left(\hat{a}_{i}^{\dagger} \hat{a}_{i+1}+\hat{a}_{i+1}^{\dagger} \hat{a}_{i}\right)+\sum_{i=1}^{N}-\mu \hat{n}_{i}+\frac{U}{2} \hat{n}_{i}\left(\hat{n}_{i}-1\right)+V \sum_{i=1}^{N} \hat{n}_{i} \hat{n}_{i+1}-\frac{\varepsilon}{2} \sum_{i=1}^{N}\left(\hat{a}_{i}^{2}+\hat{a}_{i+1}^{+2}\right)
```
# Arguments
- `N::Int`: 一维链长度
- `Lp::Int`: 每个格点最大容纳粒子数
- `L::Int`: 最多粒子数考虑 所有 0 到 L 个粒子的模型；最大可设置为 Lp*N，

"""
function give_a_BH_Model(;N::Int        = 4
                        , Lp::Int       = 2
                        , L::Int        = N*Lp
                        , μ::Float64    = 1.0
                        , U::Float64    = 1.0
                        , J::Float64    = 1.0
                        , V::Float64    = 1.0
                        , ε::Float64    = 1.0
                        )
    states = String[]
    for i in 0:L
        append!( states, give_all_states(N, i, Lp) )
    end
    # states = Array{Int}( undef, 2, length(states) )
    # states_n[1, :] .= [i for i in 1:length(states)]
    # states_n[2, :] .= map(states) do x
    #     parse(Int, x, base = Lp + 1)
    # end
    BH_Model(states, N, L, Lp, μ, U, J, V, ε)
    # BH_Model(states, states_n, N, L, Lp, μ, J, V, ε)

end


```
给出一个BH模型的哈密顿矩阵
```
function give_H(the_model::BH_Model)
    nhil = length(the_model.states)
    H = spzeros(Float64, nhil, nhil)

    # if nhil > 10000
    #     H = spzeros(Float64, nhil, nhil)
    # else
    #     H = zeros(Float64, nhil, nhil)
    # end

    if the_model.μ != 0.0
        fill_H_μ_U!(H, the_model)
    end

    if the_model.J != 0.0
        fill_H_J!(H, the_model)
    end

    if the_model.V != 0.0
        fill_H_V!(H, the_model)
    end

    if the_model.ε != 0.0
        fill_H_ε!(H, the_model)
    end

    H
end

# ```
# 给出一个BH模型的哈密顿矩阵, 并行计算，最大可利用 4 核心
# ```
# function give_H_p(the_model::BH_Model)
#     nhil = length(the_model.states)
#     H = spzeros(Float64, nhil, nhil)
#     h1 = spzeros(Float64, nhil, nhil)
#     h2 = spzeros(Float64, nhil, nhil)
#     h3 = spzeros(Float64, nhil, nhil)
#     h4 = spzeros(Float64, nhil, nhil)
#     call_result =[]
#     # if nhil > 10000
#     #     H = spzeros(Float64, nhil, nhil)
#     # else
#     #     H = zeros(Float64, nhil, nhil)
#     # end

#     if the_model.μ != 0.0
#         # @spawn push!( call_state, fill_H_μ_U!(h1, the_model) )
#         push!( call_result, remotecall(BH_ED.fill_H_μ_U!, workers()[1], h1, the_model) )
#     end

#     if the_model.J != 0.0
#         # @spawn call_state[2] = fill_H_J!(h2, the_model)
#         push!( call_result, remotecall(BH_ED.fill_H_J!, workers()[2], h2, the_model) )
#     end

#     if the_model.V != 0.0
#         # @spawn call_state[3] = fill_H_V!(h3, the_model)
#         push!( call_result, remotecall(BH_ED.fill_H_V!, workers()[3], h3, the_model) )
#     end

#     if the_model.ε != 0.0
#         # @spawn call_state[4] = fill_H_ε!(h4, the_model)
#         push!( call_result, remotecall(BH_ED.fill_H_ε!, workers()[4], h4, the_model) )
#     end

#     for r in result
#         H = H + fetch(r)
#     end

#     H
# end

@doc raw"""
    BH_solver(; N::Int        = 4
              , Lp::Int       = 2
              , L::Int        = N*Lp
              , μ::Float64    = 1.0
              , U::Float64    = 1.0
              , J::Float64    = 1.0
              , V::Float64    = 1.0
              , ε::Float64    = 1.0
              , screen_p::Bool= true 
                    )

定义并求解 Bose–Hubbard model, 返回最小特征值与对应的特征向量

```math
\hat{H}=-J \sum_{i=1}^{N}\left(\hat{a}_{i}^{\dagger} \hat{a}_{i+1}+\hat{a}_{i+1}^{\dagger} \hat{a}_{i}\right)+\sum_{i=1}^{N}-\mu \hat{n}_{i}+\frac{U}{2} \hat{n}_{i}\left(\hat{n}_{i}-1\right)+V \sum_{i=1}^{N} \hat{n}_{i} \hat{n}_{i+1}-\frac{\varepsilon}{2} \sum_{i=1}^{N}\left(\hat{a}_{i}^{2}+\hat{a}_{i+1}^{+2}\right)
```
# Arguments
- `N::Int`: 一维链长度
- `Lp::Int`: 每个格点最大容纳粒子数
- `L::Int`: 最多粒子数考虑 所有 0 到 L 个粒子的模型；最大可设置为 Lp*N
- `screen_p::Bool`: 是否屏幕打印信息

"""
function BH_solver(;N::Int      = 4
                , Lp::Int       = 2
                , L::Int        = N*Lp
                , μ::Float64    = 1.0
                , U::Float64    = 1.0
                , J::Float64    = 1.0
                , V::Float64    = 1.0
                , ε::Float64    = 1.0
                , screen_p::Bool= true 
    )::Tuple{Vector{Float64}, Vector{Vector{Float64}}, Float64}

    the_bhm = give_a_BH_Model(N = N, Lp = Lp, L = L, μ = μ, U = U, J = J, V = V, ε = ε)
    screen_p && println("- N: $(the_bhm.N) with Lp: $(the_bhm.Lp)")
    screen_p && println("- We got all required states.")
    screen_p && println("- Dimension of Hilbert space: $(length(the_bhm.states))")
    H = give_H(the_bhm)
    screen_p && println("- We get Hamiltonion Matrix.")
    screen_p && print("- Number of nonzero elements: ")
    screen_p && print(nnz(H))
    screen_p && print(" (",round((100.0*nnz(H))/(H.m*H.n), digits = 3))
    screen_p && println("%) of total.")
    if H == H'
        screen_p && println("- Hamiltonion Matrix is a Real Symmetric Matrix.")
    else
        throw("Hamiltonion Matrix is NOT a Real Symmetric Matrix!")
        # return H
    end
    eigsol = eigsolve(H, 1, :SR, eager = true)
    n_ev = give_the_num_of_particles_of_eigvector(eigsol[2][1], the_bhm)
    screen_p && println("- Ground energy is $(eigsol[1][1])")
    screen_p && println("- The corresponding number of particles is $n_ev")
    screen_p && println("- Then We will return the eigenvalue, the corresponding eigenvector and the number of particles")
    eigsol[1], eigsol[2], n_ev
    # H
end

end
