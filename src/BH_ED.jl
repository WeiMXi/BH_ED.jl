# baiaiyu
# 2021,Jul,31
# Repetition of https://stone-zeng.github.io/2019-10-03-exact-diagonalization/
# 统摄不懈，只为破敌，岂可擅退

module BH_ED

using SparseArrays
using LinearAlgebra
using KrylovKit
using Memoize

struct BH_Model
    states::Vector{String}
    # states_n::Matrix{Int}
    N::Int
    L::Int
    Lp::Int
    μ::Float64
    J::Float64
    V::Float64
    ε::Float64
end

include("./bh_hamiltonian.jl")

export give_a_BH_Model, give_H, BH_solver


@doc raw"""
    give_a_BH_Model(; N::Int        = 4
                    , Lp::Int       = 2
                    , L::Int        = N*Lp
                    , μ::Float64    = 1.0
                    , J::Float64    = 1.0
                    , V::Float64    = 1.0
                    , ε::Float64    = 1.0
                    )

定义一个Bose–Hubbard model

```math
\hat{H}=-J \sum_{\langle i j\rangle}\left(\hat{a}_{i}^{\dagger} \hat{a}_{j}+\hat{a}_{j}^{\dagger} \hat{a}_{i}\right)+\frac{\mu}{2} \sum_{i=1}^{N} ( \hat{n}_{i}\left( \hat{n}_{i}-1\right) - 2 \hat{n}_{i} )  + V \sum_{i=1}^{N} \hat{n}_{i}\hat{n}_{i + 1} - \frac{\varepsilon}{2} \sum_{i=1}^{N}( \hat{a}_{i}^{2} + \hat{a}_{j}^{\dagger 2})
```
# Arguments
- `N::Int`: 一维链长度
- `Lp::Int`: 每个格点最大容纳粒子数
- `L::Int`: 最多粒子数考虑 所以 0 到 L 个粒子的模型；最大可设置为 Lp*N，

"""
function give_a_BH_Model(;N::Int        = 4
                        , Lp::Int       = 2
                        , L::Int        = N*Lp
                        , μ::Float64    = 1.0
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
    BH_Model(states, N, L, Lp, μ, J, V, ε)
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
    fill_H_μ!(H, the_model)
    fill_H_J!(H, the_model)
    fill_H_V!(H, the_model)
    fill_H_ε!(H, the_model)
    H
end

@doc raw"""
    BH_solver(; N::Int        = 4
              , Lp::Int       = 2
              , L::Int        = N*Lp
              , μ::Float64    = 1.0
              , J::Float64    = 1.0
              , V::Float64    = 1.0
              , ε::Float64    = 1.0
                    )

定义并求解 Bose–Hubbard model, 返回最小特征值与对应的特征向量

```math
\hat{H}=-J \sum_{\langle i j\rangle}\left(\hat{a}_{i}^{\dagger} \hat{a}_{j}+\hat{a}_{j}^{\dagger} \hat{a}_{i}\right)+\frac{\mu}{2} \sum_{i=1}^{N} ( \hat{n}_{i}\left( \hat{n}_{i}-1\right) - 2 \hat{n}_{i} )  + V \sum_{i=1}^{N} \hat{n}_{i}\hat{n}_{i + 1} - \frac{\varepsilon}{2} \sum_{i=1}^{N}( \hat{a}_{i}^{2} + \hat{a}_{j}^{\dagger 2})
```
# Arguments
- `N::Int`: 一维链长度
- `Lp::Int`: 每个格点最大容纳粒子数
- `L::Int`: 最多粒子数考虑 所以 0 到 L 个粒子的模型；最大可设置为 Lp*N，

"""
function BH_solver(;N::Int       = 4
                , Lp::Int       = 2
                , L::Int        = N*Lp
                , μ::Float64    = 1.0
                , J::Float64    = 1.0
                , V::Float64    = 1.0
                , ε::Float64    = 1.0
    )

    the_bhm = give_a_BH_Model(N = N, Lp = Lp, L = L, μ = μ, J = J, V = V, ε = ε)
    println("N: $(the_bhm.N) with Lp: $(the_bhm.Lp)")
    println("We got all required states.")
    println("Dimension of Hilbert space: $(length(the_bhm.states))")
    H = give_H(the_bhm)
    println("We get Hamiltonion Matrix.")
    print("Number of nonzero elements: ")
    print(nnz(H))
    print(" (",round((100.0*nnz(H))/(H.m*H.n), digits = 3))
    println("%) of total.")
    if H == H'
        println("Hamiltonion Matrix is a Real Symmetric Matrix.")
    else
        throw("Hamiltonion Matrix is NOT a Real Symmetric Matrix!")
        # return H
    end
    eigsol = eigsolve(H, 1, :SR, eager = true)
    println("Ground energy is $(eigsol[1][1])")
    println("Then We will return the eigenvalue and the corresponding eigenvector")
    eigsol
    # H
end

end
