# BH_ED

[![Build Status](https://github.com/WeiMXi/BH_ED.jl/workflows/CI/badge.svg)](https://github.com/WeiMXi/BH_ED.jl/actions)
[![Coverage](https://codecov.io/gh/WeiMXi/BH_ED.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/WeiMXi/BH_ED.jl)

It's a Julia program for repetition of some work in [Exact diagonalization: Bose–Hubbard model | stone-zeng.github.io](https://stone-zeng.github.io/2019-10-03-exact-diagonalization/)

But our work is a bit different with the reference.

Done by Dr. Ning and me.

Our Hamiltonion formula:

![](./equation.svg)

## For start, just

```Julia
]add https://github.com/WeiMXi/BH_ED.jl
using BH_ED
#]up BH_ED # for update the pkg
```

### Get the solver

```julia
- N: 4 with Lp: 2
- We got all required states.
- Dimension of Hilbert space: 81
- We get Hamiltonion Matrix.
- Number of nonzero elements: 584 (8.901%) of total.
- Hamiltonion Matrix is a Real Symmetric Matrix.
- Ground energy is -7.068791063296066
- The corresponding number of particles is 2.8845941217640876
- Then We will return the eigenvalue, the corresponding eigenvector and the number of particles
([-7.068791063296066], [[1.317044712859158e-14, 0.1395381558487617, 0.1395381558487586, 0.13953815584876542, 0.13953815584876827, 2.903609550038398e-14, 8.067861312949678e-15, 4.287966882712928e-15, -1.6260138595209274e-15, 8.188930552422516e-15  …  -6.203944800148619e-16, -4.788552859264462e-16, -1.6613819558997184e-15, 3.070263591561661e-15, 4.486288336861432e-15, 0.0028170873328648574, 0.002817087332865875, 0.0028170873328649403, 0.0028170873328656407, -2.7849159609571244e-17]], 2.8845941217640876)

julia> sol[1]   # the eigenvaluse
1-element Vector{Float64}:
 -7.068791063296066

julia> sol[2]   #  the corresponding eigenvector
1-element Vector{Vector{Float64}}:
 [1.317044712859158e-14, 0.1395381558487617, 0.1395381558487586, 0.13953815584876542, 0.13953815584876827, 2.903609550038398e-14, 8.067861312949678e-15, 4.287966882712928e-15, -1.6260138595209274e-15, 
8.188930552422516e-15  …  -6.203944800148619e-16, -4.788552859264462e-16, -1.6613819558997184e-15, 3.070263591561661e-15, 4.486288336861432e-15, 0.0028170873328648574, 0.002817087332865875, 0.0028170873328649403, 0.0028170873328656407, -2.7849159609571244e-17]

julia> sol[3]   # the number of particles of eigenvector
2.8845941217640876
```

### Get the Hamiltonion Matrix

```julia
julia> them = give_a_BH_Model(N = 6, Lp = 3);

julia> theh = give_H(them)
4096×4096 SparseArrays.SparseMatrixCSC{Float64, Int64} with 56313 stored entries:
⣿⣿⣛⠳⢦⣴⣄⡀⢀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⢿⡘⠿⣧⣀⡀⠉⠛⠾⣷⣄⡀⢀⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⢈⣷⠀⠸⢿⣷⡄⠀⠀⠀⠙⠳⣦⣙⢦⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠹⣧⠀⠀⠉⣻⣾⣦⡀⠀⠀⠀⠉⠹⣷⣀⠳⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠐⢾⣧⠀⠀⠈⠻⢿⣷⡄⠀⠀⠀⠀⠈⠛⢶⣌⡃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠹⢷⡀⠀⠀⠀⠉⣿⣿⣦⡀⠀⠀⠀⠀⠈⠻⣦⡙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠰⣌⢻⡄⠀⠀⠀⠈⠻⣿⣿⣆⠀⠀⠀⠀⠀⠙⢷⣦⡙⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠈⠓⢷⣦⡀⠀⠀⠀⠈⠙⢻⣶⣦⡀⠀⠀⠀⠀⠈⠻⣆⠲⣄⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⢤⡘⢻⣄⠀⠀⠀⠀⠈⠻⣟⣽⣦⡀⠀⠀⠀⠀⠙⣧⡌⠓⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠙⠦⠹⣦⡀⠀⠀⠀⠀⠈⠻⠿⣧⣄⡀⠀⠀⠀⠈⠻⢷⢤⡀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣌⠻⢷⣄⠀⠀⠀⠀⠀⠹⣿⣿⣦⡀⠀⠀⠀⠘⣧⡙⠆⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠳⣌⠻⣦⡀⠀⠀⠀⠀⠈⠻⣿⣿⣀⠀⠀⠀⠈⢷⣆⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢨⡙⠷⣤⡀⠀⠀⠀⠀⠘⢿⣷⣦⡀⠀⠀⢻⡷⠄⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⠉⢿⣆⣀⠀⠀⠀⠈⠻⡿⣯⣀⠀⠀⢻⣆⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠳⣍⠻⢦⣄⠀⠀⠀⠘⢿⣷⡆⠀⢿⡁
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠁⠈⠙⢿⡶⣤⣀⠈⠉⢻⣶⡌⣷
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠁⠈⠙⠟⠳⢦⣭⣿⣿

```

## Help
```
?BH_soler
?give_a_BH_Model
?give_give_H
```

![200x200](N4_Lp2_all1_G.svg)

