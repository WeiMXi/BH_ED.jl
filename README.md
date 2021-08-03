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
 julia> sol = BH_solver()
- N: 4 with Lp: 2
- We got all required states.
- Dimension of Hilbert space: 81
- We get Hamiltonion Matrix.
- Number of nonzero elements: 580 (8.84%) of total.
- Hamiltonion Matrix is a Real Symmetric Matrix.
- Ground energy is -5.060635855555107
- The corresponding number of particles is 2.1658167162979147
- Then We will return the eigenvalue, the corresponding eigenvector and the number of particles
([-5.060635855555107], [[0.15442257743686513, -5.639815515016311e-15, -9.608770168492102e-15, 5.0204336273891484e-15, -5.742905228043529e-16, 0.27629364230280484, 0.30101692306920486, 0.27629364230281034, 0.31790851466867864, 0.3010169230692081  …  0.006707879845250392, 0.0063690847741138725, 0.008191147222244531, 0.006707879845246998, 0.008191147222239592, -3.567169362373543e-15, 7.916709677476003e-15, -3.62708759182321e-15, 7.369954646503184e-15, 0.0005642402389404083]], 2.1658167162979147)

julia> sol[1]   # the eigenvaluse
1-element Vector{Float64}:
 -5.060635855555107

julia> sol[2]   #  the corresponding eigenvector
1-element Vector{Vector{Float64}}:
 [0.15442257743686513, -5.639815515016311e-15, -9.608770168492102e-15, 5.0204336273891484e-15, -5.742905228043529e-16, 0.27629364230280484, 0.30101692306920486, 0.27629364230281034, 0.31790851466867864, 0.3010169230692081  …  0.006707879845250392, 0.0063690847741138725, 0.008191147222244531, 0.006707879845246998, 0.008191147222239592, -3.567169362373543e-15, 7.916709677476003e-15, -3.62708759182321e-15, 7.369954646503184e-15, 0.0005642402389404083]

julia> sol[3]   # the number of particles of eigenvector
2.1658167162979147
```

### Get the Hamiltonion Matrix

```julia
julia> theh = give_H(them)
4096×4096 SparseArrays.SparseMatrixCSC{Float64, Int64} with 56293 stored entries:
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

