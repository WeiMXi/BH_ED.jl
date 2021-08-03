using JLD
using BH_ED
using Gadfly
using DataFrames
using Cairo

#]add JLD, Gadfly, DataFrames, Cairo    # 安装依赖

# Tuple 格式数据
# result = Tuple{Int, Int, Float64, Float64, Float64, Float64, Float64, Float64, Float64}[]

# DataFrame格式数据
result_df = DataFrame(N = Int[], Lp = Int[], μ = Float64[], U = Float64[], J = Float64[], V = Float64[], ε = Float64[], E = Float64[], np = Float64[])

for N_ in 6, Lp_ in 3, μ_ in 1.0, U_ in 1.0, J_ in 1.0, V_ in 0.1:0.1:2.0, ε_ in 1.0
# for N_ in 2:4, Lp_ in 1:4, μ_ in 0.1:0.2:1.0, U_ in 0.1:0.2:1.0, J_ in 0.1:0.2:1.0, V_ in 0.1:0.2:1.0, ε_ in 0.1:0.2:1.0
    sol = BH_solver(N = N_, Lp = Lp_, μ = μ_, U = U_, J = J_, V = V_, ε = ε_, screen_p = false)
    # push!( result, (N_, Lp_, μ_, U_, J_, V_, ε_, sol[1][1], sol[3])) )
    push!(result_df, ( N_, Lp_, μ_, U_, J_, V_, ε_, sol[1][1], sol[3]) )
end

# 保存数据
# save("data.jld", "result", result_df)

# 绘图
mypl1 = plot(result_df, x = :V, y = :E, Geom.point, Geom.line)
# 保存图片
# img = SVG("N = 6, Lp = 3.svg", 12inch, 8inch) # SVG 格式
# img = PDF("N = 6, Lp = 3.pdf", 12inch, 8inch) # PDF 格式
# img = SVG("N = 6, Lp = 3.png", 12inch, 8inch) # PNG 格式
# draw(img, mypl1)  # 输出文件