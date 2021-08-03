function give_all_states(N::Int, L::Int, Lp::Int)

    # states = Array{String}( undef, (Lp + 1)^N )
    states = String[]

    # 以 Lp + 1 进制表示，生成最大值
    max_s = string(Lp, base = Lp + 1)
    for i in 2:N
        max_s = max_s * string(Lp, base = Lp + 1)
    end

    # 从0到max_s，找出所有符合粒子总数条件的状态
    for i in 0:parse(Int, max_s, base = Lp + 1)
        # 进制转换为 Lp + 1 
        the_state = string(i, base = Lp + 1, pad = N)
        # 粒子数符合条件则记录
        if total_num_particle(the_state, N, Lp) == L
            push!(states, the_state)
        end
    end
    # parse(Int, "444", base =5)
    # string(4444, base = 5)
    states
    # max_s
end

@memoize function total_num_particle(the_state::String, N::Int, Lp::Int)
    num = 0
    for i in 1:N
        # num = num + parse(Int, the_state[i], base = Lp+1)
        num = num + num_of_position(the_state[i], Lp)
    end
    num
end

@memoize function num_of_position(the_state_of_i::Char, Lp::Int)
    parse(Int, the_state_of_i, base = Lp + 1)
end


function fill_H_μ_U!(H, the_model::BH_Model)
    for i in 1:length(the_model.states)
        the_num = total_num_particle(the_model.states[i], the_model.N, the_model.Lp)
        H[i, i] = H[i, i] - (0.5*the_model.U + the_model.μ)*the_num + 0.5*the_model.U*(the_num^2)
        # H[i, i] = H[i, i] + the_model.μ*( -1.5*the_num + 0.5*the_num^2)
    end
    true
end

function fill_H_J!(H, the_model::BH_Model)
    for j in 1:length(the_model.states)
        for i in 1:the_model.N
            
            index, quan = give_index_cr_an(the_model.states, j, i, i + 1, the_model.N, the_model.Lp)
            if index != 0
                H[index, j] = H[index, j] - the_model.J*quan
                # println("$i-cr, $(i + 1)-an $(the_model.states[j]) -> $(the_model.states[index]) | ($j, $index), $quan")  # debug
            end
            index, quan = give_index_cr_an(the_model.states, j, i + 1, i, the_model.N, the_model.Lp)
            if index != 0
                H[index, j] = H[index, j] - the_model.J*quan
                # println("$(i + 1)-cr, $i-an $(the_model.states[j]) -> $(the_model.states[index]) |($j, $index), $quan")  # debug
            end
        end

    end
    true

end


# i产生, j湮灭
@memoize function give_index_cr_an(states::Vector{String}, label::Int, i::Int, j::Int, N::Int, Lp::Int)
    the_state_s = states[label]
    i, j = adjust_index_for_pbc(i, j, N)
    i_num = num_of_position(the_state_s[end - i + 1], Lp)
    j_num = num_of_position(the_state_s[end - j + 1], Lp)
    if j_num == 0 || i_num == Lp
        return 0, 0.0
    end
    the_state_n = parse(Int, the_state_s, base = Lp + 1)
    the_state_n = the_state_n + (Lp + 1)^(i - 1) - (Lp + 1)^(j - 1)
    the_state_s = string(the_state_n, base = Lp + 1, pad = N)
    # findall(x -> x == the_state_s, states), i_num, j_num, states[label], the_state_s # debug

    findall(x -> x == the_state_s, states)[1], sqrt(i_num + 1)*sqrt(j_num) # release

end

function fill_H_V!(H, the_model::BH_Model)

    for j in 1:length(the_model.states)
        the_sum = 0.0
        for i in 1:the_model.N
            x, y = adjust_index_for_pbc(i, i + 1, the_model.N)
            the_sum = the_sum + the_model.V*num_of_position(the_model.states[j][x], the_model.Lp)*num_of_position(the_model.states[j][y], the_model.Lp)
        end
        H[j, j] = H[j, j] + the_sum
    end

    # for i in 1:length(the_model.states)
    #     i, j = adjust_index_for_pbc(i, i + 1, the_model.N)
    #     H[i, i] = H[i, i] + the_model.V*total_num_particle(the_model.states[i], the_model.N, the_model.Lp)*total_num_particle(the_model.states[j], the_model.N, the_model.Lp)
    # end
    true

end

function fill_H_ε!(H, the_model::BH_Model)
    for j in 1:length(the_model.states)

        for i in 1:the_model.N
            index_1 , v1, index_2, v2 = give_index_cr_and_an_2times(the_model.states, j, i, the_model.Lp, the_model.N)
            if index_1 != 0
                H[index_1, j] = H[index_1, j] - 0.5*the_model.ε*v1
            end
            if index_2 != 0
                H[index_2, j] = H[index_2, j] - 0.5*the_model.ε*v2
            end
        end

    end
    true

end

@memoize function give_index_cr_and_an_2times(states::Vector{String}, label::Int, i::Int, Lp::Int, N::Int)
    the_state_s = states[label]
    i_num = num_of_position(the_state_s[end - i + 1], Lp)
    index_1 , v1, index_2, v2 = 0, 0.0, 0, 0.0
    # 产生两次
    if i_num <= Lp - 2
        the_state_1_n = parse(Int, the_state_s, base = Lp + 1)
        the_state_1_n = the_state_1_n + 2*(Lp + 1)^(i - 1)
        the_state_1_s = string(the_state_1_n, base = Lp + 1, pad = N)
        index_1, v1 = findall(x -> x == the_state_1_s, states)[1], sqrt(i_num + 1)*sqrt(i_num + 2)
        # println("cr 2 times, $the_state_s -> $the_state_1_s, ($index_1, $label) = $v1") # debug
    end
    # 湮灭两次
    if i_num >= 2
        the_state_2_n = parse(Int, the_state_s, base = Lp + 1)
        the_state_2_n = the_state_2_n - 2*(Lp + 1)^(i - 1)
        the_state_2_s = string(the_state_2_n, base = Lp + 1, pad = N)
        index_2, v2 = findall(x -> x == the_state_2_s, states)[1], sqrt(i_num)*sqrt(i_num - 1)
        # println("an 2 times, $the_state_s -> $the_state_2_s, ($index_2, $label) = $v2") # debug
    end
    index_1 , v1, index_2, v2
end

# 根据周期性边界条件调整临近作用坐标
@memoize function adjust_index_for_pbc(i::Int, j::Int, max::Int)
    if i > max
        # print("adjust_index_for_pbc, ($i, $j) -> ") #debug
        i = 1
    end
    if j > max
        # print("adjust_index_for_pbc, ($i, $j) -> ") #debug
        j = 1
    end
    # println("($i, $j)") #debug
    i, j
end

# 求特征向量的粒子数
function give_the_num_of_particles_of_eigvector(vec::Vector{Float64}, the_model::BH_Model)
    probability = vec.^2
    nhil = length(the_model.states)
    num_of_p = Array{Float64}(undef, nhil)
    for i in 1:nhil
        num_of_p[i] = total_num_particle(the_model.states[i], the_model.N, the_model.Lp)
    end
    sum(probability .* num_of_p)
end