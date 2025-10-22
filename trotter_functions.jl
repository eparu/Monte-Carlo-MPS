function trotter_gates_second_order_xxz(sites, tau, gamma, Δ; J=1)
    N = length(sites)
    gates = ITensor[]
    for j in 1:N-1
        # make h_j,j+1 term
        s1 = sites[j]
        s2 = sites[j + 1]
        hjj_1 = 4*J * (op("Sx", s1) * op("Sx", s2) + 
                    op("Sy", s1) * op("Sy", s2) + 
                     Δ * op("Sz", s1) * op("Sz", s2))  # XXZ interaction
        G = exp(-im * (tau / 2) * hjj_1)  # Half-step evolution
        push!(gates, G)
    end
    # Full-step two-site interaction terms
    for j in 1:N
        s = sites[j]
        hloc = -im/2 * gamma * op("S- * S+", s) # Local collapse operator
        G = exp(-im * tau * hloc)  # Full-step evolution
        push!(gates, G)
    end

    append!(gates, reverse(gates[1:N-1]))
    return gates
end


function trotter_gates_second_order_XYZz(sites, tau, gamma, hz; Jx=1, Jy=0.75, Jz=0.5)
    N = length(sites)
    gates = ITensor[]
    for j in 1:N-1
        # make h_j,j+1 term
        s1 = sites[j]
        s2 = sites[j + 1]
        hjj_1 = 4* (Jx*op("Sx", s1) * op("Sx", s2) + 
                    Jy*op("Sy", s1) * op("Sy", s2) + 
                    Jz*op("Sz", s1) * op("Sz", s2))  # XXZ interaction
        G = exp(-im * (tau / 2) * hjj_1)  # Half-step evolution
        push!(gates, G)
    end
    # Full-step two-site interaction terms
    for j in 1:N
        s = sites[j]
        hloc = 2*hz*op("Sz", s) - im/2 * gamma * op("S- * S+", s) # Local collapse operator
        G = exp(-im * tau * hloc)  # Full-step evolution
        push!(gates, G)
    end

    append!(gates, reverse(gates[1:N-1]))

    return gates
end


function trotter_gates_second_order_XXZz(sites, tau, gamma; Δ=0.5, hz=1.0, J=1.0)
    N = length(sites)
    gates = ITensor[]


    for j in 1:N-1
        # make h_j,j+1 term
        s1 = sites[j]
        s2 = sites[j + 1]
        hjj_1 = 4*J * (op("Sx", s1) * op("Sx", s2) + 
                    op("Sy", s1) * op("Sy", s2) + 
                     Δ * op("Sz", s1) * op("Sz", s2))  # XXZ interaction
        G = exp(-im * (tau / 2) * hjj_1)  # Half-step evolution
        push!(gates, G)
    end
    # Full-step two-site interaction terms
    for j in 1:N
        s = sites[j]
        hloc = 2*hz*op("Sz", s) - im/2 * gamma * op("S- * S+", s) # Local collapse operator
        G = exp(-im * tau * hloc)  # Full-step evolution
        push!(gates, G)
    end

    append!(gates, reverse(gates[1:N-1]))

    return gates
end