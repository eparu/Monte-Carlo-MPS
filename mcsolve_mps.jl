using Distributed
using ITensors, ITensorMPS
using Random
using LinearAlgebra
using NPZ

function entanglement_entropy(psi)
    # Compute the reduced density matrix of the left half of the system
    b = div(length(psi), 2)
    psi = orthogonalize(psi, b)
    #SVD
    U,S,V = svd(psi[b], (linkinds(psi, b-1)..., siteinds(psi, b)...))
    s = diag(S)[diag(S).>1e-10]
    entropy = -sum(s.^2 .* log2.(s.^2))
    return entropy
end




function MC_MPS_evolution(
    sites, 
    make_trotter_gates,     # Function to create Trotter gates
    t_list, 
    initial_state,          # MPS initial state
    gamma;                  # Collapse operator: we consider that collapse operators have the same form but applied to different sites
    seed=42,
    n_traj=10, 
    cutoff=1e-9, 
    Dmax=100, 
    save=true,
    folder="data/"
)
    
    # Number of spins
    N = length(sites)
    # Calculate time step from the time list
    t_step = t_list[2] - t_list[1]
    # Trotter gates
    gates = make_trotter_gates(sites, t_step, gamma)

    # Initialize the random number generator and seeds
    rng = Random.MersenneTwister(seed)
    seeds = map(i -> rand(rng, UInt64), 1:n_traj)

    # Define the function to execute for each trajectory
    function evolve_trajectory(traj)
        ee_list = zeros(length(t_list))
        z_list = zeros(length(t_list))
        bond_dimension_list = zeros(length(t_list))

        println("Processing trajectory: ", traj)
        rng_traj = Random.MersenneTwister(seeds[traj])

        norm_prev = 1.0
        psi = initial_state

        r_1 = rand(rng_traj)
        

        for (t_idx, t) in enumerate(t_list)
            
            #save entanglement entropy, z expectation value and bond dimension
            psi_norm = normalize(psi)
            ee_list[t_idx] = entanglement_entropy(psi_norm)
            z_list[t_idx] = expect(psi, "Sz"; sites=1)
            bond_dimension_list[t_idx] = maximum(linkdims(psi))

            # Apply Trotter gates using the provided function
            psi_prev = psi
            psi = apply(gates, psi; cutoff=cutoff, maxdim=Dmax)

            norm_ = norm(psi)^2
            if norm_ < r_1

                # Calculate tau for norm matching
                delta_norm = norm_prev - norm_
                tau = t_step * (norm_prev - r_1) / delta_norm

                # Generate gates for the fractional step `tau`
                gates_tau = make_trotter_gates(sites, tau, gamma)
                psi_jump = apply(gates_tau, psi_prev; cutoff=cutoff, maxdim=Dmax)

                # Compute probabilities for jump
                r_2 = rand(rng_traj)
                probabilities = [expect(psi_jump, "S- * S+"; sites=j) for j in 1:N]
                sum_probs = sum(probabilities)

                if sum_probs > 0
                    probabilities /= sum_probs
                    cum_sums = cumsum(probabilities)

                    # Select the jump operator
                    n = findfirst(x -> x >= r_2, cum_sums)
                    if n !== nothing
                        L_n = sqrt(gamma) * op("S+", sites[n])
                        psi = apply(L_n, psi_jump; cutoff=cutoff, maxdim=Dmax)
                        psi = normalize(psi)
                    end
                end

                # Continue evolution for the remainder of the time step
                gates_remainder = make_trotter_gates(sites, t_step - tau, gamma)
                psi = apply(gates_remainder, psi; cutoff=cutoff, maxdim=Dmax)

                # Update `r_1`
                r_1 = rand(rng_traj)
            end
            norm_prev = norm_
        end
        if save
            # Save data for this trajectory
            npzwrite("$(folder)/ee_list_traj$(traj).npz", ee_list)
            npzwrite("$(folder)/z_list_traj$(traj).npz", z_list)
            npzwrite("$(folder)/bond_dim_traj$(traj).npz", bond_dimension_list)
        else
            return ee_list, z_list, bond_dimension_list
        end
    end

    # Use `pmap` to execute the function in parallel
    results = map(evolve_trajectory, 1:n_traj)
    if save
        return nothing
    else
        ee_array = zeros(n_traj, length(t_list))
        z_array = zeros(n_traj, length(t_list))
        bond_array = zeros(n_traj, length(t_list))
        for (i, (ee_list, z_list, bond_dimension_list)) in enumerate(results)
            ee_array[i, :] = ee_list
            z_array[i, :] = z_list
            bond_array[i, :] = bond_dimension_list
        end
        return ee_array, z_array, bond_array
    end
end


function run_instance_with_pmap(N, gamma, n_traj, t_list, trotter_function; folder_name="XXZ", Dmax=100, seed=42)
    sites = siteinds("Qubit", N)
    print("run_instance")
    if !isdir("data/$(folder_name)")
        mkdir("data/$(folder_name)")
    end
    if !isdir("data/$(folder_name)/N$(N)_gamma$(gamma)_ntraj$(n_traj)_tf$(t_list[end])")
        mkdir("data/$(folder_name)/N$(N)_gamma$(gamma)_ntraj$(n_traj)_tf$(t_list[end])")
    end
    repo = "data/$(folder_name)/N$(N)_gamma$(gamma)_ntraj$(n_traj)_tf$(t_list[end])"
    print(repo)
    #initial state
    psi0 = productMPS(sites, vcat(["1"], ["0" for i in 1:N-1]))

    # Run the Monte Carlo simulation
    MC_MPS_evolution_pmap(
            sites, 
            trotter_function, 
            t_list, 
            psi0, 
            gamma; 
            seed=seed,
            n_traj=n_traj, 
            Dmax=Dmax, 
            save=true,
            folder=repo
        )
end

