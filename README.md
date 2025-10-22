# Monte-Carlo-MPS

This repository provides a numerical toolkit for simulating the dynamics of open quantum systems. It unravels the master equation using the quantum jumps method and uses the Matrix Product State (MPS) formalism to represent the states.

The algorithm proceeds as follows:
1. Initialize the system in a specified Matrix Product State (MPS) $|\psi(0)\rangle$
2. For each quantum trajectory:
   - **Time Evolution**: Evolve the state under the non-Hermitian effective Hamiltonian:
     $$H_{\text{eff}} = H - \frac{i}{2}\sum_k L_k^\dagger L_k
     $$
     using Time-Evolving Block Decimation (TEBD).
   
   - **Jump Condition**: Monitor the norm of the state $\langle\psi(t)|\psi(t)\rangle$
     - When the $\langle\psi(t)|\psi(t)\rangle$ is below a randomly generated threshold $\epsilon \in [0,1]$:
       - Randomly select a Lindblad operator $L_k$ with probability proportional to $\langle L_k^\dagger L_k\rangle$
       - Apply the quantum jump: $|\psi'\rangle = L_k|\psi\rangle$
       - Renormalize the state: $|\psi(t)\rangle = |\psi'\rangle/\||\psi'\rangle\|$

 3. **Ensemble Averaging**
   - Repeat for $N$ independent trajectories
   - Compute the evolution of observables as the ensemble average:
     
     $$\langle O \rangle = \frac{1}{N}\sum_{i=1}^N \langle\psi_i(t)|O|\psi_i(t)\rangle
     $$
   - Compute the evolution of average bipartite entanglement entropy:
     
     $$S_e(t) = \frac{1}{N}\sum_{i=1}^N S_e(|\psi_i(t)\rangle)
     $$,
     
     where $S_e(|\psi_i(t)\rangle) = -Tr(\rho_A \log \rho_A)$, $\rho_A = Tr_B \ |\psi_i(t)\rangle \langle \psi_i(t)|$, $dim A = dim B$.
