# -*- coding: utf-8 -*-
"""
Particle Swarms -Boids function optimization

Python Implementation DescriptionCore Logic: 
    
A high-performance Python port of the Kennedy & Eberhart (1995) 
Particle Swarm Optimization, optimized via NumPy vectorization to
eliminate the overhead of IDL-style loops.

Dynamic Adaptation: Features a dynamic inertia weight ($w$) and
cognitive/social coefficient ($c_1, c_2$) reduction scheme, triggered
when the swarm fails to improve or when boundary exceedance is high.

Boids Integration: Includes an optional Bio-Inspired Steering mode
based on Craig Reynolds’ (1987) model. When enabled, it applies local
Separation, Alignment, and Cohesion forces to the velocity vector.

Diversity Preservation: Unlike standard PSO which can collapse into
local optima, the Boids-enhanced version uses a "short-range repulsion
sphere" to maintain swarm diversity and improve global exploration.

Robustness: Implements "Random Re-entry" for particles exceeding boundaries
and supports Function Deflation to steer the search away from previously
discovered or undesired coordinate sets.

 Input
    ErrGoal = desired accuracy (float or double)
    func    = name of the function (string)
    upbnd   = upper bound (float or double), array
    lwbnd   = lower bound (float or double), array
  
 Output
    success = 0 (unsuccessful), 1 (successful)
    xmin    = best parameters (array)
    fxmin   = evaluation (merit, fitness) at the best position
    fevals  = number of function evaluations

 Algorithm
    Particle Swarm Optimization, dynamic inertia reduction variant

    In the event that a particle falls outside the allowed domain, it is
    given a new position randomly chosen inside the parameter domain

 Reference

  Kennedy, J.; Eberhart, R. (1995). "Particle Swarm Optimization".
  Proceedings of IEEE International Conference on Neural Networks.
  IV. pp. 1942–1948.

  Clerc, M. (2012). "Standard Particle Swarm Optimisation" (PDF).
  HAL open access archive.  
  https://hal.archives-ouvertes.fr/file/index/docid/764996
  /filename/SPSO_descriptions.pdf
  
  Craig W. Reynolds 
  "Flocks, Herds, and Schools: A Distributed Behavioral Model""
  ACM SIGGRAPH Computer Graphics, Volume 21, Issue 4, July 1987, pp. 25–34.
  DOI: 10.1145/37401.37406
  
  Couzin, I. D., et al. (2002). Collective Memory and Spatial 
  Self-Organization in Animal Groups. Journal of Theoretical Biology.

  Vicsek, T., et al. (1995). Novel Type of Phase Transition in a System of
  Self-Driven Particles. Physical Review Letters.

 Author
    Wing-Fai Thi

 Version 1 22/3/2018 based on an old version from 2007
 Version 2 23/4/2026 Translated to Python with the help of Gemina AI and
                      add the Boids model
"""               
import numpy as np
import scipy.stats.qmc as qmc
import time


def particle_swarm(err_goal, func, lwbnd, upbnd, pop_size=32, max_it=1000, 
                   deflation=None, verbose=False, boids=False, print_results=False,
                   w_sep_init=0.5, w_align_init=0.1, w_coh_init=0.1, 
                   coverage_threshold=0.4, min_it_param=25, quasi=True):
    """
    Highly optimized PSO with vectorized Boids logic and Area Coverage switching.
    
    Logic:
    1. Forced Exploration: Boids stay active until min_it is reached.
    2. Conditional Exploration: Boids stay active until coverage_threshold is met.
    3. High-Precision Descent: Boids disable, allowing particles to collapse.
    """
    lwbnd, upbnd = np.array(lwbnd), np.array(upbnd)
    nparam = len(upbnd)
    # Scaling the exploration floor by dimensionality
    min_it = min_it_param * nparam 
    bound = upbnd - lwbnd
    
    # PSO Parameters
    max_w, min_w = 1.0, 0.1
    c1, c2 = 1.0, 1.0
    w, alpha = max_w, 0.9
    
    # Boids Parameters
    r_sep = np.linalg.norm(bound) * 0.05 
    r_neigh = np.linalg.norm(bound) * 0.2
    
    # Initialization
    if ('quasi'):
        # 1. Initialize Sobol Engine
        sampler = qmc.Sobol(d=nparam, scramble=True)
    
        # 2. Generate Quasi-Random Population
        # Sobol generates values in [0, 1]. We scale them to our bounds.
        # Using pop_size=16 or 32 is recommended to avoid the UserWarning
        sample = sampler.random(n=pop_size) # Shape (pop_size, nparam)
        popul = qmc.scale(sample, lwbnd, upbnd).T # Transpose to (nparam, pop_size)
    
        # 3. Initialize Velocity (Using Sobol or Standard Normal)
        # Usually, velocity is still initialized as a small random nudge
        vel = (sampler.random(n=pop_size).T - 0.5) * 0.1 * bound[:, np.newaxis]
    else:
        popul = np.random.rand(nparam, pop_size) * bound[:, np.newaxis] + lwbnd[:, np.newaxis]
        vel = np.random.randn(nparam, pop_size) * 0.1

    fpopul = np.apply_along_axis(func, 0, popul)
    
    best_pos = np.copy(popul)
    f_best_pos = np.copy(fpopul)
    g_idx = np.argmin(f_best_pos)
    f_best_part = f_best_pos[g_idx]
    
    fevals = pop_size
    iter_count = 0
    success = 0
    stagnation_counter = 0

    # Grid Coverage Tracking
    grid_res = 20 
    visited_cells = set()
    total_cells = grid_res ** nparam
    
    # SAFETY FIX 1: Cap coverage tracking for high dimensions (no perf impact on 2D)
    MAX_CELLS_TRACKED = 100000 if nparam <= 5 else 10000
    coverage_capped = (total_cells > MAX_CELLS_TRACKED)
    if coverage_capped:
        total_cells = MAX_CELLS_TRACKED
    
    boids_active = boids

    while (success == 0) and (iter_count < max_it):
        iter_count += 1
        
        # 1. Boids Component (Vectorized)
        v_boids = np.zeros_like(vel)
        if boids_active:
            diff_tensor = popul.T[:, np.newaxis, :] - popul.T[np.newaxis, :, :]
            dist_matrix = np.linalg.norm(diff_tensor, axis=2)
            decay = 1.0 - (iter_count / max_it)
            
            # Separation, Cohesion, Alignment
            sep_mask = (dist_matrix > 0) & (dist_matrix < r_sep)
            v_sep = np.sum(sep_mask[:, :, np.newaxis] * diff_tensor, axis=1).T
            
            neigh_mask = (dist_matrix > 0) & (dist_matrix < r_neigh)
            neigh_counts = np.sum(neigh_mask, axis=1, keepdims=True) + 1e-9
            v_coh = (np.sum(neigh_mask[:, :, np.newaxis] * -diff_tensor, axis=1) / neigh_counts).T
            v_align = (neigh_mask @ vel.T / neigh_counts).T 
            
            v_boids = (v_sep * w_sep_init + v_coh * w_coh_init + v_align * w_align_init) * decay

        # 2. Velocity & Position Update
        r1, r2 = np.random.rand(nparam, pop_size), np.random.rand(nparam, pop_size)
        vel = (w * vel + 
               c1 * r1 * (best_pos - popul) + 
               c2 * r2 * (best_pos[:, [g_idx]] - popul) + 
               v_boids)
        popul += vel

        # 3. Boundary Handling
        out_mask = np.any((popul < lwbnd[:, None]) | (popul > upbnd[:, None]), axis=0)
        if np.any(out_mask):
            popul[:, out_mask] = np.random.rand(nparam, np.sum(out_mask)) * bound[:, None] + lwbnd[:, None]
            vel[:, out_mask] = 0

        # 4. Evaluation & Deflation
        fpopul = np.apply_along_axis(func, 0, popul)
        
        # SAFETY FIX 2: Prevent division by zero in deflation (critical for NaN prevention)
        if deflation is not None:
            for d_pos in np.atleast_2d(deflation):
                dists_sq = np.sum((popul - d_pos[:, None])**2, axis=0)
                dists_sq = np.maximum(dists_sq, 1e-20)  # ← Only change: prevents Inf/NaN
                fpopul = (fpopul + 1e-30) / dists_sq
        fevals += pop_size

        # 5. Coverage Monitoring
        # SAFETY FIX 3: Only track if not capped (no perf impact on 2D)
        if not coverage_capped or len(visited_cells) < MAX_CELLS_TRACKED:
            for i in range(pop_size):
                norm_pos = (popul[:, i] - lwbnd) / (bound + 1e-30)
                grid_idx = tuple(np.clip((norm_pos * (grid_res - 1)).astype(int), 0, grid_res - 1))
                visited_cells.add(grid_idx)
        
        current_coverage = len(visited_cells) / total_cells
        if boids_active and (current_coverage >= coverage_threshold) and (iter_count >= min_it):
            boids_active = False
            if verbose: print(f"Iter {iter_count}: Coverage {current_coverage:.2%} reached. Boids disabled.")

        # 6. Best Tracking & Stagnation (YOUR LOGIC - KEEP IT!)
        improved = fpopul < f_best_pos
        f_best_pos[improved] = fpopul[improved]
        best_pos[:, improved] = popul[:, improved]
        new_g_idx = np.argmin(f_best_pos)
        
        if f_best_pos[new_g_idx] < f_best_part:
            if abs(f_best_part - f_best_pos[new_g_idx]) < err_goal: success = 1
            f_best_part, g_idx, stagnation_counter = f_best_pos[new_g_idx], new_g_idx, 0
        else:
            stagnation_counter += 1

        if stagnation_counter > 10:
            w = max(w * alpha, min_w)
            
        # Add this inside your loop for expensive runs:
        if verbose and iter_count % 5 == 0:
            print(f"[{iter_count}] Best: {f_best_part:.6e} | Coverage: {current_coverage:.1%}")
        
    if print_results:
        print(f"Best: {f_best_part:.4e} | Evals: {fevals} | Success: {success}")
        
    return success, best_pos[:, g_idx], f_best_part, fevals

# --- TEST SUITE ---

def run_benchmarks():
    # Test Functions
    def rastrigin(x):
        return 10 * len(x) + np.sum(x**2 - 10 * np.cos(2 * np.pi * x))

    def rosenbrock(x):
        return np.sum(100 * (x[1:] - x[:-1]**2)**2 + (1 - x[:-1])**2)

    def ackley(x):
        arg1 = -0.2 * np.sqrt(0.5 * np.sum(x**2))
        arg2 = 0.5 * np.sum(np.cos(2 * np.pi * x))
        return -20 * np.exp(arg1) - np.exp(arg2) + 20 + np.e
    
    def step_func(x):
        return np.sum(np.floor(x))
    
    def griewank_2d(x):
        """Griewank function (2-D).
        A multimodal function with a high number of local minima and
        a single global minimum.

        This is the Griewank Function (2-D or 10-D)
        Bound: X(i)=[-600,600], for i=1,2,...,10
        Global Optimum: 0, at origin

        Args:
            x:  A list or array of size 2 representing the input variables.

        Returns:
            The value of the Griewank function at the given point.
        """
        d = 200.0
        u1 = 0.0
        u2 = 1.0
        for j in range(2):
            u1 += x[j] * x[j] / d
            u2 *= np.cos(x[j] / np.sqrt(j + 1))
        return u1 - u2 + 1.0


    tests = [
        {"name": "Rastrigin (2D)", "func": rastrigin, "lw": [-5.12]*2, "up": [5.12]*2, "goal": 1e-6},
        {"name": "Rosenbrock (2D)", "func": rosenbrock, "lw": [-2.048]*2, "up": [2.048]*2, "goal": 1e-8},
        {"name": "Ackley (2D)", "func": ackley, "lw": [-32.768]*2, "up": [32.768]*2, "goal": 1e-8},
        {"name": "Step Function (5D)", "func": step_func, "lw": [-5.12]*5, "up": [5.12]*5, "goal": -25},
        {"name": "Griewank (2D)", "func": griewank_2d, "lw": [-600.]*5, "up": [600.]*5, "goal": 1e-8}
    ]

    for t in tests:
        print(f"\n--- {t['name']} ---")
        for mode in [False, True]:
            start = time.time()
            # Using pop_size=32 to satisfy Sobol power-of-2 and ensure coverage
            _, _, val, evals = particle_swarm(
                t['goal'], t['func'], t['lw'], t['up'], 
                boids=mode, 
                pop_size=32,
                verbose=False # Set to True for expensive-run monitoring
            )
            elapsed = time.time() - start
            mode_str = "Boids-Decay" if mode else "Standard PSO"
            print(f"{mode_str:12}: {val:.4e} in {evals} evals ({elapsed:.3f}s)")

def run_expensive_experiment(sim, lw, up):
    print("Initializing Sobol-Boids Swarm for High-Stakes Optimization...")
    
    # 16 particles (Power of 2) is often the 'Sweet Spot' for 
    # expensive 2D/3D functions using Sobol.
    success, best_coords, best_val, total_evals = particle_swarm(
        err_goal=1e-8, 
        func=sim, 
        lwbnd=lw, upbnd=up, 
        pop_size=16, 
        boids=True,
        verbose=True, # Shows the Iteration/Coverage logs
        coverage_threshold=0.3, # Exit Boids early to save money
        min_it_param=15 # Scale for efficiency
    )
    
    print("-" * 30)
    print(f"Optimization Complete!")
    print(f"Final Sim Result: {best_val:.8e}")
    print(f"Total Compute Cost: {total_evals} evaluations")
    
    
if __name__ == "__main__":
    run_benchmarks()

