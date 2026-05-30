import numpy as np

def particle_swarm(func, lwbnd, upbnd, err_goal=1e-3, pop_size=20, max_it=1000, 
                   verbose=False, deflation=None):
    # --- Initialization
    nparam = len(upbnd)
    lwbnd, upbnd = np.array(lwbnd), np.array(upbnd)
    bound = upbnd - lwbnd
    
    # PSO Hyperparameters from your IDL code
    w, alpha = 1.0, 0.9
    c1, c2 = 1.0, 1.0
    minw, minc1, minc2 = 0.1, 0.2, 0.2
    h = 5
    
    # Initialize swarm
    popul = np.random.uniform(lwbnd, upbnd, (pop_size, nparam))
    vel = np.random.uniform(-1, 1, (pop_size, nparam))
    
    # Evaluate initial population
    fpopul = np.array([func(p) for p in popul])
    fevals = pop_size
    
    # Personal bests
    bestpos = np.copy(popul)
    fbestpos = np.copy(fpopul)
    
    # Global best
    g = np.argmin(fbestpos)
    fbest_history = [fbestpos[g]]
    mean_fitness_past = 1e10
    
    success = 0
    iter_count = 0
    nexceed = 0

    while success == 0 and iter_count < max_it:
        iter_count += 1
        
        # --- Dynamic Parameter Update (Your alpha-reduction logic)
        if iter_count > h:
            if fbest_history[-1] >= fbest_history[-h]:
                w, c1, c2 = w * alpha, c1 * alpha, c2 * alpha
        
        if (nexceed / pop_size) > 0.2:
            w, c1, c2 = w * alpha, c1 * alpha, c2 * alpha

        w = max(w, minw)
        c1 = max(c1, minc1)
        c2 = max(c2, minc2)

        # --- Velocity and Position Update
        r1 = np.random.rand(pop_size, nparam)
        r2 = np.random.rand(pop_size, nparam)
        
        vel = (w * vel + 
               c1 * r1 * (bestpos - popul) + 
               c2 * r2 * (bestpos[g] - popul))
        
        popul += vel
        
        # --- Boundary Check (Random Reset logic from IDL)
        nexceed = 0
        for i in range(pop_size):
            if np.any(popul[i] < lwbnd) or np.any(popul[i] > upbnd):
                popul[i] = np.random.uniform(lwbnd, upbnd)
                vel[i] = np.random.uniform(-0.1, 0.1, nparam)
                nexceed += 1

        # --- Evaluation & Deflation
        fpopul = np.array([func(p) for p in popul])
        fevals += pop_size
        
        if deflation is not None:
            # Penalize areas near known minima
            for d_point in deflation:
                dist_sq = np.sum((popul - d_point)**2, axis=1)
                fpopul = (fpopul + 1e-30) / dist_sq

        # --- Update Bests
        better_mask = fpopul < fbestpos
        fbestpos[better_mask] = fpopul[better_mask]
        bestpos[better_mask] = popul[better_mask]
        
        g = np.argmin(fbestpos)
        fbest_history.append(fbestpos[g])
        
        mean_fitness = np.median(fpopul)
        
        # --- Stopping Criterion
        criterion = abs((mean_fitness - mean_fitness_past) / (mean_fitness_past + 1e-15))
        if criterion <= err_goal:
            success = 1
        
        mean_fitness_past = mean_fitness
        
        if verbose and iter_count % 10 == 0:
            print(f"Iter {iter_count}: Best {fbestpos[g]:.5f}, w={w:.3f}")

    return bestpos[g], fbestpos[g], fevals
