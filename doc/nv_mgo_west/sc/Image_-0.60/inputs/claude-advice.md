I don't think the algorithm explains -0.60. In WEST TDDFT, the cost of one iteration is fixed by things that are essentially the same in -0.60 and -0.20. Something in the environment slowed that one run down.

What one wbse.x Davidson iteration costs

Each iteration applies the Liouvillian to the new correction vectors, one per root that hasn't converged yet. That's WBSE_dim in the output. For each vector, wbse.x does three things:

1. Applies H_KS − ε to every occupied orbital's component, including the ACE exchange operator built from n_exx_lowrank: 1294 orbitals. The cost depends on system size, ecut and number of occupied bands, which are identical across images.
2. Hartree and xc kernel response. These are FFTs on a fixed grid, so the cost is the same everywhere.
3. Exact-exchange kernel. This uses the Coulomb integrals precomputed by wbse_init.x, only for Wannier pairs with overlap above overlap_thr: 0.001. The pair count is the one quantity that could change with geometry.

The subspace diagonalization (Diago_dim up to 48) takes under 0.1 s.

So the time per iteration should scale roughly as WBSE_dim × the pair count. Here's how the two runs compare:

┌──────────────────────────────────┬────────────────┬────────────┐
│                                  │     -0.20      │   -0.60    │
├──────────────────────────────────┼────────────────┼────────────┤
│ WBSE_dim (vectors per iteration) │ 6              │ 6          │
├──────────────────────────────────┼────────────────┼────────────┤
│ Wannier pairs (nglob)            │ 7924           │ 7904       │
├──────────────────────────────────┼────────────────┼────────────┤
│ Wannier spread                   │ 2492.95        │ 2492.82    │
├──────────────────────────────────┼────────────────┼────────────┤
│ Memory estimate                  │ 1643.66 Mb     │ 1643.66 Mb │
├──────────────────────────────────┼────────────────┼────────────┤
│ wbse_init wall time              │ 27m56s         │ 28m12s     │
├──────────────────────────────────┼────────────────┼────────────┤
│ "starting" iteration             │ 2m32s          │ 25m37s     │
├──────────────────────────────────┼────────────────┼────────────┤
│ Later iterations                 │ 35 s – 1.5 min │ 8 – 57 min │
└──────────────────────────────────┴────────────────┴────────────┘

The work per iteration is the same, but -0.60 takes about 20× longer. The only real per-iteration variation inside WEST is visible in -0.20 itself: once two roots converge, WBSE_dim drops from 6 to 4 and iterations get shorter (iteration 35 took 35 s).

Two separate effects in your runs

1. Number of iterations (physics, expected). Near Q = 0 the runs need about 50 iterations. At |Q| ≥ 0.5 they need 100–260 (+0.60: 221, +0.80: 263, +1.40: 242). Along the configuration coordinate the excited states get closer together and start to cross, so Davidson converges 6 roots to trev_liouville = 1e-6 more slowly. Those runs still take only 1–2 min per iteration.

2. Time per iteration in -0.60 (not physics). The evidence points to the environment:
- In the same job, on the same 4 nodes, pw.x (about 2.1 min per SCF step, same as -0.20) and wbse_init.x ran at normal speed. Only wbse.x was slow, from its very first step.
- The iteration times jump around erratically (25, 22, 24, … 8, 16, 9, 7.6, 20, 52, 57 min). Compute cost that's tied to the physics would be steady. Random swings like this usually come from a degraded GPU or network card on one node, or from interconnect or filesystem contention.
- Image +1.40 was running at the same time on other nodes at a normal 1.6 min per iteration.

I can't name the exact cause. The job hit the 10 h limit after 15 iterations, so WEST never printed its per-routine timing breakdown. Its nodes were nid001116, nid001117, nid001785 and nid001788. No restart was saved, because that only happens every n_steps_write_restart: 50 iterations.

Suggested next step

Resubmit only the wbse.x step for -0.60. The pw, nscf and wbse_init outputs are fine to reuse.
- Delete the leftover pwscf.wfc* files first.
- Add #SBATCH --exclude=nid001116,nid001117,nid001785,nid001788.
- Set n_steps_write_restart to about 10 so a timeout doesn't lose everything.

If the first iteration then takes about 2 min, the old nodes were the problem. If it's still about 25 min on new nodes, the slowdown really does depend on the geometry, and the timing report from a clean finish would show which routine is responsible.

✻ Baked for 21s · done 3:34 PM

※ recap: We're working out why your WBSE runs differ: large-displacement images need more iterations, while -0.60 alone runs about 20× slower per iteration on the same workload. Next, resubmit only wbse.x for -0.60, excluding its old nodes. (disable recaps in /config)

❯ so lets say in one calculation energy eigen values are very degenerate or nearly degenerate. does this mean it each convergence iteration would talk longer, or each iteration should take the same amount of time but more time steps:

Near-degeneracy gives you more iterations. It doesn't make each iteration take longer.

Why the time per iteration doesn't change: each Davidson iteration applies the Liouvillian once to each unconverged correction vector (WBSE_dim). Applying it to one vector costs the same no matter what the eigenvalues are. That cost depends on the system size, cutoff, number of occupied bands, ACE rank (n_exx_lowrank) and the Wannier pair count. The subspace diagonalization is tiny (under 0.1 s for Diago_dim ≤ 48), so the spectrum barely affects it.

Why there are more iterations: how fast Davidson converges depends on how well separated the wanted eigenvalues are from each other and from the next root you didn't ask for. Nearly degenerate roots are hard to separate, so the residuals shrink slowly. Two things make it worse:
- If roots 6 and 7 are close, root 6 converges especially slowly, because the solver keeps mixing in the unwanted 7th state.
- Each time the subspace hits its maximum size (Diago_dim 48) it restarts back to 12 and throws away some information. Slow convergence means more restarts.

One indirect effect on time per iteration: the cost scales with how many roots are still unconverged. In a well-behaved run, WBSE_dim drops from 6 to 4 and so on as roots converge, and iterations get cheaper. You can see this in -0.20 from iteration 34 on, where one iteration took 35 s. With near-degeneracy, roots stay unconverged longer, so the run stays at the full cost of 6 vectors for more iterations. But a run with all 6 roots unconverged is the most an iteration can ever cost, and that's about 1–2 min for your system.

So degeneracy can turn 50 iterations into 250, as in your +0.80 and +1.40 runs. It can't turn a 1-minute iteration into a 25-minute one. -0.60 had WBSE_dim = 6, the same as every other image, so its slow iterations must come from the run environment rather than the eigenvalue spectrum.



You should rerun -0.60, but as a fresh wbse.x run. A Davidson restart (wbse_calculation: d) isn't possible.

Why a restart won't work: d needs a checkpoint, and none was written. With n_steps_write_restart: 50, the first checkpoint would have come at iteration 50, but the job was killed at iteration 15. Image_-0.60 has no restart folder, and west.wbse.save/ holds only wbse.json.

Why that costs you almost nothing: at a normal 1–2 min per iteration, those 15 iterations are about 20 minutes of work. All the expensive earlier steps finished cleanly and can be reused:
- pwscf.save, including the ACE files from nscf
- west.wbse_init.save/, the Coulomb integrals for about 7900 Wannier pairs

So rerun only the last step.

What to change before resubmitting:
1. job_submit.sh: comment out the pw.x, cp/rm/mv and wbse_init.x run_cmd lines, and keep only the wbse.x line. Don't leave the rm -rf pwscf.save pwscf.wfc* line active, because it would delete the ground state you need.
2. Add #SBATCH --exclude=nid001116,nid001117,nid001785,nid001788 so the job avoids the nodes the slow run used.
3. wbse.in: keep wbse_calculation: D and set n_steps_write_restart: 10. If this run gets cut off, you can then resubmit with d, using the same task count and -nb 16.
4. Delete the leftover pwscf.wfc1–pwscf.wfc16 in Image_-0.60 (16 × 1.5 GB). wbse.x writes these when it starts, and they're stale from the killed run.

Expect about 150–250 iterations, similar to +0.60 (221) and +0.80 (263). At about 1 min per iteration that's 2.5–4 h, well within the 10 h limit. Check the "starting" iteration time in the new wbse.out. If it's about 2 min, the old nodes were the problem. If it's still about 25 min, something specific to this geometry is causing it, and we'd look further.

I can make the edits to job_submit.sh and wbse.in and delete the stale wfc files if you'd like. I'll leave submitting the job to you.

✻ Brewed for 19s · done 3:42 PM