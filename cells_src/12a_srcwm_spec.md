# SRCWM latent-dynamics engine (E4 core)

The smoke-scale engine below makes the model *precise*, not vague. It is a
first-order, linear, action-conditioned dynamics on a cell-type program space:

* **Latent space.** Cell-type centroid profiles $W_c \in \mathbb{R}^{G\times C}$
  from the single-cell reference (Census slice in E3; synthetic here) define the
  gene programs. A drug signature $z_d$ is mapped to a $C$-dim action vector
  $v_d = W_c^{\!\top} z_d$.
* **Action-conditioned dynamics.** Cell types interact through a mixing matrix
  $M$ (row-stochastic, derived from centroid correlations). One latent step is
  $v^{(t+1)} = (1-\eta)\, v^{(t)} + \eta\, M v^{(t)}$; the rollout depth $T$,
  chosen from $\{1,2,3\}$, and the scale-matching constant $\eta$ are the named
  design constants from E0's deliverable list.
* **Lifting ($\Lambda$).** The disease bulk delta $\delta_k$ is attributed to
  cell types by NNLS on $W_c$ (validated in E3), giving the latent target
  $\hat{\delta}_k = W_c^{\!\top} a^{(k)}$.
* **Rescue score.** $s(k,d) = -\cos\!\big(v_d^{(T)}, \hat{\delta}_k\big)$ —
  a drug rescues disease $k$ if its latent action anti-aligns with the lifted
  pathology direction.
* **Shared objective.** Panel-level ranking uses $\min_k$ over the *training*
  diseases only, per LOGO fold — so held-out performance tests transfer of the
  shared-rescue hypothesis, not memorization.

Known-answer validation (this cell, executed): a synthetic world with planted
shared-rescue drugs, single-disease specialists, and nulls; every gate below is
asserted, and the E5-style statistics (paired bootstrap, Holm correction) run on
the synthetic folds to validate the machinery.
