# Proposal: "Apparent complexity of quantized diffusion"

A short note (8–10 pages) on what the black-and-white frames of
`main21D.py` are, with three theorems that are elementary but exact, and
which together explain every feature of the pictures: the rise of
apparent complexity, its peak time, its collapse, the symmetric island
lattices, and why refinement destroys the texture. Numerical checks are
in `phase_theorem.py` and `results.txt`.

Honest positioning: this is not a new phenomenon in nature and not a
new automaton. It is a clean, self-contained piece of mathematics about
*observables* of a diffusing field, with a cautionary message for anyone
who measures "complexity" of a thresholded simulation. arXiv: math.NA
primary, cross-list math.HO and cs.GR (sampling/aliasing). Journal
targets: The Mathematical Intelligencer, European Journal of Physics,
Mathematics Magazine, Journal of Mathematics and the Arts.

## Setting

`u` solves `u_t = κ Δu` on the box `Ω = (0,L)^n` with homogeneous Neumann
data. Quantum `1/q`, phase `θ ∈ [0,1)`. The observed pattern is

    b_θ(x,t) = floor(q u(x,t) + θ) mod 2 .

The grid version samples `u` at the `N^n` cell centres of the explicit
scheme with reflected ghost values (the ternary-cube rule with `g = 0`),
under the CFL condition `Σ_m r_m ≤ 1/2`. `E` = number of grid edges,
`Δ_e u` = jump of `u` across edge `e`. Observables:

- band count `B_q(t)` = number of levels `(k−θ)/q` in `[min u, max u]`;
- total variation `TV(t) = ∫_Ω |∇u| dx`, discrete `TV_h = Σ_e |Δ_e u|`;
- parity perimeter `P_h` = number of edges whose two cells differ in
  colour; crossing count `C_h = Σ_e |floor(q u_a+θ) − floor(q u_b+θ)|`;
- block entropy `H_{q,h}` of the 2×2 patterns of `b_θ`.

## Theorem A (continuum: what is monotone)

1. `B_q(t)` is nonincreasing for every `θ`. *(Maximum principle.)*
2. `TV(t)` is nonincreasing. *Proof.* Even reflection across the faces
   turns `u` into a `2L`-periodic solution on the torus; its gradient
   components are themselves periodic heat solutions, hence
   `∇u(t) = K_{t−s} * ∇u(s)` with `K ≥ 0`, `∫K = 1`; convexity of the
   Euclidean norm gives `|∇u(t)| ≤ K * |∇u(s)|` pointwise; integrate.
   No Kato inequality, no Bakry–Émery; four lines.
3. **Phase-averaged isotherm length equals total variation exactly:**

       ∫_0^1  Σ_k  H^{n−1}({ u(·,t) = (k−θ)/q })  dθ  =  q · TV(t)

   by the coarea formula, hence it is nonincreasing. For a fixed phase the
   left side is a Riemann sum of the coarea integral, so for fixed `θ`
   it is within `o(1)` (as `q → ∞`) of a nonincreasing function.

Hence in the continuum the quantized contour map has *no* rise of
complexity, for any quantum, in the phase-averaged sense.

## Theorem B (grid: the exact perimeter formula)

1. **Discrete TVD.** For each axis `m`, `Σ_i |u_{i+e_m} − u_i|` is
   nonincreasing under the scheme. *Proof.* On the reflected periodic
   grid the scheme is a convex combination of translates (CFL); the
   axis differences satisfy the same scheme; a convex combination of
   translates is an `ℓ¹`-contraction. (Harten's TVD argument.)
2. **Random-phase identity.** With `θ` uniform on `[0,1)`,

       E_θ C_h  =  q · TV_h            (exactly),
       E_θ P_h  =  Σ_e φ( q |Δ_e u| ),  φ(x) = dist(x, 2Z)   (exactly),

   where `φ` is the triangle wave: `φ(x) = x` on `[0,1]`, `2 − x` on
   `[1,2]`, period 2. *Proof.* Under a uniform shift the number of
   quantization levels inside an interval of length `ℓ` is `floor(qℓ)`
   with probability `1 − {qℓ}` and `floor(qℓ)+1` with probability `{qℓ}`;
   the two cells differ in colour iff that number is odd; compute.
3. **Consequences.**
   - Resolved regime: if every edge has `q|Δ_e u| ≤ 1`, then
     `E_θ P_h = q · TV_h`, which is nonincreasing by (1).
   - Always `E_θ P_h ≤ min(q · TV_h, E/1)` since `φ(x) ≤ min(x, 1)`.
   - Every increase of `E_θ P_h` in time is produced by edges whose jump
     exceeds one quantum and is *decreasing through a folded branch of
     φ*. In words: the rise of apparent complexity is aliasing, and
     nothing else. The quantity the eye reads as "complexity being
     created" is `φ` unfolding while `TV_h` shrinks.
   - **Peak criterion.** `E_θ P_h` is nonincreasing from the first time
     `t*` at which `max_e q|Δ_e u| ≤ 1`. On the Neumann box the jumps
     decay like `e^{−λ₁ t}`, `λ₁ = κπ²/L²`, so `t* ≈ λ₁^{−1} log(q·c·h)`
     with `c` set by the initial data.

Numerical check (`phase_theorem.py`, 64×64, `q = 10`, Monte-Carlo over
200 phases vs the formula):

| step | q·TV_h | P(θ=0) | mean_θ P | formula | max q·jump |
|-----:|-------:|-------:|---------:|--------:|-----------:|
|   20 | 95999.7 | 1136 | 1256.2 | 1257.3 | 155.95 |
|  400 | 59437.5 | 4032 | 4004.1 | 4004.7 |  22.19 |
| 1000 | 22542.9 | 4232 | 4285.2 | 4286.4 |   5.75 |
| 2000 |  4435.4 | 4432 | 4436.8 | 4435.4 |   0.92 |
| 3000 |   870.0 |  896 |  867.3 |  870.0 |   0.17 |
| 4000 |   170.5 |  184 |  172.9 |  170.5 |   0.03 |

The formula matches to Monte-Carlo noise; `TV_h` never increased
(asserted in the script); the perimeter peaks at the step where the
largest jump crosses one quantum (0.92 at step 2000); from then on
`P = q·TV_h` to within rounding, exactly as (3) predicts.

## Theorem C (entropy: the limits do not commute)

For `t > 0`:

1. For fixed `q`, `lim_{h→0} H_{q,h}(t) = 0`: the fraction of
   non-uniform 2×2 blocks is `O(h · Λ_q(t))` where `Λ_q` is the total
   isotherm length, so the pattern is asymptotically all uniform blocks.
2. For fixed `h`, if `1, u(x_1,t), …, u(x_m,t)` are linearly independent
   over `Q`, then `(q u(x_j,t) mod 2)_j` is equidistributed on `(R/2Z)^m`
   as `q → ∞` (Weyl), so the Cesàro limit over `q` of `H_{q,h}(t)` is the
   maximum, `log₂ 16 = 4` bits.
3. Remark: symmetric initial data (as in the experiment, a centred
   square) violate the independence in (2): `u(x_j) = u(σ x_j)` for the
   symmetries `σ` of the square, and the achievable entropy is reduced
   by the orbit structure. This is why the aliased frames are
   square-symmetric "designs" rather than noise.

So "the entropy of the quantized field" is not a property of the field;
it depends on which limit is taken first. The measured entropy of 3.99
bits in the experiment is the regime (2), not a statement about the heat
equation.

## Section on the pictures (moiré)

- Locally linear `u`: `floor(q(a·x) + θ) mod 2` on the lattice is a
  two-colour digital plane; sections are Sturmian words. The striped
  regions of the frames.
- Locally quadratic `u`: a sampled Fresnel zone plate. Undersampled zone
  plates replicate their centre: adjacent cells agree in colour wherever
  `q h ∂_m u` is close to an even integer, so the fake centres lie on
  the curves `{ q h ∂_m u ∈ 2Z }` for `m = 1, 2`. *(To verify with an
  overlay figure; this is the one claim in the note not yet checked.)*
- Refinement: `q` fixed, `h → 0` removes all of this (Theorem C.1, and
  the resolution table in `ASSESSMENT.md`). The texture is a function of
  the single ratio `q h |∇u|`.

## What is deliberately NOT claimed

- Monotonicity of the number of components of the pattern, or of the
  number of local maxima of `u`: false in general in two dimensions
  (see Magnanini, arXiv:1604.00530).
- Monotonicity of the Shannon entropy of the value histogram of `u`:
  false (two-valued data have ~1 bit; spreading creates many levels).
- Any statement about "eddies": the field is analytic with smooth
  nested isotherms for `t > 0`.

## Outline

1. Introduction: the frames, the guess "complexity rises then dies".
2. The observable `b_θ`; what a cellular automaton is and why this is
   not one (two sentences, no polemics).
3. Theorem A with proofs.
4. Theorem B with proofs, the table above, and the peak-time estimate.
5. Theorem C, with the symmetry remark.
6. The pictures: digital planes, zone plates, replica centres,
   refinement test (three figures from `out/`).
7. Discussion: measuring complexity of thresholded simulations; the
   phase-averaging trick as a diagnostic (if the phase-averaged
   perimeter exceeds `q·TV_h`, you are looking at aliasing).

## Where a real research question would start (optional, separate)

- Put the mod *inside* the dynamics (Wolfram's continuous cellular
  automata, NKS ch. 4) on the box with reflecting walls in 2-D and 3-D,
  using the ternary-cube kernel. Class-4 structures interacting with the
  3^n boundary strata have not been studied; modest but genuine.
- Sharp asymptotics of `t*` and of the peak value of `E_θ P_h` for
  Gaussian data on `R^n` (closed form via the error function).
- The critical-point counting problems in Magnanini's survey, if one
  wants a hard open problem in the same landscape.
