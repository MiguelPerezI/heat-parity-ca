# Is the black-and-white pattern of `main21D.py` a new cellular automaton?

Assessment written 2026-10-22, based on the script, the screenshots of
18:27–18:31, and the simulations in `analyze.py` / `restest.py`
(results in `results.txt`, images in `out/`).

## 1. What the pattern actually is

`main21D.py` integrates the 2-D heat equation on a 64×64 grid with
real-valued cells. The display then colours each cell by

    b(i,j,t) = floor(10 · u(i,j,t)) mod 2        (white if 1, black if 0)

So the picture is **not** the state of the simulation. It is the contour
map of the temperature with isotherms every 0.1 degrees, with the bands
between successive isotherms filled alternately black and white. The 3-D
height of each facet is still `u`; only the colour is the parity.

Consequences:

- The underlying dynamics is the discrete heat equation: linear,
  real-valued, strictly smoothing. In the cellular-automaton literature
  this is the standard "continuous-valued CA" or coupled map lattice,
  and it is the *least* complex such system (every mode decays).
- The binary picture is an **observable**, not a state. Two temperature
  fields with the same parity pattern evolve to different parity
  patterns, so there is no local rule on black/white cells that
  reproduces the pictures. By the standard definition (von Neumann,
  Hedlund, Wolfram: a map on configurations over a finite alphabet
  determined by a local rule) **it is not a cellular automaton**, new
  or otherwise.

## 2. Why it looks like a CA, and why it is beautiful: aliasing

Each band is 0.1 degrees wide in temperature, hence `0.1/|∇u|` wide in
space. Wherever the gradient is steep the bands are narrower than one
cell, the lattice undersamples them, and the parity of neighbouring
cells becomes effectively uncorrelated. That is the "noise texture" in
the early screenshots (18:27, 18:31:26).

Where `u` is locally quadratic (near the peak, and in the tails of the
spreading square), `floor(10u) mod 2` is a sampled **Fresnel zone
plate**. An undersampled zone plate is a classical moiré object: it shows
replicated ring centres wherever the local fringe frequency matches a
multiple of the sampling frequency. Those replicas are exactly the
square-symmetric lattice of concentric "islands" in the 18:28–18:29
screenshots, and in `out/jacobi_N64_k02000.png` and
`out/jacobi_res_N256_k02400.png`. They are not features of the
temperature field, which has a single maximum; they are features of the
sampling. Where `u` is locally linear, the pattern is a two-colour
"digital plane" (Sturmian-type pattern), which is what the straight
striped regions are.

Two numerical tests pin this down (`results.txt`):

- **Fixed quantum 0.1, increasing resolution** (64, 128, 256 cells at the
  same physical time): perimeter per cell and 2×2 block entropy both
  *decrease* (0.879 → 0.760 edges per cell; 3.95 → 3.72 bits). Higher
  resolution makes the pattern *cleaner*, converging to the ordinary
  contour map of the smooth solution.
- **Quantum scaled with the cell** (0.1·64/N): the same measures stay
  constant (0.879, 0.894, 0.893). The texture is a function of one
  dimensionless ratio, band width / cell width, i.e. it is pure
  undersampling and is scale-free.
- **Late time, fixed quantum** (t = 0.4, when the 64×64 run shows the
  island lattice of the screenshots):

  | N | isotherm bands | perimeter / N² | 2×2 entropy | components |
  |--:|--:|--:|--:|--:|
  |  64 | 34 | 1.082 | 3.93 | 905 |
  | 128 | 36 | 0.564 | 3.45 |  88 |
  | 256 | 37 | 0.287 | 2.61 |  92 |

  The perimeter per cell halves with each doubling, which is what a
  fixed set of smooth curves does on a finer grid, and the ~900 islands
  collapse to ~90 regions, which is the band count times the few real
  rings. `out/jacobi_late_N256_k32000.png` is simply a contour map.

So the expectation "at higher resolution we would get eddies" is the
wrong way round. Eddies need vorticity, which needs a velocity field and
advection. The heat equation has neither; for `t > 0` its solution is
analytic and its isotherms are smooth nested curves. The only way to
keep the texture at higher resolution is to shrink the quantum with the
cell, which keeps the undersampling.

## 3. The complexity-over-time guess

The guess was: complexity rises after the initial state and later dies.
That is what the pictures show and the measurements confirm (64×64,
clean Jacobi scheme, Neumann on all walls):

| step | max u | isotherm bands | perimeter | 2×2 entropy (bits) | components |
|-----:|------:|---------------:|----------:|-------------------:|-----------:|
|    0 | 100.0 | 1000 |    0 | 0.00 |   1 |
|   20 | 100.0 |  999 | 1136 | 2.21 |  87 |
|  150 |  83.2 |  832 | 3600 | 3.95 | 534 |
|  400 |  49.7 |  490 | 4032 | 3.99 | 646 |
| 1000 |  24.7 |  177 | 4232 | 3.93 | 878 |
| 2000 |  16.3 |   34 | 4432 | 3.93 | 905 |
| 3000 |  14.9 |    7 |  896 | 2.33 |  20 |

But look at the third column. The number of isotherms,
`floor(10·max u) − floor(10·min u)`, **decreases monotonically from the
first step**, as the maximum principle guarantees. The pattern's
perimeter, entropy and component count rise only because bands that
were initially compressed into the one-cell-wide edge of the hot square
spread out and become resolvable. Nothing is created; the already-present
isotherms are being *revealed* as the gradient relaxes, while their
number shrinks. The peak of apparent complexity is when the steepest
band width crosses the cell width, `max|∇u| · dx ≈ 0.1`. After that the
bands widen past the cell, the texture resolves into a few rings, and
when the range of `u` falls below 0.1 the whole picture is one colour.

The entropy of the 2×2 blocks reaches 3.99 bits out of a maximum of 4,
which is why the mid-stage frames look random: at that stage the
sampling is maximally aliased, not the dynamics.

The same curve appears with the script's actual (buggy) dynamics, see
`results.txt`, only faster, because three walls leak heat.

With pure Neumann walls `u → 14.06` (the mean), so the end state is all
black (`floor(140.6)` is even), after a final phase in which the
nearly-flat field crosses a threshold and the whole picture flips. With
the script's walls `u → 0`, also all black.

## 4. What the script actually computes (bugs, for the record)

- After the first key press `u = v` makes the two names refer to the
  same array, so later steps are in-place (Gauss–Seidel-style) sweeps.
  Still stable, but order-dependent and directionally biased.
- `leftWall` uses the global `i` (left at 43 by the initialisation loop)
  instead of its parameter, so the left wall is set to `u[43] + …`.
- Only the interior and the left wall are updated. The other three walls
  stay at 0 forever, which is a Dirichlet-zero condition: heat leaks out.
- The "Total Heat" check sums all cells with equal weight, which is not
  the conserved quantity under ghost-point Neumann rows (see the ternary
  cube note, Corollary 4.5).

None of this changes the conclusions above.

## 5. Literature: has anyone done this?

Nobody has published "parity of the quantized heat equation as a CA",
and for the reasons in §1–2 nobody would present it as a CA. The
relevant prior work falls into four groups.

**Continuous-valued automata and linear automata.**
- S. Wolfram, *A New Kind of Science* (2002), ch. 4 §8 "Continuous
  cellular automata", pp. 155–160 and 243–244: real-valued cells averaged
  with neighbours, a constant added, and then `FractionalPart` taken.
  That is "the heat equation plus mod 1 inside the dynamics", and it
  genuinely generates class-3 and class-4 behaviour. Here the mod is only
  in the display, which is the essential difference. If you want a real
  automaton out of your code, this is the one-line change that gives it:
  apply the mod to the state, not to the colour.
- X.-S. Yang, Y. Young, "Cellular automata, PDEs, and pattern
  formation", arXiv:1003.1983: finite-state automata as discretized
  counterparts of the heat equation, and the relation between the two
  when the number of states is comparable.
- O. Martin, A. Odlyzko, S. Wolfram, "Algebraic properties of cellular
  automata", Comm. Math. Phys. 93 (1984): the heat stencil taken mod p is
  a linear CA (Rule 90 = Pascal's triangle mod 2 = Sierpinski). Again the
  mod is inside the dynamics.
- K. Kaneko, coupled map lattices (1980s): the general framework in
  which "real-valued cells with a local update" lives.
- B. W.-C. Chan, "Lenia and expanded universe", arXiv:2005.03742 (and
  "Lenia: biology of artificial life", 2019): what "continuous cellular
  automata" means on arXiv today (14 hits, all in this family). Lenia is
  a convolution followed by a nonlinear growth map; again the
  nonlinearity is inside the dynamics. A search for "cellular automat"
  together with "moiré" returns nothing, which is consistent with the
  pictures here never having been described as an automaton.

**Diffusion and fluids as cellular automata.** This is where eddies
really come from in the CA world.
- U. Frisch, B. Hasslacher, Y. Pomeau, "Lattice-gas automata for the
  Navier–Stokes equation", Phys. Rev. Lett. 56 (1986): the FHP automaton,
  a Boolean CA whose coarse-graining obeys Navier–Stokes and shows
  vortices.
- B. Chopard, M. Droz, *Cellular Automata Modeling of Physical Systems*
  (Cambridge, 1998), ch. on diffusion; B. Chopard, M. Droz, J. Stat.
  Phys. 64 (1991), diffusion automata.
- T. Toffoli, N. Margolus, *Cellular Automata Machines* (MIT, 1987),
  diffusion rules and the Margolus neighbourhood.

**Thresholded Laplacians that do produce intricate patterns.**
- P. Bak, C. Tang, K. Wiesenfeld, Phys. Rev. Lett. 59 (1987), the
  Abelian sandpile: integer field, discrete Laplacian, toppling threshold.
- W. Pegden, C. Smart, "Convergence of the Abelian sandpile", Duke Math.
  J. 162 (2013); L. Levine, W. Pegden, C. Smart, "Apollonian structure in
  the Abelian sandpile", GAFA 26 (2016). The sandpile pictures are the
  rigorous version of "beautiful patterns from a Laplacian", and the
  nonlinearity (the threshold inside the dynamics) is what makes them
  possible.

**Moiré, aliasing, digital geometry.** This is the correct explanation
of the pictures.
- I. Amidror, *The Theory of the Moiré Phenomenon*, 2nd ed. (Springer,
  2009), vol. I, in particular the chapters on moiré between periodic and
  aperiodic (zone-plate-like) layers and on sampling moiré.
- Undersampled Fresnel zone plates / digital Fresnel holograms:
  arXiv:1911.07997 analyses the replication of concentric zones caused
  by undersampling; the companion paper gives the well-sampling
  condition for the quadratic phase chirp.
- J.-P. Reveillès, *Géométrie discrète, calcul en nombres entiers et
  algorithmique* (thesis, Strasbourg, 1991) and the digital-plane
  literature: `floor(ax + by + c) mod m` patterns are digital planes;
  their 1-D sections are Sturmian words.
- For the monotone quantity: the maximum principle (any PDE text, e.g.
  L. C. Evans, *Partial Differential Equations*, §2.3). Note that in two
  dimensions the number of local maxima of a heat-equation solution is
  *not* monotone in general, unlike the 1-D zero-counting results of
  Sturm/Matano/Angenent, so a claim should be made only about the band
  count, not about the topology of level sets. Survey: R. Magnanini,
  "An introduction to the study of critical points of solutions of
  elliptic and parabolic equations", arXiv:1604.00530.

## 6. Verdict

- **New type of cellular automaton: no.** The pictures are the contour
  map of a diffusing field, quantized and undersampled. The dynamics
  underneath is the simplest linear smoothing system there is.
- **New result in the heat equation: no.** The rise-then-fall of apparent
  complexity is a sampling effect; the physical content (band count) is
  monotone, and the eventual death is just decay to the mean.
- **Eddies at higher resolution: no.** Higher resolution removes the
  texture. The texture survives only if the quantum shrinks with the
  cell.
- **Is there anything publishable?** Yes, if reframed honestly, and not
  on arXiv as a CA paper:
  1. A math-art / recreational piece, "moiré of the heat equation": the
     zone-plate aliasing of quantized isotherms, with the peak-complexity
     time derived from `max|∇u|·dx = 1/q` (for a Gaussian this is a
     closed-form expression in t). Natural venues: the Bridges
     conference, *Journal of Mathematics and the Arts*, *Mathematical
     Intelligencer*.
  2. A short pedagogical note on "complexity measures of an observable
     are not complexity of the system": block entropy and perimeter of a
     thresholded field peak while every norm of the field decays. This
     is a clean cautionary example for a numerical-methods or
     information-theory course.
  3. If what you actually want is a CA with eddies grown out of your
     heat-equation code: implement the FHP lattice gas, or go to a
     reaction–diffusion system (Gray–Scott, Turing patterns), where the
     structure is real and survives refinement.

Any of these can be checked against the images in `out/`: the frames
`jacobi_N64_k*.png` reproduce the sequence of your screenshots, and the
`jacobi_res_N*` and `jacobi_qscaled_N*` frames are the two resolution
tests.
