# Apparent complexity of quantized diffusion: contour maps, aliasing and moiré

Code, figures and LaTeX source for the article of the same title by
Miguel E. Pérez Ibarra and Jaziel Flores Rodríguez. The compiled article
is [`main.pdf`](main.pdf); the source is [`main.tex`](main.tex).

The article explains a set of black-and-white frames produced by a
heat-equation simulation, which look like a cellular automaton with a rich
transient. They are not one.

## The frames are aliased contour maps

Colouring the cells of a finite-difference heat simulation by the parity of
`floor(10·u)` gives frames whose apparent complexity rises from nothing,
peaks in a lattice of symmetric concentric islands, and dies:

| step 150 | step 1000 | step 2000 |
|:--:|:--:|:--:|
| ![](figures/parity_k150.png) | ![](figures/parity_k1000.png) | ![](figures/parity_k2000.png) |

They are the contour map of the temperature with isotherms every 0.1 degrees
and alternating fill. The rise is aliasing and the islands are the
replicated centres of an undersampled Fresnel zone plate. The same state on a
four times finer grid, with the same quantum, is an ordinary contour map:

| 64 × 64 cells, t = 0.4 | 256 × 256 cells, same state |
|:--:|:--:|
| ![](figures/parity_k2000.png) | ![](figures/parity_late_N256.png) |

At an earlier time on a finer grid the effect is at its purest: a central zone
plate surrounded by a lattice of replica centres, each a small zone plate of
its own, with the stripes of digital planes between them. The field has one
maximum and its isotherms never cross; everything else is the sampling.

<p align="center"><img src="figures/parity_res_N256_t003.png" width="520"></p>

Three elementary results make this exact (Theorems 2.1, 3.1, 4.1 of the
article):

- the band count and the total variation are nonincreasing under Neumann
  heat flow on a box;
- the phase-averaged perimeter of the sampled pattern is exactly
  `sum over edges of dist(q·|jump across the edge|, 2Z)`, a triangle wave of the
  edge jump, and the discrete total variation is nonincreasing, so every
  increase of the perimeter comes from edges whose jump exceeds one quantum;
- the block entropy of the pattern has two non-commuting limits: `0` as the
  grid is refined, maximal as the quantum is refined.

Scripts: [`code/analyze.py`](code/analyze.py) (complexity measures against
time, resolution tests), [`code/phase_theorem.py`](code/phase_theorem.py)
(the random-phase identity checked against Monte-Carlo), results in
[`code/results_partI.txt`](code/results_partI.txt).

## Running the code

Requirements: Python 3.12 with NumPy. No PIL is
needed, PNGs are written by [`code/pngio.py`](code/pngio.py), and
[`code/png2pdf.py`](code/png2pdf.py) converts them into vector PDFs (one
rectangle per run of equal cells) for the article, so the figures stay sharp
at any zoom.

```
cd code
python3.12 analyze.py                 # tables and frames (writes ../out/)
python3.12 phase_theorem.py           # random-phase identity, discrete TVD
python3.12 restest.py                 # resolution tests
```

`ASSESSMENT.md` is the original analysis of the frames and `PROPOSAL.md` the
theorem proposal that became the article.

## Building the article

```
pdflatex main.tex && pdflatex main.tex
```

## Related

- The reflecting walls of the simulation use the ternary-cube rule from the
  note *One case, not twenty-seven: the ternary cube as the index set for
  boundary stencils of finite-difference schemes in any dimension* by the same
  authors.
- Putting the modulus into the dynamics instead of the display gives
  Wolfram's continuous cellular automata, which do have localized structures.
  How those structures, and the gliders of rules 54 and 110, meet reflecting
  walls is the subject of the companion paper *Reflecting walls are
  self-collisions: cellular automata in a mirror box*, code at
  [MiguelPerezI/mirror-walls](https://github.com/MiguelPerezI/mirror-walls).

## License

Code under the MIT License (see `LICENSE`). The article text and figures are
© the authors.
