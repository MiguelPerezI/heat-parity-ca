# Do the walls select? An experiment (2026-10-22)

Hypothesis (from the functional-information discussion): by the reflection
principle a structure hitting a wall collides with its own mirror image, so a
box should act as a filter, removing structure types whose self-collision
annihilates and keeping those whose self-collision is elastic. If the rule
has several structure types with different self-collision outcomes, the box
population should drift toward the self-elastic types. That would be
differential persistence imposed by the environment, i.e. selection in the
sense of Wong et al.

Rule: 1-D, `u' = frac((u_{i-1}+u_i+u_{i+1})/3 + 0.398)`, nine-decimal
rounding. Scripts: `code/tracker.py`, `catalogue.py`, `catalogue2.py`,
`population.py`, `sparse_box.py`, `slowclass.py`, `slowclass2.py`,
`perstep.py`.

## 1. The structure zoo is almost a single type

Tracking active clusters on a 600-cell torus (steps 1000–4000, tracks of at
least 50 steps, 371 tracks):

| velocity (cells/step) | tracks | mean width |
|---|--:|--:|
| −0.15 … −0.05 | 204 | 3.5–3.7 |
| −0.05 … +0.05 |  19 | 7.2–7.6 |
| +0.05 … +0.15 | 148 | 3.5 |

So there is one dominant mover, speed 0.1 cells per step in either direction,
width about 3.5, and a rare slow/wide class (about 5 % of tracks, width about
7.5, mean lifetime 60 steps versus 105 for the movers), which looks like a
bound pair born in collisions. Structure count on the torus is steady
(43–48 on 600 cells).

## 2. The dominant type reflects elastically

A single random patch on a uniform background in a 200-cell mirror box emits
structures that travel to the walls and come back as V shapes
(`figures/sparse_mirror.png`); the cluster count stabilizes at 14, against 15
for the same patch on a torus (`figures/sparse_periodic.png`). The walls do
not remove movers.

## 3. Population size: torus and box agree

Three seeds, 400 cells, 20 000 steps, cluster counts at
t = 500, 1000, 2000, 4000, 8000, 12000, 16000, 20000:

- torus: 30.0, 26.3, 29.7, 31.3, 30.0, 30.0, 32.7, 28.7
- box:   30.7, 25.0, 26.0, 28.0, 30.3, 29.3, 30.3, 25.3

No decay in the box, no trend. Within noise the box carries the same number
of structures as the torus; if anything a few per cent fewer.

## 4. Population composition: a marginal difference, and a trap

Tracks with lifetime ≥ 50 steps whose mean position lies in the bulk (cells
50–350, same criterion in both geometries), five seeds, steps 3000–10 000:

| geometry | bulk tracks | slow/wide fraction | mean life slow | mean life fast |
|---|--:|--:|--:|--:|
| torus | 1826 | 0.061 (111) | 59 | 104 |
| box   | 1966 | 0.043 (85)  | 60 | 105 |

The slow fraction is lower in the box by a factor 1.4, about 1.8 standard
deviations: suggestive, not established. The lifetimes of both classes are
identical in the two geometries, which argues against the walls killing
either class; if the composition really differs, the likelier cause is a
lower birth rate of bound pairs, since pairs are born in collisions and the
box holds slightly fewer movers.

The trap: a fitted velocity classifies a reflection (a V-shaped track) as
"slow", and the periodic seam of the tracker splits structures at cell 0
into stubs that also look slow. Both artifacts pile up within five cells of
the ends (8.5 "slow" tracks per cell at the box walls, 4.2 at the torus
seam, against 0.3 in the bulk). Any claim that structures are born or die at
the walls must use a per-step measure, not fitted velocities. Section 5 does
that.

## 5. Per-step densities (no tracking)

Ten seeds, 400 cells, steps 3000–13 000, every step, bulk = cells 50–350.
A "wide" cluster has width ≥ 6 (this catches bound pairs and crossing blobs
alike; the point is that the same measure is applied to both geometries).

| geometry | clusters / cell / step | wide / cell / step | wide fraction |
|---|--:|--:|--:|
| torus | 0.0723 | 0.02937 | 0.4064 |
| box   | 0.0723 | 0.02948 | 0.4077 |

Identical to three significant figures. Within ten cells of the walls the
box has 0.0762 clusters per cell per step (5 % above the bulk, the
reflections pile up there) with a wide fraction of 0.36; the torus "ends"
show 0.094 and 0.30, which is the seam artifact and not comparable.

## 6. Verdict

The walls do not select in this rule. Three independent measurements agree:

1. The dominant structure (speed 0.1, width 3.5) reflects elastically; the
   wall collision is the same event as the bulk head-on collision of two
   opposite movers, which is also elastic (V shapes in both pictures).
2. Structure counts and per-step densities in the bulk of the box are
   indistinguishable from the torus (0.0723 vs 0.0723 clusters per cell per
   step; wide fraction 0.408 vs 0.406), and the mean lifetimes of both
   structure classes are the same in both geometries (60/59 slow, 105/104
   fast).
3. The marginal difference in tracked composition (slow fraction 0.043 vs
   0.061, 1.8 σ) is not reproduced by the per-step measure and should be
   regarded as noise or a tracking artifact.

Why the hypothesis fails here: selection by self-collision needs at least
two structure types with different self-collision outcomes. This rule has
essentially one mover whose self-collision is elastic, and a rare bound pair
whose lifetime is set by bulk collisions, not by the walls. The box is a
mirror and nothing more, which is what the article already states.

What would make the hypothesis testable: a rule with a genuine zoo of
gliders of different speeds and widths, where the outcome of a symmetric
head-on collision differs between types. Candidates are the one-dimensional
class-4 rules with several glider species (elementary rule 110 is the
classical example, with a dozen known gliders and tabulated collision
outcomes, though it is discrete rather than continuous), or Lenia-type
rules in two dimensions, where the symmetric four-body collision at a
corner is a second, independent filter. The experiment designed here
(sparse launch, population counts, per-step composition, torus vs box,
same seeds) transfers to those rules unchanged.
