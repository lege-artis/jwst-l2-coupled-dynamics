# Companion paper — placeholder

This directory will hold the source of the forthcoming arXiv preprint
companion to `lege-artis/jwst-l2-coupled-dynamics`.

## Status

**Placeholder.** The paper has not yet been drafted. This directory is
provisioned at v0.1.0 to make the arXiv-companion-repo structure visible
from the first public commit. Once drafting begins, the layout under
this directory will follow standard arXiv preprint conventions:

```
paper/
├── manuscript.tex             — preprint source (LaTeX)
├── refs.bib                   — references (extends ../shared/reference-bibliography/refs.bib)
├── figures/                   — paper-specific figure renders (regenerable from ../experiments/)
├── build.sh                   — pdflatex + bibtex runner
└── README.md                  — this file (will be updated when manuscript is in flight)
```

## Anticipated scope

The companion preprint is anticipated to cover the dynamical simulation
of probe + JWST-like model under mutual gravitational coupling: the
canonical-tier derivation chain in `../docs/canonical/en/Ch 02..06`,
worked through to numerical results from the Python and Fortran reference
backends, with the closed-form analytical limits (Kepler decoupling, free
symmetric top, gravity-gradient libration, Dzhanibekov asymmetric-top) as
validation anchors.

Specific manuscript scope, target journal vs arXiv-only, and authorship
roster TBD.

## Reproducibility commitment

When the paper lands, the figures-as-shown-in-paper will be exactly
regenerable from this repository via a `paper/build.sh` script and an
`experiments/` directory pinned to a Python environment specification.
Both are deferred until manuscript scope is firm.

## Citation

Until the paper is on arXiv, cite the software directly via
`CITATION.cff` at repository root.
