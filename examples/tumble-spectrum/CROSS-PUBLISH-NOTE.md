# Cross-publish notice

This example is also hosted in the sibling project
[lege-artis/fourier](https://github.com/lege-artis/fourier) under
`examples/jwst-l2-tumble-spectrum/`, where it demonstrates DFT applied
to non-signal-processing data.

The example reads NDJSON trajectory output from this project's
`run_first_example.py` and produces an FFT diagnostic of the
free-precession Euler-frequency peak. The mathematical primitive
(DFT) is canonical-source-of-truth at fourier; the trajectory data
is canonical-source-of-truth here.
