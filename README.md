# Random Volcanic Climate [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6799309.svg)](https://doi.org/10.5281/zenodo.6799309)

This repository contains simulation code for a research project on the long-term variability of Earth's carbon cycle and climate that I worked on with [Minmin Fu](https://minminfu.github.io/). The preprint is available [on arXiv](https://arxiv.org/abs/2208.02793) and the manuscript was [published](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2022GC010611) by the *Geochemistry, Geophysics, Geosystems* journal with open access.

-----

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> Random Volcanic Climate

It is authored by Mark Baum.

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.
