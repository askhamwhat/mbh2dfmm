
# MBHFMM2D

A Fortran library for the rapid and stable evaluation of
modified biharmonic potentials. Both high order charges
(dipoles, quadrupoles, octopoles) and high order derivatives
(up to order 3) are allowed, which facilitates Stokeslet
and stresslet-like formulations.

NOTE: feature requests, bug reports, and pull requests
will be largely ignored for this repository (with the exception
of adding proper compile flags to the Makefile and
other compatibility issues). This
repository is likely to move to a more permanent location
in the future.

## Compile instructions

You can use the configure_makefile.sh bash script to
generate a Makefile. To do so with standard options use

```
./configure_makefile.sh
```

For options use

```
./configure_makefile.sh --help
```

If this doesn't work for you, a couple of makefiles
with some standard options are included in the mkfiles
folder.

## Using the code

To run the example code, use:

```
make mbhfmmtest
```

See the driver for that routine in the test folder
for example usage.

