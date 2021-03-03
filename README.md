# nademo
Nucleid acid demo using Seamless

**Installation**

Running the initial-analysis notebook requires numpy, scipy, pandas and matplotlib installed

It also requires nglview to be installed and configured correctly with Jupyter

It may be easiest to run Jupyter inside the Seamless docker container, since nglview and all dependencies work correctly.

Finally, it requires x3dna-dssr to be available in the system path.
See https://x3dna.org/ and click License to obtain the program.
To install it in a running Seamless Docker container:
- Obtain the Linux binary and copy it into nademo
- Run `seamless-shell-existing`. This will connect into the running Docker container as root. You can then move it to /usr/bin
