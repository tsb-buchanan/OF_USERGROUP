#!/bin/bash
# Tested: R. Dwight - Jan 2020 - OF7.0

### 1. Mesh generation and solver
blockMesh > log.blockMesh
checkMesh > log.checkMesh

simpleFoam > log.simpleFoam 

# Note: `./monitor.sh` at this point, plots regularly updated convergence.

### 2. Postprocessing: Output wall-shear stresses to 'postProcessing' dir
simpleFoam -postProcess -func wallShearStress > log.wallShearStress

### 3. Postprocessing: Sample solution according to system/sampleDict, output
###    in 'postProcessing' dir
postProcess -func sampleDict > log.sampleDict

