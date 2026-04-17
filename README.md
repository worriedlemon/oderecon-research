# ODERECON-RESEARCH
This repository contains experiments of reconstructing dynamical systems using varios methods. 

## Overview
> [!NOTE]
> See original ODERECON repository at: https://github.com/worriedlemon/oderecon

This repository contains experiments of:

- Orthogonal polynomials algorithm research experiments
- IO-SINDy algorithm research experiments
- Fractional order reconstruction experiment

## Installation
This repo is designed to be used with `git`, so the installation is pretty simple. There is an alternative scenario below, if you don't want to use `git`.
1. First, clone **the repo** to arbitrary place on your PC: 
```
git clone https://github.com/worriedlemon/oderecon-research.git
```
2. Second, this repository use **a submodule ODERECON**, which should be installed as well. Switch to directory with cloned repository and **write command**:
```
git submodule update --init oderecon
```
3. Run **setup command**:
```matlab
>> setup_research
```

---

*Alternative steps:*
1. Download a zip archive of **this repo** and unpack it in the arbitrary directory on your PC (e.g. `oderecon-research`).
2. Download a zip archive of **submodule repository** (see **Overview** section with link in *Note* block), unpack it inside the location of this repository (i.e. `oderecon-research/oderecon`, probably, the directory `oderecon` will be replaced).
3. Run **setup command** from instruction above.
