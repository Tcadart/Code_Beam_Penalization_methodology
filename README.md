# Lattice Structures Study Tool
This repository contains the simulation framework and code for finding the parameters of the penalization beam methodology as described in "A beam element penalization methodology for additively manufactured lattice structures".

## Overview
This Python project is dedicated to the study of lattice structures based on density and simulation methods. Utilizing the finite element analysis capabilities of Abaqus, this code addresses the challenges associated with the simulation of Timoshenko beams and introduces a penalty methodology to facilitate the use of beam elements in the study of lattice structures. This approach allows for the exploration of lattice structures across a broad spectrum of densities and geometries.

### Key Features
- Finite element simulation of Timoshenko beams in Abaqus.
- Penalty methodology for the study of lattice structures.
- Homogenization process to compare lattice structures.

## Prerequisites
- [Abaqus](https://www.3ds.com/products-services/simulia/products/abaqus/) for finite element analysis.
- [Micromechanics plugin for Abaqus](https://www.linkedin.com/pulse/micromechanics-plugin-abaqus-ross-mclendon/) to simplify the implementation of homogenization.
- Python (tested on version 3.11)
- Python packages:
  - NumPy
  - SciPy
  - Matplotlib

## File Description
- `Parameters_optimization.py` : Main file for generating data for multiple lattice geometries, densities, behavior
- Utilities code :
  - `Material.py` : Define material properties for lattice
  - `Lattice_description.py` : Define lattice geometries
- Graph generation code :
  - `Graph_LZone.py`
  - `Graph_penalization_coefficient.py`
  - `Graph_frequency.py`
  - `Graph_plasticity.py`
