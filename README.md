# Python Finite Element Beam Analysis

An interactive structural analysis application developed using Python and Streamlit. The tool uses the finite element method (FEM) to calculate the response of simply supported Euler–Bernoulli beams subjected to various loading conditions.

## Features

- Interactive browser-based interface built with Streamlit.
- Adjustable beam length, Young's modulus and second moment of area.
- Multiple point loads at user-defined positions.
- Uniform, triangular and trapezoidal distributed loads, including overlapping loads.
- Automatic calculation of support reactions, shear forces, bending moments and beam deflections.
- Live beam loading diagram with support symbols and applied loads.
- Graphical results highlighting maximum bending moment and maximum downward deflection.
- Force and moment equilibrium checks.

## Engineering Methodology

The application uses the Euler–Bernoulli beam theory and a finite element formulation with two degrees of freedom per node: vertical displacement and rotation.

The solver assembles a global stiffness matrix using cubic Hermite beam elements. Distributed loads are converted into equivalent nodal forces using three-point Gaussian quadrature.

Simply supported boundary conditions are applied by constraining vertical displacements at both ends while allowing rotation.

The resulting system of equations is solved using NumPy.

## Validation

The program has been checked against analytical beam theory for several loading scenarios:

- Centrally applied point load.
- Off-center point load.
- Full-span uniformly distributed load.
- Full-span triangular distributed load.
- Partial-span uniformly distributed load.
- Partial-span trapezoidal distributed load.

Support reactions were also checked using force and moment equilibrium. Combined and overlapping load cases were tested.

## Technologies

- Python
- NumPy
- Matplotlib
- Streamlit

## Running the Application

1. Download or clone this repository.
2. Install the dependencies:

   `pip install -r requirements.txt`

3. Run the Streamlit application:

   `python -m streamlit run Beam_app.py`

4. Open the local address displayed in the terminal.

## Assumptions and Limitations

- Simply supported beams only.
- Constant flexural rigidity (EI) along the beam.
- Linear elastic material behavior.
- Small deflections.
- Euler–Bernoulli beam theory, neglecting shear deformation.
- Static, downward loading only.
- Results are intended for educational and portfolio demonstration purposes, not professional structural design verification.

## Project Purpose

This project was developed to apply structural mechanics and finite element theory through programming, combining numerical analysis with an interactive engineering interface.
