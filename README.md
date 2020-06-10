# Phase diagrams from classical DFT computations

This is a repository containing the source code and data for my master thesis in the *Service de Physique des Systèmes Complexes et Mécanique Statistique* at the *Université Libre de Bruxelles*. 

The code is built upon the classicalDFT library, written by my supervisor James F. Lutsko. The original repository for the library is

> https://github.com/jimlutsko/classicalDFT <br>
> commit 9b5e7e6e5f0140b4cbc87a15de89753ad89e5144 <br>

![alt text](Cedric_Memoire_Other/Illustration/diagram_WHDF_rc12.png?raw=true "Typical phase diagram for colloid-like interactions")

## Installation

1. Setup and compile the library:

> cd Lib <br>
> ../Config.sh <br>
> ../dft\_make.sh <br>

2. Compile the source code for phase diagram calculations:

> cd Cedric\_Memoire\_Code <br>
> ./compile\_everything.sh <br>

3. Other pieces of code can be found elsewhere and are generally come along scripts with explicit names indicating their purpose, e.g. 

> compile.sh <br>
> collect\_data.sh <br>
> run\_coexistence.sh <br>
> plot\_everything.sh <br>

## Directory structure

| Directory | Purpose |
| --------- | ------- |
| **Cedric\_Memoire\_Code**  | Source code for my DFT calculations |
| **Cedric\_Memoire\_Data**  | Where the code is executed and data collected |
| **Cedric\_Memoire\_Other** | Additional content not related to DFT calculations |
| **Lib** | classicalDFT library |
| **SS31-Mar-2016** | Pre-compiled weights for spherical integrations |

## Branches 

The main branch is *Cedric_Memoire*. <br>
The branches *Cedric_Memoire_preFeb2020* and  *Cedric_Memoire_preMay2020* are older versions of the code that I kept as backup. 

