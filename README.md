# PhenoCellPy

See more details at https://www.biorxiv.org/content/10.1101/2023.04.12.535625v2

PhenoCellPy is an open-source python package that defines methods for modeling changes of cell behaviors. 
It also has pre-defined sequences of behaviors.

Modeling of biological systems requires a degree of model complexity that can be daunting to new and experienced 
modelers alike. Defining cellular behaviors or states and how one cellular state leads to the next is not trivial. To 
make the development of new biological models we have created PhenoCellPy. PhenoCellPy defines methods and Python 
classes that make the definition of sequence of cell behaviors easier. The main class is the Phenotype, which is the 
container of cell behaviors. Phenotype here can mean the cell cycle, the stages of necrosis, the fact that the cell is 
alive, the fact that it is dead, _etc_. The Phenotype is made of one or more Phases, each Phase defines the target 
volume for the cell and the volume change rates it displays. It also defines which is the next Phase of the Phenotype, 
what conditions trigger Phase change, if cell division occurs when exiting the Phase, what behaviors occur immediately 
on Phase entry and just before Phase exit (_e.g._, changing the target volume). The cell volume dynamics are 
handled by the Cell Volume class. The cell volume is subdivided among the solid and fluid cytoplasm, solid and fluid 
nucleus, and a calcified fraction.

By default, the Phase transition can be either stochastic (with a set Phase transition rate) or deterministic (with a 
set Phase period). The user can also define a custom transition function. One of our pre-built Phenotypes, the Necrotic 
Standard Model, uses a custom transition function for its first Phase (which represents the osmotic swelling of a 
necrotic cell). The custom transition function monitors the cell volume and changes from the swelling Phase to the 
ruptured cell Phase when the cell reaches its rupturing volume. The Phase class can, optionally, check if a cell should 
exit the Phenotype and enter senescence.

As mentioned, a Phenotype can be any sequence of Phases. For instance, the cell cycle is a Phenotype. The modeler using 
PhenoCellPy must call the division methods of the modeling framework being used to divide the cell. 

PhenoCellPy is intended to be used with other python-based modeling framework, _e.g._, CompuCell3D, Tissue Forge, as 
an embedded model. PhenoCellPy is inspired by the phenotype definitions of PhysiCell [1]. 

[1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, PhysiCell: an Open Source Physics-Based 
Cell Simulator for Multicellular Systems, PLoS Comput. Biol. 14(2): e1005991, 2018. DOI: 10.1371/journal.pcbi.1005991

## Requirements

PhenoCellPy only requirements are

* Python 3 support
* NumPy
* SciPy

## Installation

PhenoCellPy's _alpha_ version does not have any installer. To use it you should clone or download its GitHub 
repository and add its folder to the simulation's system path. E.g., 

```
import sys
sys.path.extend(['C:\\PhenoCellPy_Dowload_Folder', 
                'C:/PhenoCellPy_Dowload_Folder'])
```
