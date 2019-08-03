# Gaussian process modifier adaptation
The code in this repository is based on the works in [[1]](#1)[[2]](#2) and is written in Matlab. It is a variation of a modifier adapation algorithm using Gaussian processes (GP-MA). To cite GP-MA please use the publications [[1]](#1)[[2]](#2).

## Getting started
Download the entire folder containing all Matlab files and add the entire folder and subfolders to the Matlab path. Next run [RTO_MA_GP](/RTO_MA_GP.m), which should run the pre-defined problem as defined in [[2]](#2). Once this works the "true" plant model can be modified in [Plant_model](Plant_model.m) with the corresponding objective defined in [Plant_data](Plant_data.m). The "approximate" plant model can be modified in [Approx_model](Approx_model.m) with the corresponding approximate objective defined in [Approx_data](Approx_data.m). Note the plotting only works for two dimensional problems. 

## References
[1] E.A. del Rio Chanona, J.E. Alves Graciano, E.Bradford, and B.Chachuat, [Modifier-Adaptation Schemes Employing Gaussian Processes and Trust Regions for Real-Time Optimization](https://www.sciencedirect.com/science/article/pii/S2405896319301211), IFAC-PapersOnLine, vol. 52, no. 1, pp. 52-57, 2019. 
<a name="1">
</a>

[2] T.A. Ferreira, H.A. Shukla, T. Faulwasser, C.N. Jones, and D. Bonvin, [Real-time optimization of uncertain process systems via modifier adaptation and Gaussian process](https://ieeexplore.ieee.org/abstract/document/8550397), In European Control Conference (ECC) 2018, pp. 465-470. 
<a name="2">
</a>

## Legal information
This project is licensed under the MIT license – see LICENSE.md in the repository for details.
