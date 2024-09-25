# LinearPowerFlow.jl
A simple, easy to use fixed-point linear power flow model based on [[1]](@ref).

!!! NOTE this package is in the process of deprecation.
The model will be a part of [BusInjectionModel.jl](https://github.com/NLaws/BusInjectionModel.jl), where usage will look like:
```julia
m = JuMP.Model(My.Optimizer)
net = CommonOPF.Network(<from-yaml-or-opendss-file-path>)
BranchFlowModel.build_bim!(m, net, BranchFlowModel.FixedPointLinear)
JuMP.optimize!(m)
# a standard results-getter is a work in progress, for now have to use JuMP.value
```
where none of the module names are necessary, just added for illustration.


### [1]
Bernstein, Andrey, and Emiliano Dall'Anese. "Linear power-flow models in multiphase distribution networks." 2017 IEEE PES Innovative Smart Grid Technologies Conference Europe (ISGT-Europe). IEEE, 2017.
