mutable struct Inputs
    Ntimesteps::Int
    Y::AbstractArray{<:Complex, 2}
    Z::AbstractArray{<:Complex, 2}
    Nnodes::Int
    Sbase::Real
    vlo::Real
    vhi::Real
    Pref::AbstractArray{Float64, 2}  # (node, time)
    Qref::AbstractArray{Float64, 2}
    vref::AbstractArray{<:Complex, 2}
    V0::Real
    YLL
    Y0L
    YL0
    invYLL
    Y₀₀
    pf::Float64
    Nequality_cons::Int
    Nlte_cons::Int
    P0lo
    P0hi
    Q0lo
    Q0hi

    function Inputs(;
        Y=missing,  # can provide one or both of Y, Z
        Z=missing,
        Pref,
        Qref,
        vref,
        Sbase=1.0,
        Vbase=1.0,
        V0=1.0,
        vlo=0.95,
        vhi=1.05,
        P0lo=-10,
        P0hi=10,
        Q0lo=-10,
        Q0hi=10
    )

        if !ismissing(Y) & !ismissing(Z)
            # do nothing
        elseif ismissing(Y) & !ismissing(Z)
            Y = inv(Z)
            # TODO probably need inf2ero or nan2zero here
        elseif ismissing(Z) & !ismissing(Y)
            Z = inv(Y) 
            nan2zero!(Z)
        else
            error("Must provide one or both of Y, Z")
        end

        Ibase = Sbase / (Vbase * sqrt(3))
        Zbase = Vbase / (Ibase * sqrt(3))
        Ybase = 1 / Zbase
        Y = Y ./ Ybase
        Z = Z ./ Zbase
        # pu loads
        Nnodes = size(Pref, 1)
        Ntimesteps = size(Pref, 2)
        @assert Nnodes+1 == size(Y, 1)  # b/c Y includes substation

        YLL = Y[2:end, 2:end]
        Y0L = Y[1:1, 2:end]
        YL0 = Y[2:end, 1:1]
        invYLL = inv(YLL)
        Y₀₀ = Y[1:1,1:1]

        Nequality_cons = 0 
        Nlte_cons = 0

        return new(
            Ntimesteps,
            Y,
            Z,
            Nnodes,
            Sbase,
            vlo,
            vhi,
            Pref,
            Qref,
            vref,
            V0,
            YLL,
            Y0L,
            YL0,
            invYLL,
            Y₀₀,
            0.1,
            Nequality_cons,
            Nlte_cons,
            P0lo,
            P0hi,
            Q0lo,
            Q0hi
        )
    end
end


"""

parse all inputs from an OpenDSS model, including:
- nodal admittance matrix
- Sbase
- Vbase
- Pload
- Qload
"""
function Inputs(dss_path::String)
    dss("Redirect $dss_path")
    dss("Solve")  # is calcv enough to get Y ?

    if !(OpenDSSDirect.Solution.Converged() == true)
        @warn "OpenDSS solution is not converged for $dss_path"
    end

end