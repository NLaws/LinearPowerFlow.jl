function build_lpf(m::JuMP.AbstractModel, p::Inputs)
    add_variables(m, p)
    constrain_voltage(m, p)
    constrain_substation_power(m, p)
    constrain_bounds(m, p)
end


function add_variables(m::JuMP.AbstractModel, p::Inputs)
    N = p.Nnodes
    T = p.Ntimesteps

    @variable(m, P₀[1:T])
    @variable(m, Q₀[1:T])
    # TODO Pj Qj variables and change p.P/Qref to P/Qload dicts ?

    @variable(m, v[1:N, 1:T] >=0)  # vlo and vhi applied in constraints s.t. they end up in C matrix
end


function constrain_voltage(m::JuMP.AbstractModel, p::Inputs)
    v = m[:v]

    # Bernstein & Dall'anese 2017 equ. 5b, 9b, 9c
    for t = 1:p.Ntimesteps
        # TODO should we update the M matrices with each time step? Depends on if Vref is updated with time
        invDiagConjVref = inv(Diagonal(vec(conj(p.vref[:, t]))))  # TODO mv to Inputs so only calculated once (and Mwye)
        Mwye = [p.invYLL * invDiagConjVref  -im * p.invYLL * invDiagConjVref]
        Kwye = inv(Diagonal(abs.(p.vref[:,t]))) * real( Diagonal(conj(p.vref[:,t])) * Mwye )
        b = abs.(p.vref[:,t]) - Kwye * vec([p.Pref[:,t] / p.Sbase; p.Qref[:,t]] / p.Sbase)

        KL = Kwye[:, 1:p.Nnodes]
        KR = Kwye[:, p.Nnodes+1:end]

        @constraint(m, [j = 1:p.Nnodes],
            v[j,t] == sum( KL[j,n] * p.Pref[n,t] / p.Sbase for n in p.Nnodes) 
                    + sum( KR[j,n] * p.Qref[n,t] / p.Sbase for n in p.Nnodes) 
                    + b[j]
        )
    end
    #= 
    equ (10) distributable voltage estimation, less accurate

    v₀ = [complex(p.V0, 0)]
    w = -p.invYLL * p.YL0 * v₀
    W = Diagonal(w)
    big_one = vec(ones(p.Nnodes))

    abs.(W) * (big_one + real(inv(W) * Mwye) * vec([p.Pref[:,t] / p.Sbase; p.Qref[:,t]] / p.Sbase))
    =# 
    p.Nequality_cons += p.Nnodes * p.Ntimesteps
    nothing
end



"""
    constrain_substation_power(m::JuMP.AbstractModel, params::YPQVarraysNodeByTime)

P₀ & Q₀ definitions ∀ t ∈ T
"""
function constrain_substation_power(m::JuMP.AbstractModel, p::Inputs)
    P₀ = m[:P₀]
    Q₀ = m[:Q₀]
    v₀ = [complex(p.V0, 0)]
    # parameters using notation from Bernstein and Dall'anese 2017 equ.s 5c, 13a, 13b
    w = -p.invYLL * p.YL0 * v₀  # can be replaced with ntwk.vnom?
    c = Diagonal(v₀) * (conj(p.Y₀₀) * conj(v₀) + conj(p.Y0L) * conj(w))
    creal = real(c)[1]
    cimag = imag(c)[1]
    
    # s₀ = G*x + c
    for t = 1:p.Ntimesteps
        invDiagConjVref = inv(Diagonal(vec(conj(p.vref[:, t]))))
        Mwye = [p.invYLL * invDiagConjVref  -im * p.invYLL * invDiagConjVref]
        gwye = Diagonal(v₀) * conj(p.Y0L) * conj(Mwye)

        @constraint(m, 
            P₀[t] == 1/p.Sbase * dot(real(gwye), [p.Pref[:,t]; p.Qref[:,t]]) + creal
        )
        @constraint(m, 
            Q₀[t] == 1/p.Sbase * dot(imag(gwye), [p.Pref[:,t]; p.Qref[:,t]]) + cimag
        )
    end
    
    p.Nequality_cons += 2 * p.Ntimesteps
    nothing
end


function constrain_bounds(m::JuMP.AbstractModel, p::Inputs)

    @constraint(m, con_vhi[n=1:p.Nnodes, t=1:p.Ntimesteps],
        m[:v][n, t] ≤ p.vhi
    )
    p.Nlte_cons += p.Nnodes * p.Ntimesteps

    @constraint(m, con_vlo[n=1:p.Nnodes, t=1:p.Ntimesteps],
        p.vlo - m[:v][n, t] ≤ 0
    )
    p.Nlte_cons += p.Nnodes * p.Ntimesteps

    @constraint(m, con_P0hi[t=1:p.Ntimesteps],
        m[:P₀][t] ≤ p.P0hi
    )
    p.Nlte_cons += p.Ntimesteps

    @constraint(m, con_P0lo[t=1:p.Ntimesteps],
        p.P0lo - m[:P₀][t] ≤ 0
    )
    p.Nlte_cons += p.Ntimesteps

    @constraint(m, con_Q0hi[t=1:p.Ntimesteps],
        m[:Q₀][t] ≤ p.Q0hi
    )
    p.Nlte_cons += p.Ntimesteps

    @constraint(m, con_Q0lo[t=1:p.Ntimesteps],
        p.Q0lo - m[:Q₀][t] ≤ 0
    )
    p.Nlte_cons += p.Ntimesteps
    nothing
end


function check_existence_condition(m,p)
    # TODO: still applies? if so update the P,Q vectors
    upperbound = p.V0^2 / (4 * matrixnormstar(p.Z))
    Snorm = zeros(Float64, p.Ntimesteps)
    for t in 1:p.Ntimesteps
        P = value.(m[:P][:,t])
        Q = value.(m[:P][:,t])
        Snorm[t] = LinearAlgebra.norm( sqrt(sum( P[j]^2 + Q[j]^2 for j in 0:p.Nnodes)) )
    end
    return all(Snorm[t] < upperbound for t in 1:p.Ntimesteps)
end


function matrixnormstar(A::AbstractArray{<:Complex{Float64}, 2})
    rownorms = zeros(Float64, size(A, 1))
    for r = 1:size(A,1)
        rownorms[r] += LinearAlgebra.norm(abs.(A[r,:]))
    end
    return maximum(rownorms) 
end
