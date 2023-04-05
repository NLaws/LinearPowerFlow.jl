function get_constraints_by_variable_name(m, v::String)
    ac = ConstraintRef[]
    for tup in list_of_constraint_types(m)
        append!(ac, all_constraints(m, tup[1], tup[2]))
    end
    filter( cr -> occursin(v, string(cr)), ac )
end

function nan2zero!(A::AbstractArray{<:Complex})
    for i in eachindex(A)
        @inbounds A[i] = Complex(
            ifelse(isnan(real(A[i])), 0, real(A[i])),
            ifelse(isnan(imag(A[i])), 0, imag(A[i]))
        )
    end
end