
"Compute bus connection indices for power system components based on a continuous mapping"
function _calc_connection_indices(
    data::Dict{Symbol, _DF.DataFrame},
    bus_to_idx::Vector{Pair{Int, Int}},
    label::Symbol
)

    C = Dict{Symbol, Matrix{Int}}()
    for (key, df) in data

        if contains(string(key), "bus")
            continue
        end
        vars = filter(x -> contains(string(x), "bus"), names(df))
        sort!(vars)
        bus_ids = Matrix{Int64}(df[!, vars])
        replace!(bus_ids, bus_to_idx...)
        C[Symbol(key, label)] = bus_ids
    end
    return C
end

"Compute bus and branch admittance matrices"
function _calc_admittance_matrices(
    n_bus::Int,
    edge_data::_DF.DataFrame,
    shunt_data::_DF.DataFrame,
    indices::Dict{Symbol, Matrix{Int}}
)

    data = shunt_data
    # Shunt admittance matrix (n_bus x n_bus)
    ys = zeros(ComplexF64, n_bus)
    if !isempty(data)
        name = first(filter(x -> contains(string(x), "shunt"), keys(indices)))
        ys[indices[name]] = data.gs + data.bs * im
    end
    Ys = spdiagm(n_bus, n_bus, ys)

    mask = Bool.(edge_data.br_status)
    data = edge_data[mask, :]
    # Branch series admittance
    y = pinv.(data.br_r + data.br_x * im)
    # From and to shunt leg terms of branch
    lg_fr = (data.g_fr + data.b_fr * im)
    lg_to = (data.g_to + data.b_to * im)
    # Branch transformer tap
    tr = data.tap .* cos.(data.shift)
    ti = data.tap .* sin.(data.shift)
    t = tr + ti * im
    # Branch admittance diagonal matrices (n_edge x n_edge)
    n_edge = size(data, 1)
    Y_ff = spdiagm(n_edge, n_edge, (y + lg_fr) ./ abs2.(t))
    Y_ft = spdiagm(n_edge, n_edge, -y ./ conj.(t))
    Y_tf = spdiagm(n_edge, n_edge, -(y ./ t))
    Y_tt = spdiagm(n_edge, n_edge, y + lg_to)
    # Get branch connection matrices (n_edge x n_bus)
    name = first(filter(x -> contains(string(x), "branch"), keys(indices)))
    edges = indices[name][mask, :]
    C_f = sparse(1:n_edge, edges[:, 1], 1, n_edge, n_bus)
    C_t = sparse(1:n_edge, edges[:, 2], 1, n_edge, n_bus)

    Yf = sparse(Y_ff * C_f + Y_ft * C_t)
    Yt = sparse(Y_tf * C_f + Y_tt * C_t)
    Ybus = transpose(C_f) * Yf + transpose(C_t) * Yt + Ys

    return (Ybus = Ybus, Yf = Yf, Yt = Yt)
end

function calc_connection_indices(pm::_PM.AbstractPowerModel; label::Symbol = :_index)

    bus = vec(get_pm_value(pm, :bus, ["bus_i"], Array{Any, 2}))
    bus_to_idx = bus .=> 1:length(bus)
    
    data = Dict{Symbol, _DF.DataFrame}()
    for key in keys(filter(p->isa(p.second, Dict{Int, Any}) && !isempty(p.second), _PM.ref(model)))
        vars = get_pm_key(pm, key)
        data[key] = get_pm_value(pm, key, vars, _DF.DataFrame)
    end
    return _calc_connection_indices(data, bus_to_idx, label)
end

function calc_connection_indices(data::Dict{Symbol, _DF.DataFrame}; label::Symbol = :_index)
    bus = data[:bus].bus_i
    bus_to_idx = bus .=> 1:length(bus)
    return _calc_connection_indices(data, bus_to_idx, label)
end

function calc_admittance_matrices(pm::_PM.AbstractPowerModel, indices::Dict{Symbol, Matrix{Int}})
    n_bus = length(_PM.ref(pm, :bus))
    vars = ["b_fr", "b_to", "br_r", "br_x", "g_fr", "g_to", "tap", "shift"]
    mask = findall(Bool.(vec(get_pm_value(pm, :branch, ["br_status"], Array{Any, 2}))))

    shunt_data = get_pm_value(pm, :shunt, ["gs", "bs"], _DF.DataFrame)
    branch_data = get_pm_value(pm, :branch, vars, _DF.DataFrame; mask = mask)
    return _calc_admittance_matrices(n_bus, branch_data, shunt_data, indices)
end

function calc_admittance_matrices(data::Dict{Symbol, _DF.DataFrame}, indices::Dict{Symbol, Matrix{Int}})
    n_bus = size(data[:bus], 1)
    branch_data = data[:branch]
    if haskey(data, :shunt)
        shunt_data = data[:shunt]
    else
        shunt_data = _DF.DataFrame()
    end
    return _calc_admittance_matrices(n_bus, branch_data, shunt_data, indices)
end


"""
    calc_susceptance_matrices(data::Dict{Symbol, _DF.DataFrame}, indices::Dict{Symbol, Matrix{Int}})

Compute `B'` and `B''` matrices for Fast Decoupled Power Flow in BX version
"""
function calc_susceptance_matrices(data::Dict{Symbol, _DF.DataFrame}, indices::Dict{Symbol, Matrix{Int}})

    ref = :branch
    # Compute B'
    data_temp = deepcopy(data)
    if !isempty(get(data_temp, :shunt, _DF.DataFrame()))
        data_temp[:shunt][!, "bs"] .= 0.0
    end
    for (key, value) in zip(["b_fr", "b_to", "tap"], [0.0, 0.0, 1.0])
        data_temp[ref][!, key] .= value
    end
    Bp = -1.0 .* imag(calc_admittance_matrices(data_temp, indices).Ybus)
    # Compute B''
    data_temp = deepcopy(data)
    for (key, value) in zip(["br_r", "shift"], [0.0, 0.0])
        data_temp[ref][!, key] .= value
    end
    Bpp = -1.0 .* imag(calc_admittance_matrices(data_temp, indices).Ybus)

    return (Bp = Bp, Bpp = Bpp)
end
