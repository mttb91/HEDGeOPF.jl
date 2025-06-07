
function _init_channel(type::Type, size::Int)
    return _DC.RemoteChannel(() -> Channel{type}(size))
end

"Record convergence for the value of total load active power of given OPF case"
function _record_info!(data::Vector{Vector{Float64}}, case::Vector{Float64})
    push!(data, case)
    return Int(last(case))
end

"Update count for number of OPF iteration and feasible cases"
function _update_counter!(data::Vector{Int}, status::Int)
    data[end] += 1
    data[end-1] += status
    return nothing
end

"Monitor if all input samples have been processed"
function _convergence_check(data::Vector{Int}, limit::Int)
    return (data[end] < limit)
end

"Update meter defining OPF convergence rate and the remaining number of iterations."
function _update_meter!(p, counter::Vector{Int}, n_iter::Int, limit::Int)

    # Get convergence rate
    @views value = counter[end-1] / (counter[end] + n_iter) * 100
    _PMT.update!(p,
        limit - minimum(counter[[begin, end]]);
        showvalues = () -> [("Convergence rate", round(value; digits=2))]
    )
end

"Initialize remote workers importing relevant modules"
function init_workers!(settings::NamedTuple)
    path = pwd()
    base_seed = settings.CASE.baseseed

    _DC.nprocs() > 1 && _DC.rmprocs(_DC.workers())
    # Create worker processes and propagate current active environment to workers
    _DC.addprocs(
        ceil(Int, Sys.CPU_THREADS * (settings.PARALLEL.cpu_ratio / 100));
        exeflags = "--project=$(Base.active_project())"
    )
    # Run worker initialization
    init_block = quote
        using HEDGeOPF
        # Import LP and NLP solvers
        for name in getproperty.(Ref($(settings.SOLVER)), [:lp, :nlp])
            Core.eval(Main, :(import $(Symbol(name))))
        end
        cd($path)
        HEDGeOPF._RND.seed!(HEDGeOPF._DC.myid() + $base_seed)
    end
    _DC.remotecall_eval(Main, _DC.workers(), init_block)

    return nothing
end

"Sample the polytope uniformly in terms of total load active power"
function sample_polytope_uniformly(bounds::Matrix{Float64}, A::_SA.SharedMatrix{Float64}, b::Vector{Float64}, model::Tuple,
    n_samples::Int, idx::Int, rng::_RND.AbstractRNG
)

    @sync for w in _DC.workers()
        _DC.@spawnat w begin
            # Create references to mutable objects on workers
            !isdefined(Main, :GA) && (global GA = Ref{Union{_SA.SharedMatrix{Float64}, Nothing}}(nothing))
            !isdefined(Main, :GB) && (global GB = Ref{Union{Vector{Float64}, Nothing}}(nothing))
            !isdefined(Main, :GMODEL) && (global GMODEL = Ref{Union{Tuple, Nothing}}(nothing))
            # Broadcast relevant variables to workers
            GA[] = A
            GB[] = deepcopy(b)
            GMODEL[] = deepcopy(model)
        end
    end

    @views samples = _DC.pmap(x-> begin
            m, vars, cons = GMODEL[]
            # Update polytope and model right hand side
            GB[][idx-1:idx] .= x[2]
            JuMP.set_normalized_rhs.(cons[idx-1:idx], x[2])
            # Compute chebyshev centre
            x0 = chebyshev_centre(m, vars)
            return sample_polytope(GA[], GB[], x0, n_samples + 1; seed = x[1])[:, begin+1:end]
        end,
        zip(_RND.rand(rng, Int32, size(bounds, 2)), eachcol(bounds))
    )

    # Free up space
    @sync for w in _DC.workers()
        _DC.@spawnat w begin
            GA[] = nothing
            GB[] = nothing
            GMODEL[] = nothing
        end
    end
    return reduce(vcat, [copy.(eachcol(h)) for h in samples])
end

"Generate OPF cases in parallel for given load instances"
function generate_opf_samples!(counter::Vector{Int}, model::_PM.AbstractPowerModel, database::NamedTuple, samples::Vector{Dict{String, <:Any}}, 
    config::String)

    nw = _DC.nworkers()
    limit = length(samples)
    chan_to = _init_channel(Vector{Float64}, nw * 4)
    chan_fr = _init_channel(Dict{String, <:Any}, nw * 4)
    meter = _PMT.ProgressThresh(0; dt=30, showspeed=true)

    n_iter = counter[end]
    counter[end] -= n_iter
    data = Vector{Float64}[]

    @sync begin
        [put!(chan_fr, popfirst!(samples)) for _ in 1:minimum([nw * 2, length(samples)])]

        for w in _DC.workers()
            _DC.@spawnat w begin
                while true
                    try
                        put!(chan_to, solve_model!(model, database, take!(chan_fr)))
                        export_batch!(database, w, config)
                    catch
                        break
                    end
                end
                export_batch!(database, w, config; final_batch=true)
            end
        end
        while _convergence_check(counter, limit)

            status = _record_info!(data, take!(chan_to))
            _update_counter!(counter, status)
            _update_meter!(meter, counter, n_iter, limit)
            if isempty(samples)
                isopen(chan_fr) && close(chan_fr)
            else
                put!(chan_fr, popfirst!(samples))
            end
        end
    end
    isopen(chan_fr) && close(chan_fr)
    close(chan_to)
    counter[end] += n_iter
    return permutedims(reduce(hcat, data))
end
