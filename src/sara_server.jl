using Oxygen
using HTTP
using JSON
using SARA
using SARA: TemperatureProfile
using SARA_Crystal
using SARA_Crystal: STGSettings
using GaussianDistributions: Gaussian
using CrystalShift
using CrystalShift: OptimizationSettings, FixedPseudoVoigt
using CrystalTree: TreeSearchSettings
using Plots

mutable struct ServerState
    cs
    SARA_Tprofile# ::TemperatureProfile
    SARA_Crystal_Tprofile# ::TemperatureProfile
    kernel# ::Gaussian # = GradientModel and XshiftModel
    relevant_T
    log_composition_proximity_likelihood


    gradient_conditional_model
    xshift_conditional_model
    char_kernel # Inner loop kernel
    xshift_char_kernel

    gradient_policy
    xshift_policy

    opt_settings
    tree_search_settings
    stripe_to_global_settings
    iter
end

global server_state = ServerState("", "", "", "", "",
                                  "", "", "", "", "",
                                  "", "", "", "", "", 1)

@put "/set_crystal_phases" function (req::HTTP.Request)
    dict = JSON.parse(String(req.body))
    path = dict["path"]
    width = dict["width"]
    α = dict["α"]
    
    open(path, "r") do f
        server_state.cs = CrystalPhase(f, width, FixedPseudoVoigt(α))
    end
    println(server_state.cs)
    println(length(server_state.cs))

    return
end

@put "/Sara_get_temperature_profile_CHESS_Spring_2024" function (req::HTTP.Request)
    dict = JSON.parse(String(req.body))
    dT = convert(Float64, dict["dT"])
    dx = convert(Float64, dict["dx"])
    server_state.SARA_Tprofile = SARA.get_temperature_profile_CHESS_Spring_2024(dT, dx)

    return
end


@put "/SARA_Crystal_get_temperature_profile_CHESS_Spring_2024" function (req::HTTP.Request)
    dict = JSON.parse(String(req.body))
    dT = convert(Float64, dict["dT"])
    dx = convert(Float64, dict["dx"])
    server_state.SARA_Crystal_Tprofile = SARA_Crystal.get_temperature_profile_CHESS_Spring_2024(dT, dx)

    return 
end


@put "/Sara_get_temperature_profile_CHESS_Spring_2025" function (req::HTTP.Request)
    dict = JSON.parse(String(req.body))
    dT = convert(Float64, dict["dT"])
    dx = convert(Float64, dict["dx"])
    server_state.SARA_Tprofile = SARA.get_temperature_profile_CHESS_Spring_2025(dT, dx)

    return
end


@put "/SARA_Crystal_get_temperature_profile_CHESS_Spring_2025" function (req::HTTP.Request)
    dict = JSON.parse(String(req.body))
    dT = convert(Float64, dict["dT"])
    dx = convert(Float64, dict["dx"])
    server_state.SARA_Crystal_Tprofile = SARA_Crystal.get_temperature_profile_CHESS_Spring_2025(dT, dx)

    return 
end


@put "/Sara_get_relevant_T" function (req::HTTP.Request)
    dict = JSON.parse(String(req.body))
    c_min = convert(Float64, dict["c_min"])
    c_max = convert(Float64, dict["c_max"])
    nout = convert(Int64, dict["nout"])
    server_state.relevant_T = SARA.get_relevant_T(c_min, c_max, nout)

    return 
end


@put "/Sara_oSARA" function (req::HTTP.Request)
    dict = JSON.parse(String(req.body))
    length_scales = convert.(Float64, dict["length_scales"])
    server_state.kernel = Gaussian(SARA.oSARA(length_scales))
    
    return
end


@put "/Sara_get_log_data_proximity_likelihood" function (req::HTTP.Request)
    dict = JSON.parse(String(req.body))
    mesh_points = [convert.(Float64, i) for i in dict["mesh_points"]]
    sigma2 = convert(Float64, dict["sigma2"])
    thresh = convert(Float64, dict["thresh"])
    alpha = convert(Float64, dict["alpha"])
    server_state.log_composition_proximity_likelihood = SARA.get_log_data_proximity_likelihood(mesh_points,
                                               sigma2,
                                               thresh,
                                               alpha)
    return 
end


@put "/Sara_CHESS_Spring_2024_sampling_gradient" function (req::HTTP.Request)
    dict = JSON.parse(String(req.body))
    lower_bounds = convert.(Float64, dict["lower_bounds"])
    upper_bounds = convert.(Float64, dict["upper_bounds"])
    maxiter = convert(Int64, dict["maxiter"])
    verbose = convert(Bool, dict["verbose"])
    server_state.gradient_policy = SARA.CHESS_Spring_2024_sampling(lower_bounds, upper_bounds,
                  server_state.relevant_T,
                  log_composition_proximity_likelihood = server_state.log_composition_proximity_likelihood,
                  maxiter=maxiter,
                  verbose=verbose,
                  acquisition=SARA.log_uncertainty_acquisition(),
                  )

    return
end


@put "/Sara_CHESS_Spring_2024_sampling_xshift" function (req::HTTP.Request)
    dict = JSON.parse(String(req.body))

    lower_bounds = convert.(Float64, dict["lower_bounds"])
    upper_bounds = convert.(Float64, dict["upper_bounds"])
    best_f = convert(Float64, dict["best_f"])
    maximize = convert(Bool, dict["maximize"])
    maxiter = convert(Int64, dict["maxiter"])
    verbose = convert(Bool, dict["verbose"])
    acquisition = SARA.log_ei_acquisition(best_f, maximize=maximize)
    server_state.xshift_policy = SARA.CHESS_Spring_2024_sampling(
                lower_bounds, upper_bounds, server_state.relevant_T,
                maxiter = maxiter,
                verbose = verbose,
                acquisition = acquisition,
                log_composition_proximity_likelihood = server_state.log_composition_proximity_likelihood)
    return 
end


@put "/Sara_CHESS_Spring_2025_sampling_gradient" function (req::HTTP.Request)
    dict = JSON.parse(String(req.body))
    lower_bounds = convert.(Float64, dict["lower_bounds"])
    upper_bounds = convert.(Float64, dict["upper_bounds"])
    maxiter = convert(Int64, dict["maxiter"])
    verbose = convert(Bool, dict["verbose"])
    composition_mesh_points = [convert.(Float64, i) for i in dict["composition_mesh_points"]]
    server_state.gradient_policy = SARA.CHESS_Spring_2025_sampling(lower_bounds, upper_bounds,
                  server_state.relevant_T,
                  log_composition_proximity_likelihood = server_state.log_composition_proximity_likelihood,
                  maxiter=maxiter,
                  verbose=verbose,
                  acquisition=SARA.log_uncertainty_acquisition(),
                  composition_mesh_points = composition_mesh_points
                  )

    return
end



@put "/Sara_CHESS_Spring_2025_sampling_xshift" function (req::HTTP.Request)
    dict = JSON.parse(String(req.body))

    lower_bounds = convert.(Float64, dict["lower_bounds"])
    upper_bounds = convert.(Float64, dict["upper_bounds"])
    best_f = convert(Float64, dict["best_f"])
    maximize = convert(Bool, dict["maximize"])
    maxiter = convert(Int64, dict["maxiter"])
    verbose = convert(Bool, dict["verbose"])
    composition_mesh_points = [convert.(Float64, i) for i in dict["composition_mesh_points"]]
    acquisition = SARA.log_ei_acquisition(best_f, maximize=maximize)
    server_state.xshift_policy = SARA.CHESS_Spring_2025_sampling(
                lower_bounds, upper_bounds, server_state.relevant_T,
                maxiter = maxiter,
                verbose = verbose,
                acquisition = acquisition,
                log_composition_proximity_likelihood = server_state.log_composition_proximity_likelihood,
                composition_mesh_points = composition_mesh_points
               )
    return 
end


@put "/Sara_get_conditional_model_XShiftModel" function (req::HTTP.Request)
    dict = JSON.parse(String(req.body))
    xshift_conditions_all = [convert.(Float64, i) for i in dict["xshift_conditions_all"]]
    xshift_data_all = convert.(Float64, dict["xshift_data_all"])
    xshift_uncert_all = convert.(Float64, dict["xshift_uncert_all"])
    server_state.xshift_conditional_model = SARA.get_conditional_model(
                                                            server_state.kernel,
                                                            xshift_conditions_all,
                                                            xshift_data_all,
                                                            xshift_uncert_all
                                                           )
    return 
end


@get "/Sara_evaluate_XShiftConditionalModel" function (req::HTTP.Request)
    dict = JSON.parse(String(req.body))
    xshift_evaluate_conditions = [convert.(Float64, i) for i in dict["xshift_evaluate_conditions"]]
    μ, σ = SARA.evaluate(server_state.xshift_conditional_model,
                         xshift_evaluate_conditions)
    d = Dict{String, Vector{Float64}}()
    d["mean"] = μ
    d["std"] = σ
    return d
end


@put "/Sara_evaluate_GradientConditionalModel" function (req::HTTP.Request)
    dict = JSON.parse(String(req.body))
    gradient_evaluate_conditions = dict["gradient_evaluate_conditions"]
    μ, σ =  SARA.evaluate(server_state.gradient_conditional_model,
                          gradient_evaluate_conditions)
    d = Dict{String, Vector{Float64}}()
    d["mean"] = μ
    d["std"] = σ
    return d
end


@put "/Sara_iSARA_April2022_gradient" function (req::HTTP.Request)
    dict = JSON.parse(String(req.body))

    length_scales = convert(Float64, dict["gradient_lengthscale"])
    rescale_parameters = Tuple(convert.(Float64, dict["gradient_rescale_parameters"]))
    server_state.char_kernel = SARA.iSARA_April2022(length_scales, rescale_parameters)

    return
end


@put "/Sara_iSARA_April2022_xshift" function (req::HTTP.Request)
    dict = JSON.parse(String(req.body))
    length_scales =      convert(Float64, dict["xshift_lengthscale"])
    rescale_parameters = Tuple(convert.(Float64, dict["xshift_rescale_parameters"]))
    server_state.xshift_char_kernel = SARA.iSARA_April2022(length_scales, rescale_parameters)

    return
end


@get "/Sara_xrd_to_global" function (req::HTTP.Request)
    dict = JSON.parse(String(req.body))

    x = convert.(Float64, dict["x"])
    y = reduce(hcat, [convert.(Float64, i) for i in dict["Y"]])'
    sigma = convert(Float64, dict["sigma"])
    condition = [convert.(Float64, i) for i in dict["condition"]]
    verbose = convert(Bool, dict["verbose"])
    conditions, gradient, gradient_uncert = SARA.xrd_to_global(x, y, sigma,
                                                    server_state.char_kernel,
                                                    server_state.SARA_Tprofile,
                                                    condition,
                                                    server_state.relevant_T,
                                                    Val(verbose))
    d = Dict{String, Any}()
    d["conditions"] = conditions
    d["gradient"] = gradient
    d["gradient_uncert"] = gradient_uncert

    return d
end


@get "/Sara_next_stripe_GradientModel" function (req::HTTP.Request)
    dict = JSON.parse(String(req.body))
    conditions         = [convert.(Float64, i) for i in dict["conditions_all"]]
    graident           = convert.(Float64, dict["gradient_all"])
    gradient_uncert    = convert.(Float64, dict["gradient_uncert_all"])
    allowed_conditions = [convert.(Float64, i) for i in dict["allowed_conditions"]]
    # println(server_state.kernel)
    # println(conditions)
    # println(graident)
    # println(gradient_uncert)
    # println(allowed_conditions)
    # println(server_state.gradient_policy)
    next_cond, server_state.gradient_conditional_model = SARA.next_stripe(server_state.kernel,
                                                                 conditions,
                                                                 graident,
                                                                 gradient_uncert,
                                                                 allowed_conditions,
                                                                 server_state.gradient_policy)
    
    d = Dict{String, Any}()
    d["next_cond"] = next_cond
    return d
end


@get "/Sara_next_stripe_XShiftModel" function (req::HTTP.Request)
    dict = JSON.parse(String(req.body))
    conditions         = [convert.(Float64, i) for i in dict["conditions_all"]]
    xshift             = convert.(Float64, dict["xshift_all"])
    xshift_uncert      = convert.(Float64, dict["xshift_uncert_all"])
    allowed_conditions = [convert.(Float64, i) for i in dict["allowed_conditions"]]
    # println(typeof(server_state.kernel))
    # println(typeof(conditions))
    # println(typeof(xshift))
    # println(typeof(xshift_uncert))
    # println(typeof(allowed_conditions))
    # println(typeof(server_state.gradient_policy), server_state.gradient_policy)
    next_cond, server_state.xshift_conditional_model = SARA.next_stripe(
                                                  server_state.kernel,
                                                  conditions,
                                                  xshift,
                                                  xshift_uncert,
                                                  allowed_conditions,
                                                  server_state.xshift_policy)

    println("Next cond: ", next_cond)
    d = Dict{String, Any}()
    d["next_cond"] = next_cond
    d
end


@get "/Sara_project_on_available_compositions" function (req::HTTP.Request)
    dict = JSON.parse(String(req.body))
    next_cond = convert.(Float64, dict["next_cond"])
    allowed_compositions = reduce(hcat, [convert.(Float64, i) for i in dict["allowed_compositions"]])'
    next_cond = SARA.project_on_available_compositions(next_cond, allowed_compositions)

    d = Dict{String, Any}()
    d["next_cond"] = next_cond
    d
end


@put "/SARA_Crystal_OptimizationSettings" function (req::HTTP.Request)
    dict = JSON.parse(String(req.body))

    std_noise = convert(Float64, dict["std_noise"])
    mean = convert.(Float64, dict["mean"])
    std = convert.(Float64, dict["std"])
    maxiter = convert.(Int64, dict["maxiter"])
    regularization = convert(Bool, dict["regularization"])
    opt_method = dict["opt_method"]
    loss_func = dict["loss_func"]
    opt_mode = dict["opt_mode"]
    em_loop_num = convert(Int64, dict["em_loop_num"])

    server_state.opt_settings = OptimizationSettings{Float64}(std_noise,
                                                             mean,
                                                             std,
                                                             maxiter,
                                                             regularization,
                                                             opt_method_helper(opt_method),
                                                             loss_func,
                                                             opt_mode_helper(opt_mode),
                                                             em_loop_num)
    return
end

function opt_mode_helper(s::String)
    if s == "Simple"
        return Simple
    end
    if s == "EM"
        return EM
    end
    throw(ArgumentError("Unknown optimization mode string: $(s)"))
end


function opt_method_helper(s::String)
    if s == "LM"
        return LM
    end
    if s == "Newton"
        return Newton
    end
    if s == "bfgs"
        return bfgs
    end
    if s == "l_bfgs"
        return l_bfgs
    end
    throw(ArgumentError("Unknown optimization method string: $(s)"))
end


@put "/SARA_Crystal_TreeSearchSettings" function (req::HTTP.Request)
    dict = JSON.parse(String(req.body))

    depth            = convert(Int64, dict["depth"])
    if dict["k"] isa Vector
        k            = convert.(Int64, dict["k"])
    else
        k            = convert.(Int64, dict["k"])
    end
    check_amorphous  = convert(Bool, dict["tree_check_amorphous"])
    model_background = convert(Bool, dict["tree_model_background"])
    bg_length        = convert(Float64, dict["tree_bg_length"])

    server_state.tree_search_settings = TreeSearchSettings{Float64}(depth,
                                                                    k,
                                                                    check_amorphous,
                                                                    model_background,
                                                                    bg_length,
                                                                    server_state.opt_settings
                                                                )
    return
end


@put "/SARA_Crystal_STGSettings" function (req::HTTP.Request)
    dict = JSON.parse(String(req.body))

    rank            = convert(Int64, dict["rank"])
    h_threshold     = convert(Float64, dict["h_threshold"])
    frac_threshold  = convert(Float64, dict["frac_threshold"])
    n_top_node      = convert(Int64, dict["n_top_node"])
    bg_length       = convert(Float64, dict["length"])
    norm_constant   = convert(Float64, dict["norm_constant"])
    sigma           = convert(Float64, dict["xshift_sigma"])
    condition       = Tuple(convert.(Float64, dict["condition"]))
    check_amorphous = convert(Bool, dict["check_amorphous"])
    additional      = dict["additional"]

    server_state.stripe_to_global_settings = STGSettings(rank,
                                                         h_threshold,
                                                         frac_threshold,
                                                         n_top_node,
                                                         bg_length,
                                                         norm_constant,
                                                         server_state.char_kernel,
                                                         sigma,
                                                         server_state.SARA_Crystal_Tprofile,
                                                         condition,
                                                         check_amorphous,
                                                         Val(true),
                                                         additional) # save_plot filename which does not work
    return
end

@get "/SARA_Crystal_expected_fraction_to_global" function (req::HTTP.Request)
    dict = JSON.parse(String(req.body))
    x = convert.(Float64, dict["x"])
    Qval = convert.(Float64, dict["Qval"])
    Y = reduce(hcat, [convert.(Float64, i) for i in  dict["Y"]])'
    phase_idx = convert.(Int64, dict["phases_idx"]) .+ 1
    # println(server_state.cs)
    # println(server_state.cs[phase_idx])
    conditions, expected_fraction, expected_fraction_uncertainty, phase_fraction = SARA_Crystal.expected_fraction_to_global(
                                                                                                                            x, Qval, Y, server_state.cs[phase_idx],
                         server_state.tree_search_settings,
                         server_state.stripe_to_global_settings,
                         server_state.relevant_T)
    # heatmap(1:size(expected_fraction,2),
    #         getindex.(conditions, 1),
    #         expected_fraction, dpi = 300,
    #         xticks=(collect(1:length(phase_idx)),
    #                 getproperty.(server_state.cs[phase_idx], :name)))

    t = reverse(getindex.(conditions, 1))
    heatmap(t,# size(expected_fraction,1),
            1:size(expected_fraction,2),
            expected_fraction', dpi = 300,
            yticks=(collect(1:length(phase_idx)+1),
                    vcat(getproperty.(server_state.cs[phase_idx], :name), "Amorphous")),
           yflip=true, xflip=true)
    plot!(xlabel="Temperature (C)")
    savefig("iter_$(lpad(server_state.iter, 3, "0")).png")
    server_state.iter += 1

    println("conditions", conditions)
    println("expected_fraction", expected_fraction)

    d = Dict{String, Any}() 
    d["conditions"] = conditions
    d["xshift_activation"] = [v for v in eachrow(expected_fraction)]
    d["xshift_activation_uncert"] = [v for v in eachrow(expected_fraction_uncertainty)]
    d["xshift_phasefraction"] = [v for v in eachrow(phase_fraction)]
    return d
end

serve(port=4200)
# end # module sara_server
