# Tunnel coupling dependence of the particle and heat current
# for different average DQD levels
begin
    system = deepcopy(dqd_leads)
    system.dqd.blockade = false
    system.dqd.U = 0.0

    eV = system.leads.Δμ
    N_range = 1000
    tc_range = range(1e-4, 1.25 * eV, N_range)
    ϵ_avg_list = [0.0, eV / 4., eV / 2.]

    # Figure template
    f = Figure()
    Label(
        f[0, 1:2],
        latexstring("\\Gamma = $(system.ΓL), \\quad eV = $eV, \\quad U = $(system.dqd.U), \\quad \\Delta \\epsilon = $(system.dqd.Δϵ), \\quad \\bar{T} = $(system.leads.TL)"),
        # fontsize = 26,
        tellwidth = false,  # ensures proper centering
    )
    ax1 = Axis(
        f[1, 1],
        xlabel = L"t_c",
        ylabel = L"I"
    )
    ax2 = Axis(
        f[1, 2],
        xlabel = L"t_c",
        ylabel = L"J"
    )
    lines1_gl = Lines[]
    lines2_gl = Lines[]
    lines1_neqfg = Lines[]
    lines2_neqfg = Lines[]
    
    # lines_palette = [(c, 0.6) for c in ColorSchemes.Dark2_4]
    lines_palette = [(get(ColorSchemes.gist_earth, i), 0.8) for i in range(0., 0.75, length=4)]
    labels = LaTeXString[]

    for (j, ϵ_avg) in enumerate(ϵ_avg_list)
        print(N_range - j, " steps left\n")
        system.dqd.ϵ_avg = ϵ_avg
        local I_list_gl = Vector{Float64}(undef, N_range)
        local J_list_gl = Vector{Float64}(undef, N_range)
        local I_list_neqgf = Vector{Float64}(undef, N_range)
        local J_list_neqgf = Vector{Float64}(undef, N_range)
        for (i, tc) in enumerate(tc_range)
            system.dqd.tc = tc
            # global approach
            local H_dqd = build_H_dqd_LR(system.dqd)
            local L_ops = build_L_ops_dqd_gl(system)
            local ρss = steadystate(H_dqd, L_ops)

            I_list_gl[i] = get_particle_current_num(system, ρss, L_ops[1:4])
            J_list_gl[i] = get_heat_current_num(system, ρss, L_ops[1:4], H_dqd)
            # NEQGFs approach
            I_list_neqgf[i] = get_particle_current_neqgf(system, int_lims = (-1e2 * eV, 1e2 * eV))
            J_list_neqgf[i] = get_heat_current_neqgf(system, int_lims = (-1e2 * eV, 1e2 * eV))
        end
        # --- Plot lines ---
        # global approach
        label = L"\bar\epsilon = %$ϵ_avg"
        push!(labels, label)
        l1 = lines!(ax1, tc_range, I_list_gl, color = lines_palette[j])
        push!(lines1_gl, l1)
        l2 = lines!(ax2, tc_range, J_list_gl, color = lines_palette[j])
        push!(lines2_gl, l2)
        # NEQGFs approach
        l1 = lines!(ax1, tc_range, I_list_neqgf, color = lines_palette[j], linestyle = :dash)
        push!(lines1_neqfg, l1)
        l2 = lines!(ax2, tc_range, J_list_neqgf, color = lines_palette[j], linestyle = :dash)
        push!(lines2_neqfg, l2)
    end
    Legend(
        f[2, 1:2], lines2_gl, labels, 
        orientation = :horizontal,
        halign = :center
    )
    display(f)
end

# Particle and heat current for different average temperatures
begin
    system = deepcopy(dqd_leads)
    eV = system.leads.Δμ
    N_range = 1000
    # tc_range = range(-eV - 10., eV + 10., N_range)
    tc_range = range(1e-4, 1.25 * eV, N_range)
    T_avg_list = [1e-3, 1e1, 1e2/2., 1e2]

    # Figure template
    f = Figure()
    Label(
        f[0, 1:2],
        latexstring("\\Gamma = $(system.ΓL), \\quad eV = $eV, \\quad U = $(system.dqd.U), \\quad \\Delta \\epsilon = $(system.dqd.Δϵ), \\quad \\bar{\\epsilon} = $(system.dqd.ϵ_avg)"),
        # fontsize = 26,
        tellwidth = false,  # ensures proper centering
    )
    ax1 = Axis(
        f[1, 1],
        xlabel = L"t_c",
        ylabel = L"I"
    )
    ax2 = Axis(
        f[1, 2],
        xlabel = L"t_c",
        ylabel = L"J"
    )

    lines1_gl = Lines[]
    lines2_gl = Lines[]
    lines1_neqfg = Lines[]
    lines2_neqfg = Lines[]

    lines_palette = [(get(ColorSchemes.thermal, i), 0.6) for i in range(0, 0.75, length=4)]
    labels = LaTeXString[]

    for (j, T_avg) in enumerate(T_avg_list)
        print(N_range - j, " steps left\n")
        system.leads.TL = T_avg
        system.leads.TR = T_avg
        local I_list_gl = Vector{Float64}(undef, N_range)
        local J_list_gl = Vector{Float64}(undef, N_range)
        local I_list_neqgf = Vector{Float64}(undef, N_range)
        local J_list_neqgf = Vector{Float64}(undef, N_range)
        for (i, tc) in enumerate(tc_range)
            system.dqd.tc = tc

            local H_dqd = build_H_dqd_LR(system.dqd)
            local L_ops = build_L_ops_dqd_gl(system)
            # local L_ops = build_L_ops_dqd_gl_int(system)
            local ρss = steadystate(H_dqd, L_ops)
            
            # global approach (non-interacting or Coulomb blockade)
            I_list_gl[i] = get_particle_current_num(system, ρss, L_ops[1:4])
            J_list_gl[i] = get_heat_current_num(system, ρss, L_ops[1:4], H_dqd)
            # global approach (finite Coulomb interaction)
            # I_list_gl[i] = get_particle_current_num(system, ρss, L_ops[1:8])
            # J_list_gl[i] = get_heat_current_num(system, ρss, L_ops[1:8], H_dqd)
            # NEQGFs approach
            I_list_neqgf[i] = get_particle_current_neqgf(system, int_lims = (-1e2 * eV, 1e2 * eV))
            J_list_neqgf[i] = get_heat_current_neqgf(system, int_lims = (-1e2 * eV, 1e2 * eV))

        end
        # --- Plot lines ---
        label = L"\bar{T} = %$T_avg"
        push!(labels, label)
        # global approach
        l1 = lines!(ax1, tc_range, I_list_gl, color = lines_palette[j])
        push!(lines1_gl, l1)
        l2 = lines!(ax2, tc_range, J_list_gl, color = lines_palette[j])
        push!(lines2_gl, l2)
        # NEQGFs approach
        l1 = lines!(ax1, tc_range, I_list_neqgf, color = lines_palette[j], linestyle = :dash)
        push!(lines1_neqfg, l1)
        l2 = lines!(ax2, tc_range, J_list_neqgf, color = lines_palette[j], linestyle = :dash)
        push!(lines2_neqfg, l2)
    end
    Legend(
        f[2, 1:2], lines2_gl, labels, 
        orientation = :horizontal,
        halign = :center
    )
    display(f)
end

# Particle and heat current heatmaps (tunnel coupling vs average DQD energy)
begin
    system = deepcopy(dqd_leads)
    eV = system.leads.Δμ
    N_range = 500
    tc_range = range(- 3 * eV / 2., 3 *eV / 2., N_range)
    ϵ_avg_range = range(-eV / 2., eV, N_range)

    local I_matrix = Matrix{Float64}(undef, N_range, N_range)
    local J_matrix = Matrix{Float64}(undef, N_range, N_range)
    for (j, ϵ_avg) in enumerate(ϵ_avg_range)
        print(N_range - j, " steps left\n")
        system.dqd.ϵ_avg = ϵ_avg
        for (i, tc) in enumerate(tc_range)
            system.dqd.tc = tc
            # global approach
            local H_dqd = build_H_dqd_LR(system.dqd)
            local L_ops = build_L_ops_dqd_gl(system)
            local ρss = steadystate(H_dqd, L_ops)

            I_matrix[i, j] = get_particle_current_num(system, ρss, L_ops[1:4])
            J_matrix[i, j] = get_heat_current_num(system, ρss, L_ops[1:4], H_dqd)
        end

    end
    # --- Figure ---
    f = Figure(; size = (1000, 450))
    # --- Particle current heatmap ---
    ax1 = Axis(
        f[1, 1],
        xlabel = L"t_c",
        ylabel = L"\bar\epsilon"
    )
    hm1 = heatmap!(
        ax1,
        tc_range, ϵ_avg_range,
        I_matrix;
        colormap = :viridis
    )

    Colorbar(f[1, 2], hm1, label = "I")
    # --- Heat current heatmap ---
    ax2 = Axis(
        f[1, 3],
        xlabel = L"t_c",
        ylabel = L"\bar\epsilon"
    )
    hm2 = heatmap!(
        ax2,
        tc_range, ϵ_avg_range,
        J_matrix;
        colormap = :plasma
    )

    Colorbar(f[1, 4], hm2, label = "I")
    display(f)
end
