
using CairoMakie, LaTeXStrings

begin
    local eV = dqd_leads.leads.Δμ
    local f = Figure()
    ax = Axis(
        f[1, 1];
        xlabel = L"\Delta \epsilon",
        ylabel = L"\mathcal{C}",
        title = latexstring("\\mathcal{C} \\text{ vs } \\Delta \\epsilon \\quad (eV = $(eV))"),
        xscale = log10
    )
    lines!(ax, Δϵ_logrange, C_list_Δϵ)
    display(f)
end

begin
    local Δϵ = dqd_leads.dqd.Δϵ
    local f = Figure()
    ax = Axis(
        f[1, 1];
        xlabel = L"eV",
        ylabel = L"\mathcal{C}",
        title = latexstring("\\mathcal{C} \\text{ vs } eV \\quad (\\Delta\\epsilon = $Δϵ)"),
        xscale = log10
    )
    lines!(ax, eV_logrange, C_list_eV)
    display(f)
end

begin
    local Δϵ = dqd_leads.dqd.Δϵ
    local eV = dqd_leads.leads.Δμ
    local tc = dqd_leads.dqd.tc

    local f = Figure()
    ax = Axis(
        f[1, 1],
        xlabel = L"\Delta \epsilon",
        ylabel = L"eV",
        title = latexstring("\\text{Analytical Concurrence Heatmap - Global} (t_c=$tc)"),
        xscale = log10,
        yscale = log10
    )
    hm = heatmap!(
        ax,
        Δϵ_logrange, eV_logrange,
        C_matrix_ana;
        colormap = :viridis
    )

    vlines!(ax, [Δϵ], color = :red, linewidth = 1, linestyle = :dash)
    hlines!(ax, [eV], color = :red, linewidth = 1, linestyle = :dash)

    Colorbar(f[1, 2], hm, label = "C")
    display(f)
end

begin
    local Δϵ = dqd_leads.dqd.Δϵ
    local eV = dqd_leads.leads.Δμ
    local tc = dqd_leads.dqd.tc

    local f = Figure()
    ax = Axis(
        f[1, 1],
        xlabel = L"\Delta \epsilon",
        ylabel = L"eV",
        title = latexstring("\\text{Numerical Concurrence - Global} (t_c=$tc)"),
        xscale = log10,
        yscale = log10
    )
    hm = heatmap!(
        ax,
        Δϵ_logrange, eV_logrange,
        C_matrix_num;
        colormap = :viridis
    )

    vlines!(ax, [Δϵ], color = :red, linewidth = 1, linestyle = :dash)
    hlines!(ax, [eV], color = :red, linewidth = 1, linestyle = :dash)

    Colorbar(f[1, 2], hm, label = "C")
    display(f)
end

begin
    local Δϵ = dqd_leads.dqd.Δϵ
    local eV = dqd_leads.leads.Δμ
    local tc = dqd_leads.dqd.tc

    local f = Figure()
    ax = Axis(
        f[1, 1],
        xlabel = L"\Delta \epsilon",
        ylabel = L"eV",
        title = latexstring("\\text{Numerical Coherence - Global} (t_c=$tc)"),
        xscale = log10,
        yscale = log10
    )
    hm = heatmap!(
        ax,
        Δϵ_logrange, eV_logrange,
        α_matrix_num;
        colormap = :plasma
    )

    vlines!(ax, [Δϵ], color = :red, linewidth = 1, linestyle = :dash)
    hlines!(ax, [eV], color = :red, linewidth = 1, linestyle = :dash)

    Colorbar(f[1, 2], hm, label = L"|\alpha|")
    display(f)
end

begin
    local Δϵ = dqd_leads.dqd.Δϵ
    local eV = dqd_leads.leads.Δμ
    local tc = dqd_leads.dqd.tc

    local f = Figure()
    ax = Axis(
        f[1, 1],
        xlabel = L"\Delta \epsilon",
        ylabel = L"eV",
        title = latexstring("\\text{Numerical Double Population - Global} (t_c=$tc)"),
        xscale = log10,
        yscale = log10
    )
    hm = heatmap!(
        ax,
        Δϵ_logrange, eV_logrange,
        pD_matrix_num;
        colormap = :cividis
    )

    vlines!(ax, [Δϵ], color = :red, linewidth = 1, linestyle = :dash)
    hlines!(ax, [eV], color = :red, linewidth = 1, linestyle = :dash)

    Colorbar(f[1, 2], hm, label = L"p_D")
    display(f)
end