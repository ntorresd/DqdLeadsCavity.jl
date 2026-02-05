using CairoMakie, LaTeXStrings

# global approach plot
begin
    f = Figure(size = (900, 600))
    ax = Axis(
        f[1, 1];
        xlabel = L"t_c",
        ylabel = L"I_L",
        title = L"\text{Particle current (global) vs } t_c",
        xscale = log10
    )
    lines!(ax, tc_range, I_L_num_gl; label=L"I_L^{gl,num}", color=:green, linestyle=:solid, alpha=0.7)
    lines!(ax, tc_range, I_L_ana_gl; label=L"I_L^{gl,ana}", color=:green, linestyle=:dashdot, alpha=0.7)
    lines!(ax, tc_range, I_L_neqgf; label=L"I_L^{ana}", color=:black, linestyle=:dash, alpha=0.7)
    axislegend(ax, position=:lc)
    if save_fig
        save(joinpath(@__DIR__, "I_gl.png"), f)
    else
        display(f)
    end
end

# local approach plot
begin
    f = Figure(size = (900, 600))
    ax = Axis(
        f[1, 1];
        xlabel = L"t_c",
        ylabel = L"I_L",
        title = L"\text{Particle current (local) vs } t_c",
        xscale = log10
    )
    lines!(ax, tc_range, I_L_num_loc; label=L"I_L^{loc,num}", color=:red, linestyle=:solid, alpha=0.7)
    lines!(ax, tc_range, I_L_ana_loc; label=L"I_L^{loc,ana}", color=:red, linestyle=:dashdot, alpha=0.7)
    lines!(ax, tc_range, I_L_neqgf; label=L"I_L^{ana}", color=:black, linestyle=:dash, alpha=0.7)
    axislegend(ax, position=:lc)
    if save_fig
        save(joinpath(@__DIR__, "I_loc.png"), f)
    else
        display(f)
    end
end

# local and global approaches (left)
begin
    f = Figure(size = (900, 600))
    ax = Axis(
        f[1, 1];
        xlabel = L"t_c",
        ylabel = L"I_L",
        title = L"\text{Particle current (left) vs } t_c",
        xscale = log10
    )
    lines!(ax, tc_range, I_L_num_gl; label=L"I_L^{gl,num}", color=:green, linestyle=:solid, alpha=0.7)
    lines!(ax, tc_range, I_L_num_loc; label=L"I_L^{loc,num}", color=:red, linestyle=:solid, alpha=0.7)
    lines!(ax, tc_range, I_L_neqgf; label=L"I_L^{ana}", color=:black, linestyle=:dash, alpha=0.7)
    axislegend(ax, position=:lc)
    if save_fig
        save(joinpath(@__DIR__, "I_L_gl_loc.png"), f)
    else
        display(f)
    end
end
