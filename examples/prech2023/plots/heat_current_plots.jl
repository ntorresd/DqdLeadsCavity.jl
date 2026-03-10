using CairoMakie, LaTeXStrings

# global approach plot
begin
    f = Figure(size = (900, 600))
    ax = Axis(
        f[1, 1];
        xlabel = L"t_c",
        ylabel = L"J_L",
        title = L"\text{Heat current (global) vs } t_c",
        xscale = log10
    )
    lines!(ax, tc_range, JL_num_gl; label=L"J_L^{num}", color=:green, linestyle=:solid, alpha=0.7)
    lines!(ax, tc_range, JL_ana_gl; label=L"J_L^{ana}", color=:green, linestyle=:dashdot, alpha=0.7)
    lines!(ax, tc_range, JL_neqgf; label=L"J_L^{ana}", color=:black, linestyle=:dash, alpha=0.7)
    axislegend(ax, position=:lc)
    if save_fig
        save(joinpath(@__DIR__, "J_gl.png"), f)
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
        ylabel = L"J_L",
        title = L"\text{Heat current (local) vs } t_c",
        xscale = log10
    )
    lines!(ax, tc_range, JL_num_loc; label=L"J_L^{num}", color=:red, linestyle=:solid, alpha=0.7)
    lines!(ax, tc_range, JL_ana_loc; label=L"J_L^{ana}", color=:red, linestyle=:dashdot, alpha=0.7)
    lines!(ax, tc_range, JL_neqgf; label=L"J_L^{ana}", color=:black, linestyle=:dash, alpha=0.7)
    axislegend(ax, position=:lc)
    if save_fig
        save(joinpath(@__DIR__, "J_loc.png"), f)
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
        ylabel = L"J_L",
        title = L"\text{Heat current (left) vs } t_c",
        xscale = log10
    )
    lines!(ax, tc_range, JL_num_gl; label=L"J_L^{num}", color=:green, linestyle=:solid, alpha=0.7)
    lines!(ax, tc_range, JL_num_loc; label=L"J_L^{num}", color=:red, linestyle=:solid, alpha=0.7)
    lines!(ax, tc_range, JL_neqgf; label=L"J_L^{ana}", color=:black, linestyle=:dash, alpha=0.7)
    axislegend(ax, position=:lc)
    if save_fig
        save(joinpath(@__DIR__, "JL_gl_loc.png"), f)
    else
        display(f)
    end
end

# local and global approaches (right)
begin
    f = Figure(size = (900, 600))
    ax = Axis(
        f[1, 1];
        xlabel = L"t_c",
        ylabel = L"J_L",
        title = L"\text{Heat current (right) vs } t_c",
        xscale = log10
    )
    lines!(ax, tc_range, J_num_R_gl; label=L"J_R^{num}", color=:green, linestyle=:solid, alpha=0.7)
    lines!(ax, tc_range, J_num_R_loc; label=L"J_R^{num}", color=:red, linestyle=:solid, alpha=0.7)
    lines!(ax, tc_range, J_R_neqgf; label=L"J_R^{ana}", color=:black, linestyle=:dash, alpha=0.7)
    axislegend(ax, position=:lc)
    if save_fig
        save(joinpath(@__DIR__, "JR_gl_loc.png"), f)
    else
        display(f)
    end
end
