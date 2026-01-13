### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ 4f695af2-76ea-4ecf-bcb6-81e173e35fb8
begin
    import Pkg
    # Activate the notebooks project (expects Project.toml in ./notebooks)
    Pkg.activate(@__DIR__)
    # Ensure the local package is available even if the Manifest is missing/outdated
    Pkg.develop(path=joinpath(@__DIR__, ".."))
    Pkg.instantiate()
end

# ╔═╡ 46697b63-8b91-4777-9526-02089d46c5ad
begin
	using QuantumToolbox
	using Plots, LaTeXStrings
	using DqdLeadsCavity
end

# ╔═╡ 5e64e2ba-4560-4335-950c-7bb6e8557a40
md"""
## Structures
"""

# ╔═╡ a93bbb84-afde-46f2-9a25-114978bb8d75
md"""
### DQD
"""

# ╔═╡ 09b48a59-3d3f-4bd4-937c-86d2138312f3
begin
	# DQD parameters
	Δϵ = -17.57
	ϵ_avg = 0.0
	tc = 6.78
	γm = 0.5 		# Relaxation rate
	γphi = 3.92 	# Dephasing rate
	# DQD object
	dqd = Dqd(Δϵ, ϵ_avg, tc, γm, γphi)
end

# ╔═╡ 64f53092-6413-4630-b664-46a79d247e8e
md"""
### Leads
"""

# ╔═╡ e84b023c-0765-44af-b5ee-bc20e06049c0
begin
	# Leads' properties
	TL = 1e-7
	TR = 1e-7
	μL = 0.
	μR = 0.
	# Leads object
	leads = Leads(TL, TR, μL - μR, (μL + μR)/2.)
end

# ╔═╡ e65884cd-9bb4-4e16-864a-47c42d49453d
md"""
### DQD-Leads
"""

# ╔═╡ eecbf29e-c2d4-4d08-a0d1-a13ac1c54b58
begin
	# Transition rates
	ΓL = 1.
	ΓR = 1.
	# DQD-Leads object
	dqd_leads = DqdLeads(dqd, leads, ΓL, ΓR)	
end

# ╔═╡ 38ec0012-e7b0-4669-ab51-e9e12963f0b7
md"""
### Cavity
"""

# ╔═╡ fedc5569-5042-4521-8aca-958728f4162e
begin
	# Cavity properties
	ωc = get_Ω(dqd)
	κ = 0.337
	Tc = 1e-7
	# Coherent driving
	κin = 0.094
	Ndot = 0.01
	ωd = get_Ω(dqd)
	gJC = 0.457
	# Cavity object
	cavity = Cavity(ωc, κ, Tc, ωd, κin, Ndot, gJC, 3)
end

# ╔═╡ 4dd064f2-d017-47ba-912c-34bf17e09411
begin
	# DQD-Cavity object
	dqd_leads_cavity = DqdLeadsCavityObj(dqd_leads, cavity)
end

# ╔═╡ 450be9b9-fdb5-412f-9ba5-16a1fd843288
dqd_leads_cavity

# ╔═╡ 6db4fe9f-1277-4628-93dd-963f1fa02ebe
md"""
## Setup
"""

# ╔═╡ 85bcdf7e-7ec7-4ed3-af48-6b166ad94728
"""
Hamiltonian used in [zenelaj2022]
"""
function build_H_z2022(dqd_leads_cavity::DqdLeadsCavityObj)
	H_dqd_ge = build_H_dqd_ge(dqd_leads_cavity)
	H_drive_coh = build_H_cav_drive_coh(dqd_leads_cavity)
	H_JC = build_H_JC(dqd_leads_cavity)

	return H_dqd_ge + H_drive_coh + H_JC
end

# ╔═╡ cdb0e8d5-644e-4160-9226-3a58dadda702
"""
Global Lindblad operators for the DQD used in [zenelaj2022]
"""
function build_L_ops_z2022(dqd_leads_cavity::DqdLeadsCavityObj; side::Any = false)
	# Derived parameters
    Γg0, Γe0, Γ0g, Γ0e = get_transition_rates_ge(dqd_leads_cavity.dqd_leads; side = side)
	γm = dqd_leads_cavity.dqd_leads.dqd.γm
	γϕ = dqd_leads_cavity.dqd_leads.dqd.γϕ
	κ = dqd_leads_cavity.cavity.κ
	nth_cav = get_nth_cav(dqd_leads_cavity.cavity)

	# Derived operators
	sg, se = build_dqd_vladder_ops_ge(dqd_leads_cavity)
	σm, σp = build_dqd_ladder_ops_ge(dqd_leads_cavity)
	σz = build_dqd_σz_op(dqd_leads_cavity)
	a = build_cav_a_op(dqd_leads_cavity)

	# Dissipator
    L_ops = [
        sqrt(Γg0) * sg', sqrt(Γ0g) * sg, sqrt(Γ0e) * se, sqrt(Γe0) * se', 	# Leads
		sqrt(γm) * σm, 		# Relaxation
		sqrt(γϕ)/2 * σz, 	# Dephasing
		sqrt(κ * (nth_cav + 1)) * a, sqrt(κ * nth_cav) * a' # Cavity environment
    ]
    return(L_ops)
end

# ╔═╡ f0fe0c89-ed04-4f86-b5f1-9506c189fb6f
"Current through the DQD"
function get_particle_current(
	ρ::QuantumObject,
	dqd_leads_cavity::DqdLeadsCavityObj;
	side::Any = "left"
)
	# Derived parameters
    Γg0, Γe0, Γ0g, Γ0e = get_transition_rates_ge(dqd_leads_cavity.dqd_leads; side = side)

	# Derived operators
	sg, se = build_dqd_vladder_ops_ge(dqd_leads_cavity)
	nL, nR = build_dqd_number_ops_LR(dqd_leads_cavity)

	# Dissipator
	Lrho = Γg0 * D_sop(sg', ρ) + Γ0g * D_sop(sg, ρ) + Γe0 * D_sop(se', ρ) + Γ0e * D_sop(se, ρ);

	return real(expect(nL + nR, Lrho))
end

# ╔═╡ c17b9e83-1300-4ff7-a349-8b37b549065f
md"""
## Particle current
"""

# ╔═╡ 755dc2d4-e7cf-42de-9ed5-ed117148e5a8
md"""
### Large-drive limit
"""

# ╔═╡ e699d93d-5902-4ba6-a7e3-113aa04cc9a6
"""
Current in the large drive limit and the zero temperature limit
"""
function I_large_drive_ana(dqd_leads::DqdLeads)
	# Parameters
	ΓL, ΓR = dqd_leads.ΓL, dqd_leads.ΓR
    Ω = get_Ω(dqd_leads.dqd);
    θ = get_θ(dqd_leads.dqd);
	
    ΓLin = ΓL * cos(θ/2)^2;
    ΓRin = ΓR * sin(θ/2)^2;
    ΓLout = ΓL * sin(θ/2)^2;
    ΓRout = ΓR * cos(θ/2)^2;

    Γg0 = ΓLin + ΓRin;
    Γ0e = ΓLout + ΓRout;

    return (ΓLin * ΓRout - ΓRin * ΓLout)/(Γ0e + 2. * Γg0)
end

# ╔═╡ 417bb6df-e57b-4e90-9aeb-5facc0c3fbae
let
	Ndots = logrange(1e-3, 1000, 200);
	IL = [];
	IL_ld = []; 	# large drive
	IL_ld_ana = [];

	dqdObj = deepcopy(dqd_leads_cavity)
	dqdObj_ld = deepcopy(dqd_leads_cavity)
	dqdObj_ld.cavity.Ndot = 5000
	for ndot in Ndots
		dqdObj.cavity.Ndot = ndot
		# Numeric current
		ρss = steadystate(build_H_z2022(dqdObj), build_L_ops_z2022(dqdObj))
		push!(IL, get_particle_current(ρss, dqdObj; side = "left"))

		# Large drive numeric current
		ρss_ld = steadystate(build_H_z2022(dqdObj_ld), build_L_ops_z2022(dqdObj_ld))
		push!(IL_ld, get_particle_current(ρss_ld, dqdObj_ld; side = "left"))

		# Large drive analytical current
		push!(IL_ld_ana, I_large_drive_ana(dqdObj_ld.dqd_leads))
	end

	# Plot
	current_particle_plot = plot(
		sqrt.(dqdObj.cavity.κin * Ndots), [IL, IL_ld, IL_ld_ana],
	    xlabel = L"\sqrt{κ_{in}\dot{N}}",
	    ylabel = L"I",
		label = [L"I_{L}" L"I_{L}^{ld}" L"I_{L, ana}^{ld}"],
	    linestyle = [:solid :dash :dashdot],
		xaxis = :log,
	    legend = :right,
	    dpi = 200
	)
end

# ╔═╡ c767f1d5-1137-498d-95e7-80d33bb01b5c
md"""
## Cavity average occupation
"""

# ╔═╡ 73b92325-ee3c-438f-ae6c-f61bd26d80be
"""
Steady-state photton number in the cavity without DQD
"""
function get_n_cav_ss(cavity::Cavity)
	Δc = get_Δc(cavity)
	return abs(2. * sqrt(cavity.κin * cavity.Ndot)/(cavity.κ + 2im * Δc))^2
end

# ╔═╡ 9421e4bc-6da7-4cd1-8a5f-5af1955be31d
let
	Ndots = range(0, 25, 1000);
	dqdObj = deepcopy(dqd_leads_cavity)
	
	n_cav_ss_num = []
	n_cav_ss_ana = []
	for ndot in Ndots
		dqdObj.cavity.Ndot = ndot
		ρss = steadystate(build_H_z2022(dqdObj), build_L_ops_z2022(dqdObj))	
		a = build_cav_a_op(dqdObj; disp = true)

		push!(n_cav_ss_num, expect(a' * a, ρss))
		push!(n_cav_ss_ana, get_n_cav_ss(dqdObj.cavity))
	end

	n_cav_ss_plot = plot(
		Ndots, [n_cav_ss_num, n_cav_ss_ana],
	    xlabel = L"\dot{N}",
	    ylabel = L"\left< a^\dagger a \right>",
		label = ["With DQD" "Without DQD"],
	    linestyle = [:dash :dash],
	    legend = :right,
	    dpi = 200
	)
end

# ╔═╡ Cell order:
# ╠═4f695af2-76ea-4ecf-bcb6-81e173e35fb8
# ╠═46697b63-8b91-4777-9526-02089d46c5ad
# ╟─5e64e2ba-4560-4335-950c-7bb6e8557a40
# ╟─a93bbb84-afde-46f2-9a25-114978bb8d75
# ╠═09b48a59-3d3f-4bd4-937c-86d2138312f3
# ╟─64f53092-6413-4630-b664-46a79d247e8e
# ╠═e84b023c-0765-44af-b5ee-bc20e06049c0
# ╟─e65884cd-9bb4-4e16-864a-47c42d49453d
# ╠═eecbf29e-c2d4-4d08-a0d1-a13ac1c54b58
# ╟─38ec0012-e7b0-4669-ab51-e9e12963f0b7
# ╠═fedc5569-5042-4521-8aca-958728f4162e
# ╠═4dd064f2-d017-47ba-912c-34bf17e09411
# ╠═450be9b9-fdb5-412f-9ba5-16a1fd843288
# ╟─6db4fe9f-1277-4628-93dd-963f1fa02ebe
# ╠═85bcdf7e-7ec7-4ed3-af48-6b166ad94728
# ╠═cdb0e8d5-644e-4160-9226-3a58dadda702
# ╠═f0fe0c89-ed04-4f86-b5f1-9506c189fb6f
# ╟─c17b9e83-1300-4ff7-a349-8b37b549065f
# ╟─755dc2d4-e7cf-42de-9ed5-ed117148e5a8
# ╠═e699d93d-5902-4ba6-a7e3-113aa04cc9a6
# ╠═417bb6df-e57b-4e90-9aeb-5facc0c3fbae
# ╟─c767f1d5-1137-498d-95e7-80d33bb01b5c
# ╠═73b92325-ee3c-438f-ae6c-f61bd26d80be
# ╠═9421e4bc-6da7-4cd1-8a5f-5af1955be31d
