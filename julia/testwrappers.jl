push!(LOAD_PATH, string(pwd(),"/../../code/src"))
using ModifiedStokesSolver

include("MBHFMM2D.jl")
include("ModifiedStokesWrappers.jl")

alpha = 20.0

## define source and target points

ns = 50
nt = 30

srand(1)
stokeslet_str = 1 - 2*rand(2,ns)
src = 1 - 2*rand(2,ns)
targ = 1 - 2*rand(2,nt)

# Compute using the stokeslet kernel from ModifiedStokesSolver
function mss_stokeslet_direct(src, targ, stokeslet_str, alpha)
    ns = size(src, 2)
    nt = size(targ, 2)
    u = zeros(2, nt)
    for i=1:ns
        for j=1:nt
            f = stokeslet_str[:, i]
            r = src[:,i]-targ[:,j]
            u[:,j] += stokeslet(r, f, alpha)
        end
    end
    return u
end

u1 = mss_stokeslet_direct(src, targ, stokeslet_str, alpha)
u2 = fmm_stokeslet_direct(src, targ, stokeslet_str, alpha)
u3 = fmm_stokeslet_targ(src, targ, stokeslet_str, alpha)

errdirect = norm(u1[:]-u2[:]) / norm(u1[:], Inf)
errfmm = norm(u1[:]-u3[:]) / norm(u1[:], Inf)

@show errdirect
@show errfmm

