# FMM for Modified Stokes, based on FMM for Modified Biharmonic

include("MBHFMM2D.jl")

const DEFAULT_MAXNODES = 30
const DEFAULT_MAXBOXES = 100000


## STOKESLET #########################################################

function fmm_stokeslet_direct(src, targ, str, alpha)
    ns = size(src, 2)
    nt = size(targ, 2)
    u = zeros(2, nt)

    
    ifcharge = false
    ifdipole = false
    ifquad = true
    ifoct = false
    
    ifpot = true
    ifgrad = false
    ifhess = false

    charge = Array{Float64}(0)
    dipstr = Array{Float64}(0)
    dipvec = Array{Float64}(2,0)
    quadstr = ones(ns)
    quadvec = Array{Float64}(3,ns)
    octstr = Array{Float64}(0)
    octvec = Array{Float64}(4,0)
    
    pottarg = Array{Float64}(nt)
    gradtarg = Array{Float64}(2,nt)
    hesstarg = Array{Float64}(3,nt)

    fmmpars = MBHFMM2DParams(alpha,src,targ,ifcharge, ifdipole,
                             ifquad, ifoct, charge, dipstr,
                             dipvec, quadstr, quadvec, octstr,
                             octvec, iprec = 3, ifalltarg = false)
    
    for j=1:2
        ej = zeros(2)
        ej[j] = 1.0
        for i=1:ns
            f = str[:,i]
            quadvec[:, i] = f[j]*[1.0, 0.0, 1.0] - boxfmm2d_formquadvec(ej, f)
        end
        mbhfmm2d_direct!(fmmpars,targ,ifpot,pottarg,ifgrad,
                         gradtarg,ifhess,hesstarg)
        @. u[j, :] = pottarg
    end
       
    return u
end

function fmm_stokeslet_targ(src, targ, str, alpha;
                            maxnodes::Int=DEFAULT_MAXNODES,
                            maxboxes::Int=DEFAULT_MAXBOXES)
    ns = size(src, 2)
    nt = size(targ, 2)
    u = Array{Float64}(2, nt)

    
    ifcharge = false
    ifdipole = false
    ifquad = true
    ifoct = false
    
    ifpot = true
    ifgrad = false
    ifhess = false

    charge = Array{Float64}(0)
    dipstr = Array{Float64}(0)
    dipvec = Array{Float64}(2,0)
    quadstr = ones(ns)
    quadvec = Array{Float64}(3,ns)
    octstr = Array{Float64}(0)
    octvec = Array{Float64}(4,0)
    
    pottarg = Array{Float64}(nt)
    gradtarg = Array{Float64}(2,0)
    hesstarg = Array{Float64}(3,0)

    fmmpars = MBHFMM2DParams(alpha,src,targ,ifcharge, ifdipole,
                             ifquad, ifoct, charge, dipstr,
                             dipvec, quadstr, quadvec, octstr,
                             octvec, iprec = 3, ifalltarg = false)

    # Form tree    
    tree, sorted_pts, ier = mbhfmm2d_tree(fmmpars, maxnodes=maxnodes, maxboxes=maxboxes)
    # Compute components
    for j=1:2
        ej = zeros(2)
        ej[j] = 1.0
        for i=1:ns
            f = str[:,i]
            quadvec[:, i] = f[j]*[1.0, 0.0, 1.0] - boxfmm2d_formquadvec(ej, f)
        end
        fmmstor = mbhfmm2d_form(fmmpars, tree, sorted_pts)        
        mbhfmm2d_targ!(fmmpars,fmmstor,targ,ifpot,pottarg,ifgrad,
                       gradtarg,ifhess,hesstarg)        
        @. u[j, :] = pottarg
    end
       
    return u
end

## STRESSLET #########################################################

function fmm_stresslet_direct(src, targ, fvec, nvec, alpha;
                              self=false)
    ns = size(src, 2)
    nt = size(targ, 2)
    u = zeros(2, nt)
    
    ifcharge = false
    ifdipole = true
    ifquad = false
    ifoct = true
    
    ifpot = true
    ifgrad = false
    ifhess = false

    charge = Array{Float64}(0)
    dipstr = ones(ns)
    dipvec = Array{Float64}(2,ns)
    quadstr = Array{Float64}(0)
    quadvec = Array{Float64}(3,0)
    octstr = ones(ns)
    octvec = Array{Float64}(4,ns)
    
    pottarg = Array{Float64}(nt)
    gradtarg = Array{Float64}(2,nt)
    hesstarg = Array{Float64}(3,nt)

    fmmpars = MBHFMM2DParams(alpha,src,targ,ifcharge, ifdipole,
                             ifquad, ifoct, charge, dipstr,
                             dipvec, quadstr, quadvec, octstr,
                             octvec, iprec = 3, ifalltarg = false)    
    for j=1:2
        fmm_stresslet_pack_density!(fmmpars, fvec, nvec, alpha, j)
        if self
            fmmdirect = mbhfmm2d_direct_self!(fmmpars,ifpot,pottarg,ifgrad,
                                              gradtarg,ifhess,hesstarg)
        else
            fmmdirect = mbhfmm2d_direct!(fmmpars,targ,ifpot,pottarg,ifgrad,
                                         gradtarg,ifhess,hesstarg)
        end       
        u[j, :] = pottarg
    end
    
    return u
end

## FMM TARGET

function fmm_stresslet_targ(src, targ, fvec, nvec, alpha;
                            maxnodes::Int=DEFAULT_MAXNODES,
                            maxboxes::Int=DEFAULT_MAXBOXES)
    fmmpars, tree, sorted_pts = fmm_stresslet_prep(src, targ, alpha,
                                                   maxnodes=maxnodes,
                                                   maxboxes=maxboxes)
    u = fmm_stresslet_targ(fmmpars, tree, sorted_pts, fvec, nvec, alpha)
    return u
end

function fmm_stresslet_targ(fmmpars::MBHFMM2DParams,
                            tree::BoxTree2D,
                            sorted_pts::SortedPts2D,
                            fvec, nvec, alpha)
    ns = sorted_pts.ns
    nt = sorted_pts.nt
    u = Array{Float64}(2, nt)    
    ifpot = true
    ifgrad = false
    ifhess = false
    pottarg = Array{Float64}(nt)
    gradtarg = Array{Float64}(2,0)
    hesstarg = Array{Float64}(3,0)    
    for j=1:2
        fmm_stresslet_pack_density!(fmmpars, fvec, nvec, alpha, j)
        fmmstor = mbhfmm2d_form(fmmpars, tree, sorted_pts)                
        mbhfmm2d_targ!(fmmpars,fmmstor,fmmpars.targ,ifpot,pottarg,ifgrad,
                       gradtarg,ifhess,hesstarg)        
        u[j, :] = pottarg
    end
    return u
end

## FMM SELF

function fmm_stresslet_self(src, fvec, nvec, alpha;
                            maxnodes::Int=DEFAULT_MAXNODES,
                            maxboxes::Int=DEFAULT_MAXBOXES)
    fmmpars, tree, sorted_pts = fmm_stresslet_prep(src, src, alpha,
                                                   maxnodes=maxnodes,
                                                   maxboxes=maxboxes)
    u = fmm_stresslet_self(fmmpars, tree, sorted_pts, fvec, nvec, alpha)
    return u
end

function fmm_stresslet_self(fmmpars::MBHFMM2DParams,
                            tree::BoxTree2D,
                            sorted_pts::SortedPts2D,
                            fvec, nvec, alpha)
    ns = sorted_pts.ns
    nt = sorted_pts.nt
    u = Array{Float64}(2, nt)    
    ifpot = true
    ifgrad = false
    ifhess = false
    pottarg = Array{Float64}(nt)
    gradtarg = Array{Float64}(2,0)
    hesstarg = Array{Float64}(3,0)    
    for j=1:2
        fmm_stresslet_pack_density!(fmmpars, fvec, nvec, alpha, j)
        fmmstor = mbhfmm2d_form(fmmpars, tree, sorted_pts)                
        mbhfmm2d_srcsrc!(fmmpars,fmmstor,ifpot,pottarg,ifgrad,
                         gradtarg,ifhess,hesstarg)        
        u[j, :] = pottarg
    end
    return u
end

## HELPERS

function fmm_stresslet_prep(src, targ, alpha;
                            maxnodes::Int=DEFAULT_MAXNODES,
                            maxboxes::Int=DEFAULT_MAXBOXES)
    ns = size(src, 2)
    nt = size(targ, 2)
    
    ifcharge = false
    ifdipole = true
    ifquad = false
    ifoct = true
    
    charge = Array{Float64}(0)
    dipstr = ones(ns)
    dipvec = Array{Float64}(2,ns)
    quadstr = Array{Float64}(0)
    quadvec = Array{Float64}(3,0)
    octstr = ones(ns)
    octvec = Array{Float64}(4,ns)

    fmmpars = MBHFMM2DParams(alpha,src,targ,ifcharge, ifdipole,
                             ifquad, ifoct, charge, dipstr,
                             dipvec, quadstr, quadvec, octstr,
                             octvec, iprec = 3, ifalltarg = false)
    tree, sorted_pts, ier = mbhfmm2d_tree(fmmpars, maxnodes=maxnodes, maxboxes=maxboxes)
    return fmmpars, tree, sorted_pts
end


function fmm_stresslet_pack_density!(fmmpars, fvec, nvec, alpha, component)
    j = component
    ns = size(fvec, 2)
    ej = zeros(2)
    ej[j] = 1.0
    for i=1:ns
        f = fvec[:, i]
        n = nvec[:, i]
        fdotn = f[1]*n[1] + f[2]*n[2]
        fmmpars.octvec[:,i] = ( f[j]*[n[1], n[2], n[1], n[2]] +
                                n[j]*[f[1], f[2], f[1], f[2]] +
                                fdotn*[ej; ej] ) +
                                -2*boxfmm2d_formoctvec(ej, f, n)
        fmmpars.dipvec[:,i] = -alpha^2*fdotn*ej
    end    
end
   
