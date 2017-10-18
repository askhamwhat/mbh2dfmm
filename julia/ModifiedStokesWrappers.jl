
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
                            maxnodes::Int=30, maxboxes::Int=100000)
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
    fmmstor = mbhfmm2d_form(fmmpars, maxnodes=maxnodes, maxboxes=maxboxes)
    
    for j=1:2
        ej = zeros(2)
        ej[j] = 1.0
        for i=1:ns
            # Put directly into sorted form
            idx = fmmstor.sorted_pts.isrcsort[i]            
            f = str[:,idx]
            fmmstor.quadvecsort[:, i] = f[j]*[1.0, 0.0, 1.0] - boxfmm2d_formquadvec(ej, f)
        end
        mbhfmm2d_targ!(fmmpars,fmmstor,targ,ifpot,pottarg,ifgrad,
                       gradtarg,ifhess,hesstarg)        
        @. u[j, :] = pottarg
    end
       
    return u
end

function fmm_stresslet_direct(src, targ, fvec, nvec, alpha)
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
        ej = zeros(2)
        ej[j] = 1.0
        for i=1:ns
            f = fvec[:, i]
            n = nvec[:, i]
            fdotn = f[1]*n[1] + f[2]*n[2]
            octvec[:,i] = ( f[j]*[n[1], n[2], n[1], n[2]] +
                             n[j]*[f[1], f[2], f[1], f[2]] +
                             fdotn*[ej; ej] ) +
                             -2*boxfmm2d_formoctvec(ej, f, n)
            dipvec[:,i] = -alpha^2*fdotn*ej
        end
        mbhfmm2d_direct!(fmmpars,targ,ifpot,pottarg,ifgrad,
                         gradtarg,ifhess,hesstarg)
        @. u[j, :] = pottarg
    end
    
    return u
end

function fmm_stresslet_targ(src, targ, fvec, nvec, alpha;
                            maxnodes::Int=30, maxboxes::Int=100000)
    ns = size(src, 2)
    nt = size(targ, 2)
    u = Array{Float64}(2, nt)

    
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
    fmmstor = mbhfmm2d_form(fmmpars, maxnodes=maxnodes, maxboxes=maxboxes)
    
    for j=1:2
        ej = zeros(2)
        ej[j] = 1.0
        for i=1:ns
            idx = fmmstor.sorted_pts.isrcsort[i]                        
            f = fvec[:, idx]
            n = nvec[:, idx]
            fdotn = f[1]*n[1] + f[2]*n[2]
            fmmstor.octvecsort[:,i] = ( f[j]*[n[1], n[2], n[1], n[2]] +
                             n[j]*[f[1], f[2], f[1], f[2]] +
                             fdotn*[ej; ej] ) +
                             -2*boxfmm2d_formoctvec(ej, f, n)
            fmmstor.dipvecsort[:,i] = -alpha^2*fdotn*ej
        end
        mbhfmm2d_targ!(fmmpars,fmmstor,targ,ifpot,pottarg,ifgrad,
                       gradtarg,ifhess,hesstarg)        
        @. u[j, :] = pottarg
    end
    
    return u
end
