
function fmm_stokeslet_direct(src, targ, str, lambda)
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
    quadstr = Array{Float64}(ns)
    quadvec = Array{Float64}(3,ns)
    octstr = Array{Float64}(0)
    octvec = Array{Float64}(4,0)

    dipstr = zeros(ns)
    dipvec = zeros(2,ns)
    quadstr = zeros(ns)
    quadvec = zeros(3,ns)
    octstr = zeros(ns)
    octvec = zeros(4,ns)

    
    pottarg = Array{Float64}(nt)
    gradtarg = Array{Float64}(2,nt)
    hesstarg = Array{Float64}(3,nt)

    # First component
    for i=1:ns
        quadstr[i] = 1.0
        quadvec[1,i] =  0.0
        quadvec[2,i] = -str[2,i]
        quadvec[3,i] =  str[1,i]
    end
    fmmpars = MBHFMM2DParams(lambda,src,targ,ifcharge, ifdipole,
                             ifquad, ifoct, charge, dipstr,
                             dipvec, quadstr, quadvec, octstr,
                             octvec, iprec = 3, ifalltarg = false)
    mbhfmm2d_direct!(fmmpars,targ,ifpot,pottarg,ifgrad,
                     gradtarg,ifhess,hesstarg)
    u[1, :] = pottarg

    # Second component
    for i=1:ns
        quadvec[1,i] =  str[2,i]
        quadvec[2,i] = -str[1,i]
        quadvec[3,i] =  0.0
    end

    fmmpars = MBHFMM2DParams(lambda,src,targ,ifcharge, ifdipole,
                             ifquad, ifoct, charge, dipstr,
                             dipvec, quadstr, quadvec, octstr,
                             octvec, iprec = 3, ifalltarg = false)
    mbhfmm2d_direct!(fmmpars,targ,ifpot,pottarg,ifgrad,
                     gradtarg,ifhess,hesstarg)
    u[2, :] = pottarg
       
    return u
end

function fmm_stokeslet_targ(src, targ, str, lambda; maxnodes::Int=30, maxboxes::Int=100000)
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
    quadstr = Array{Float64}(ns)
    quadvec = Array{Float64}(3,ns)
    octstr = Array{Float64}(0)
    octvec = Array{Float64}(4,0)

    dipstr = zeros(ns)
    dipvec = zeros(2,ns)
    quadstr = zeros(ns)
    quadvec = zeros(3,ns)
    octstr = zeros(ns)
    octvec = zeros(4,ns)
    
    pottarg = Array{Float64}(nt)
    gradtarg = Array{Float64}(2,0)
    hesstarg = Array{Float64}(3,0)

    fmmpars = MBHFMM2DParams(lambda,src,targ,ifcharge, ifdipole,
                             ifquad, ifoct, charge, dipstr,
                             dipvec, quadstr, quadvec, octstr,
                             octvec, iprec = 3, ifalltarg = false)

    # First component
    for i=1:ns
        quadstr[i] = 1.0
        quadvec[1,i] =  0.0
        quadvec[2,i] = -str[2,i]
        quadvec[3,i] =  str[1,i]
    end    
    fmmstor = mbhfmm2d_form(fmmpars, maxnodes=maxnodes, maxboxes=maxboxes)    
    mbhfmm2d_targ!(fmmpars,fmmstor,targ,ifpot,pottarg,ifgrad,
                     gradtarg,ifhess,hesstarg)
    u[1, :] = pottarg

    # Second component
    for i=1:ns
        quadvec[1,i] =  str[2,i]
        quadvec[2,i] = -str[1,i]
        quadvec[3,i] =  0.0
    end
    for i=1:ns
        # Sort quadvec into stored tree structure        
        idx = fmmstor.sorted_pts.isrcsort[i]
        for j=1:3
            fmmstor.quadvecsort[j,i] = quadvec[j, idx]
        end
    end    
    mbhfmm2d_targ!(fmmpars,fmmstor,targ,ifpot,pottarg,ifgrad,
                     gradtarg,ifhess,hesstarg)
    u[2, :] = pottarg
       
    return u
end
