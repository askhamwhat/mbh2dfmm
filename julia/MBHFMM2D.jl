include("BoxFMMType.jl")
include("BoxFMMUtil.jl")

# Load OpenMP library
Libdl.dlopen("libgomp", Libdl.RTLD_GLOBAL)

type MBHFMM2DParams

    # storage of parameters which
    # describe the FMM tree and call
    
    lambda::Float64
    iprec::Int32
    ifalltarg::Bool
    src::Array{Float64,2}
    targ::Array{Float64,2}
    ifcharge::Bool
    ifdipole::Bool
    ifquad::Bool
    ifoct::Bool
    charge::Array{Float64}
    dipstr::Array{Float64}
    dipvec::Array{Float64,2}        
    quadstr::Array{Float64}
    quadvec::Array{Float64,2}        
    octstr::Array{Float64}
    octvec::Array{Float64,2}
end

function MBHFMM2DParams(lambda::Float64,src::Array{Float64,2},
                        targ::Array{Float64,2},ifcharge::Bool,
                        ifdipole::Bool,ifquad::Bool,
                        ifoct::Bool,charge::Array{Float64,1},
                        dipstr::Array{Float64,1},
                        dipvec::Array{Float64,2},
                        quadstr::Array{Float64,1},
                        quadvec::Array{Float64,2},
                        octstr::Array{Float64,1},
                        octvec::Array{Float64,2};
                        iprec::Int=3,
                        ifalltarg::Bool=false)

    iprec = convert(Int32,iprec)
    
    # constructor for MBHFMM2DParams, with a couple
    # default values

    return MBHFMM2DParams(lambda,iprec,ifalltarg,
                          src,targ,ifcharge,ifdipole,
                          ifquad,ifoct,
                          charge,dipstr,dipvec,quadstr,
                          quadvec,octstr,octvec)

end

type MBHFMM2DStorage

    # storage for FMM itself
    
    tree::BoxTree2D
    ifalltarg::Array{Int32,1}
    sorted_pts::SortedPts2D
    chargesort::Array{Float64}
    dipstrsort::Array{Float64}
    dipvecsort::Array{Float64,2}        
    quadstrsort::Array{Float64}
    quadvecsort::Array{Float64,2}        
    octstrsort::Array{Float64}
    octvecsort::Array{Float64,2}
    ifcharge::Array{Int32,1}
    ifdipole::Array{Int32,1}
    ifquad::Array{Int32,1}
    ifoct::Array{Int32,1}
    isave::Array{Int32,1}
    dsave::Array{Float64,1}
    csave::Array{Complex{Float64},1}
end

function MBHFMM2DStorage_init(fmmpars::MBHFMM2DParams,
                              tree::BoxTree2D,
                              sorted_pts::SortedPts2D)

    # convert bools to fortran flags
    
    ifalltarg = boxfmm2d_booltoint32(fmmpars.ifalltarg)
    ifcharge = boxfmm2d_booltoint32(fmmpars.ifcharge)
    ifdipole = boxfmm2d_booltoint32(fmmpars.ifdipole)
    ifquad = boxfmm2d_booltoint32(fmmpars.ifquad)
    ifoct = boxfmm2d_booltoint32(fmmpars.ifoct)        

    # if charge is used, sort it 
    
    if fmmpars.ifcharge
        chargesort = fmmpars.charge[sorted_pts.isrcsort]
    else
        chargesort = zeros(Float64,1)
    end

    if fmmpars.ifdipole
        dipstrsort = fmmpars.dipstr[sorted_pts.isrcsort]
        dipvecsort = fmmpars.dipvec[:,sorted_pts.isrcsort]
    else
        dipstrsort = zeros(Float64,1)
        dipvecsort = zeros(Float64,2,1)        
    end

    if fmmpars.ifquad
        quadstrsort = fmmpars.quadstr[sorted_pts.isrcsort]
        quadvecsort = fmmpars.quadvec[:,sorted_pts.isrcsort]
    else
        quadstrsort = zeros(Float64,1)
        quadvecsort = zeros(Float64,3,1)        
    end

    if fmmpars.ifoct
        octstrsort = fmmpars.octstr[sorted_pts.isrcsort]
        octvecsort = fmmpars.octvec[:,sorted_pts.isrcsort]
    else
        octstrsort = zeros(Float64,1)
        octvecsort = zeros(Float64,4,1)
    end

    # initialize isave, dsave, csave to short arrays

    isave = Array{Int32}(1)
    dsave = Array{Float64}(1)
    csave = Array{Complex{Float64}}(1)

    return MBHFMM2DStorage(tree,ifalltarg,sorted_pts,
                           chargesort,dipstrsort,
                           dipvecsort,quadstrsort,
                           quadvecsort,octstrsort,
                           octvecsort,ifcharge,ifdipole,
                           ifquad,ifoct,isave,dsave,csave)

end
    
    

function mbhfmm2d_form(fmmpars::MBHFMM2DParams;maxnodes::Int=30,
                       maxboxes::Int=-1)

    if maxboxes < 0
        return [], 1
    end

    # form tree
    
    tree, sorted_pts, ier = BoxTree2DMaxBoxesST(src,targ,maxboxes,
                                                maxnodes=maxnodes,
                                                ifverbose=false,
                                                ifalltarg=fmmpars.ifalltarg)


    # initialize storage
    
    fmmstor = MBHFMM2DStorage_init(fmmpars,tree,sorted_pts)
                       
    # query fmm for length of isave, dsave, csave

    lambda = fmmpars.lambda
    ier1 = zeros(Int32,1)
    iprec = fmmpars.iprec
    nlev = tree.nlev
    levelbox = tree.levelbox
    iparentbox = tree.iparentbox
    ichildbox = tree.ichildbox
    icolbox = tree.icolbox
    irowbox = tree.irowbox
    nboxes = tree.nboxes
    nblevel = tree.nblevel
    iboxlev = tree.iboxlev
    istartlev = tree.istartlev
    ifalltarg = fmmstor.ifalltarg
    localonoff = tree.localonoff
    zll = tree.zll
    blength = tree.blength
    ns = sorted_pts.ns
    srcsort = sorted_pts.srcsort
    isrcladder = sorted_pts.isrcladder
    ifcharge = fmmstor.ifcharge
    chargesort = fmmstor.chargesort
    ifdipole = fmmstor.ifdipole
    dipstrsort = fmmstor.dipstrsort
    dipvecsort = fmmstor.dipvecsort
    ifquad = fmmstor.ifquad
    quadstrsort = fmmstor.quadstrsort
    quadvecsort = fmmstor.quadvecsort
    ifoct = fmmstor.ifoct
    octstrsort = fmmstor.octstrsort
    octvecsort = fmmstor.octvecsort
    isave = fmmstor.isave
    dsave = fmmstor.dsave
    csave = fmmstor.csave
    

    lisave = [convert(Int32,-1)]
    ldsave = [convert(Int32,-1)]
    lcsave = [convert(Int32,-1)]    

    ccall( (:mbhfmm2d_form_,"../bin/libmbhfmm2d"), Void,
           (Ref{Float64},Ref{Int32},Ref{Int32},Ref{Int32},
            Ref{Int32},Ref{Int32},Ref{Int32},Ref{Int32},
            Ref{Int32},Ref{Int32},Ref{Int32},Ref{Int32},
            Ref{Int32},Ref{Int32},Ref{Int32},Ref{Float64},
            Ref{Float64},Ref{Int32},Ref{Float64},Ref{Int32},
            Ref{Int32},Ref{Float64},
            Ref{Int32},Ref{Float64},Ref{Float64},
            Ref{Int32},Ref{Float64},Ref{Float64},
            Ref{Int32},Ref{Float64},Ref{Float64},
            Ref{Int32},Ref{Int32},Ref{Float64},Ref{Int32},
            Ref{Complex{Float64}},Ref{Int32}),
           lambda,ier1,iprec,nlev,levelbox,iparentbox,
           ichildbox,icolbox,irowbox,nboxes,nblevel,
           iboxlev,istartlev,ifalltarg,localonoff,zll,
           blength,ns,srcsort,isrcladder,ifcharge,chargesort,
           ifdipole,dipstrsort,dipvecsort,ifquad,
           quadstrsort,quadvecsort,ifoct,octstrsort,
           octvecsort,isave,lisave,dsave,ldsave,csave,
           lcsave)

    # allocate

    fmmstor.isave = Array{Int32}(lisave[1])
    fmmstor.dsave = Array{Float64}(ldsave[1])
    fmmstor.csave = Array{Complex{Float64}}(lcsave[1])

    isave = fmmstor.isave
    dsave = fmmstor.dsave
    csave = fmmstor.csave

    # form fmm

    ccall( (:mbhfmm2d_form_,"../bin/libmbhfmm2d"), Void,
           (Ref{Float64},Ref{Int32},Ref{Int32},Ref{Int32},
            Ref{Int32},Ref{Int32},Ref{Int32},Ref{Int32},
            Ref{Int32},Ref{Int32},Ref{Int32},Ref{Int32},
            Ref{Int32},Ref{Int32},Ref{Int32},Ref{Float64},
            Ref{Float64},Ref{Int32},Ref{Float64},Ref{Int32},
            Ref{Int32},Ref{Float64},
            Ref{Int32},Ref{Float64},Ref{Float64},
            Ref{Int32},Ref{Float64},Ref{Float64},
            Ref{Int32},Ref{Float64},Ref{Float64},
            Ref{Int32},Ref{Int32},Ref{Float64},Ref{Int32},
            Ref{Complex{Float64}},Ref{Int32}),
           lambda,ier1,iprec,nlev,levelbox,iparentbox,
           ichildbox,icolbox,irowbox,nboxes,nblevel,
           iboxlev,istartlev,ifalltarg,localonoff,zll,
           blength,ns,srcsort,isrcladder,ifcharge,chargesort,
           ifdipole,dipstrsort,dipvecsort,ifquad,
           quadstrsort,quadvecsort,ifoct,octstrsort,
           octvecsort,isave,lisave,dsave,ldsave,csave,
           lcsave)

    return fmmstor

end


function mbhfmm2d_targ!(fmmpars::MBHFMM2DParams,
                        fmmstor::MBHFMM2DStorage,
                        targ::Array{Float64,2},
                        ifpot::Bool,pottarg::Array{Float64,1},
                        ifgrad::Bool,gradtarg::Array{Float64,2},
                        ifhess::Bool,hesstarg::Array{Float64,2})

    # unpack
    
    tree = fmmstor.tree
    sorted_pts = fmmstor.sorted_pts
    
    lambda = fmmpars.lambda
    ier1 = zeros(Int32,1)
    nlev = tree.nlev
    levelbox = tree.levelbox
    iparentbox = tree.iparentbox
    ichildbox = tree.ichildbox
    icolbox = tree.icolbox
    irowbox = tree.irowbox
    nboxes = tree.nboxes
    nblevel = tree.nblevel
    iboxlev = tree.iboxlev
    istartlev = tree.istartlev
    zll = tree.zll
    blength = tree.blength
    ns = sorted_pts.ns
    srcsort = sorted_pts.srcsort
    isrcladder = sorted_pts.isrcladder

    ifcharge = fmmstor.ifcharge
    chargesort = fmmstor.chargesort
    ifdipole = fmmstor.ifdipole
    dipstrsort = fmmstor.dipstrsort
    dipvecsort = fmmstor.dipvecsort
    ifquad = fmmstor.ifquad
    quadstrsort = fmmstor.quadstrsort
    quadvecsort = fmmstor.quadvecsort
    ifoct = fmmstor.ifoct
    octstrsort = fmmstor.octstrsort
    octvecsort = fmmstor.octvecsort
    isave = fmmstor.isave
    dsave = fmmstor.dsave
    csave = fmmstor.csave

    mtemp, nt = size(targ)
    nt = convert(Int32,nt)

    ifpot1 = boxfmm2d_booltoint32(ifpot)
    ifgrad1 = boxfmm2d_booltoint32(ifgrad)
    ifhess1 = boxfmm2d_booltoint32(ifhess)

    ifder31 = zeros(Int32,1)
    der3targ = zeros(Float64,4,1)

    ccall( (:mbhfmm2d3_targ_,"../bin/libmbhfmm2d"), Void,
           (Ref{Float64},Ref{Int32},Ref{Int32},
            Ref{Int32},Ref{Int32},Ref{Int32},Ref{Int32},
            Ref{Int32},Ref{Int32},Ref{Int32},Ref{Int32},
            Ref{Int32},Ref{Float64},
            Ref{Float64},Ref{Int32},Ref{Float64},Ref{Int32},
            Ref{Int32},Ref{Float64},
            Ref{Int32},Ref{Float64},Ref{Float64},
            Ref{Int32},Ref{Float64},Ref{Float64},
            Ref{Int32},Ref{Float64},Ref{Float64},
            Ref{Int32},Ref{Float64},Ref{Complex{Float64}},
            Ref{Int32},Ref{Float64},Ref{Int32},
            Ref{Float64},Ref{Int32},Ref{Float64},
            Ref{Int32},Ref{Float64},Ref{Int32},
            Ref{Float64}),
           lambda,ier1,nlev,levelbox,iparentbox,
           ichildbox,icolbox,irowbox,nboxes,nblevel,
           iboxlev,istartlev,zll,
           blength,ns,srcsort,isrcladder,ifcharge,chargesort,
           ifdipole,dipstrsort,dipvecsort,ifquad,
           quadstrsort,quadvecsort,ifoct,octstrsort,
           octvecsort,isave,dsave,csave,nt,targ,ifpot1,
           pottarg,ifgrad1,gradtarg,ifhess1,hesstarg,
           ifder31,der3targ)
    
    return
end

function mbhfmm2d_direct!(fmmpars::MBHFMM2DParams,
                          targ::Array{Float64,2},
                          ifpot::Bool,pottarg::Array{Float64,1},
                          ifgrad::Bool,gradtarg::Array{Float64,2},
                          ifhess::Bool,hesstarg::Array{Float64,2})


    src = fmmpars.src
    mtemp, ns = size(src)
    ns = convert(Int32,ns)
    mtemp, nt = size(targ)
    nt = convert(Int32,nt)

    lambda = fmmpars.lambda
    ifcharge = boxfmm2d_booltoint32(fmmpars.ifcharge)
    ifdipole = boxfmm2d_booltoint32(fmmpars.ifdipole)
    ifquad = boxfmm2d_booltoint32(fmmpars.ifquad)
    ifoct = boxfmm2d_booltoint32(fmmpars.ifoct)

    charge = fmmpars.charge
    dipstr = fmmpars.dipstr
    dipvec = fmmpars.dipvec
    quadstr = fmmpars.quadstr
    quadvec = fmmpars.quadvec
    octstr = fmmpars.octstr
    octvec = fmmpars.octvec

    ifpot1 = boxfmm2d_booltoint32(ifpot)
    ifgrad1 = boxfmm2d_booltoint32(ifgrad)
    ifhess1 = boxfmm2d_booltoint32(ifhess)

    pottemp = zeros(Float64,1)
    gradtemp = zeros(Float64,2)
    hesstemp = zeros(Float64,3)    

    for i = 1:nt

    ccall((:mbhpotgrad2dall_cdqo_,"../bin/libmbhfmm2d"),Void,
          (Ref{Float64},Ref{Float64},Ref{Int32},Ref{Int32},
           Ref{Float64},Ref{Int32},Ref{Float64},Ref{Float64},
           Ref{Int32},Ref{Float64},Ref{Float64},
           Ref{Int32},Ref{Float64},Ref{Float64},
           Ref{Float64},Ref{Int32},Ref{Float64},
           Ref{Int32},Ref{Float64},Ref{Int32},Ref{Float64}),
          lambda,src,ns,ifcharge,charge,ifdipole,dipstr,dipvec,
          ifquad,quadstr,quadvec,ifoct,octstr,octvec,
          targ[:,i],ifpot1,pottemp,ifgrad1,gradtemp,
          ifhess1,hesstemp)

        pottarg[i] = pottemp[1]
        gradtarg[:,i] = gradtemp
        hesstarg[:,i] = hesstemp        

    end

    return
end
