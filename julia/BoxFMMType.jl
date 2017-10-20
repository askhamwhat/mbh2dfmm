type BoxTree2D

    # storage for the description of
    # the fmm tree
    
    nboxes::Int32
    nlev::Int32
    levelbox::Array{Int32,1}
    icolbox::Array{Int32,1}
    irowbox::Array{Int32,1}
    iparentbox::Array{Int32,1}    
    ichildbox::Array{Int32,2}
    istartlev::Array{Int32,1}
    iboxlev::Array{Int32,1}
    icolleagbox::Array{Int32,2}
    nblevel::Array{Int32,1}
    neighbors::Array{Int32,2}
    nnbrs::Array{Int32,1}
    localonoff::Array{Int32,1}
    zll::Array{Float64,1}
    blength::Float64
end

type SortedPts2D

    # storage for points and sorted
    # versions/ladder structures for
    # tree
    
    ns::Int32
    nt::Int32
    src::Array{Float64,2}
    targ::Array{Float64,2}
    srcsort::Array{Float64,2}
    targsort::Array{Float64,2}
    isrcsort::Array{Int32,1}
    itargsort::Array{Int32,1}    
    isrcladder::Array{Int32,2}
    itargladder::Array{Int32,2}

end

function SortedPts2D(src_in::Array{Float64,2},
                     targ_in::Array{Float64,2},
                     maxboxes)

    # constructor for sorted points storage
    # based on sources and targets
    
    src = copy(src_in)
    targ = copy(targ_in)
    srcsort = copy(src)
    targsort = copy(targ)

    m,ns = size(src)
    m,nt = size(targ)

    ns = convert(Int32,ns)
    nt = convert(Int32,nt)

    isrcsort = Array{Int32}(ns)
    itargsort = Array{Int32}(nt)

    isrcladder = Array{Int32}(2,maxboxes)
    itargladder = Array{Int32}(2,maxboxes)

    return SortedPts2D(ns,nt,src,targ,srcsort,
                       targsort,isrcsort,
                       itargsort,isrcladder,
                       itargladder)

end

function BoxTree2DMaxBoxes(maxboxes,maxlev)

    # construct an empty tree with enough
    # storage for maxboxes boxes and maxlev
    # levels
    
    nboxes = -1
    nlev = -1
    levelbox = Array{Int32}(maxboxes)
    icolbox = Array{Int32}(maxboxes)
    irowbox = Array{Int32}(maxboxes)
    iparentbox = Array{Int32}(maxboxes)
    ichildbox = Array{Int32}(4,maxboxes)
    istartlev = Array{Int32}(maxlev+1)
    iboxlev = Array{Int32}(maxboxes)    
    icolleagbox = Array{Int32}(9,maxboxes)
    nblevel = Array{Int32}(maxlev+1)
    neighbors = Array{Int32}(12,maxboxes)
    nnbrs = Array{Int32}(maxboxes)
    localonoff = Array{Int32}(maxboxes)
    zll = Array{Float64}(2)
    blength = zero(Float64)

    return BoxTree2D(nboxes,nlev,levelbox,icolbox,
                     irowbox,iparentbox,ichildbox,istartlev,
                     iboxlev,icolleagbox,nblevel,
                     neighbors,nnbrs,localonoff,zll,
                     blength)

end

function BoxTree2DMaxBoxesST(src::Array{Float64,2},
                              targ::Array{Float64,2},
                              maxboxes;maxlev=45,
                              maxnodes=40,
                              ifverbose::Bool=true,
                              ifalltarg::Bool=false)

    # build a BoxTree2D based on src and target
    # locations
    #
    # Input:
    #
    # src - source locations
    # targ - target locations
    # maxboxes - maximum number of boxes in
    #            tree (determines how much
    #            memory is allocated)
    # maxnodes - maximum number of sources and
    #            targets in a box (total)
    # maxlev - maximum depth of tree
    #
    # Output:
    #
    # tree - level-restricted tree adapted to
    #        source and target locations
    # sorted_pts - src and targets sorted into
    #              tree
    # ier - flag indicating the success of the
    #       routine

    maxnodes = convert(Int32,maxnodes)
    maxlev = convert(Int32,maxlev)
    maxboxes = convert(Int32,maxboxes)

    ier = [zero(Int32)]

    @printf "sorting points ...\n"
    
    sorted_pts = SortedPts2D(src,targ,maxboxes)

    @printf "allocating memory ...\n"
    
    tree = BoxTree2DMaxBoxes(maxboxes,maxlev)

    levelbox = tree.levelbox
    icolbox = tree.icolbox
    irowbox = tree.irowbox    
    iparentbox = tree.iparentbox
    ichildbox = tree.ichildbox
    nblevel = tree.nblevel
    iboxlev = tree.iboxlev
    istartlev = tree.istartlev
    src = sorted_pts.src
    srcsort = sorted_pts.srcsort
    isrcsort = sorted_pts.isrcsort
    isrcladder = sorted_pts.isrcladder
    ns = sorted_pts.ns
    targ = sorted_pts.targ
    targsort = sorted_pts.targsort
    itargsort = sorted_pts.itargsort
    itargladder = sorted_pts.itargladder
    nt = sorted_pts.nt
    localonoff = tree.localonoff

    xmin = min(minimum(src[1,:]),minimum(targ[1,:]))
    xmax = max(maximum(src[1,:]),maximum(targ[1,:]))    
    ymin = min(minimum(src[2,:]),minimum(targ[2,:]))
    ymax = max(maximum(src[2,:]),maximum(targ[2,:]))

    copy!(tree.zll,[xmin;ymin])
    tree.blength = max(xmax-xmin,ymax-ymin)

    zll = tree.zll
    blength = tree.blength

    itemparray = Array{Int32}(maxboxes)
    nboxes = Array{Int32}(1)
    nlev = Array{Int32}(1)

    @printf "calling tree building routine ...\n"
    
    ccall( (:lrt2d_mktst_,LIBMBHFMM2D), Void,
           (Ref{Int32},Ref{Int32},Ref{Int32},
            Ref{Int32},Ref{Int32},Ref{Int32},
            Ref{Int32},Ref{Int32},Ref{Int32},
            Ref{Int32},Ref{Int32},Ref{Int32},            
            Ref{Int32},Ref{Float64},Ref{Float64},
            Ref{Int32},Ref{Int32},Ref{Int32},
            Ref{Float64},Ref{Float64},
            Ref{Int32},Ref{Int32},Ref{Int32},
            Ref{Int32},Ref{Float64},Ref{Float64},
            Ref{Int32},Ref{Int32}),levelbox,
           icolbox,irowbox,nboxes,nlev,
           iparentbox,ichildbox,nblevel,iboxlev,
           istartlev,maxboxes,itemparray,
           maxlev,src,srcsort,isrcsort,
           isrcladder,ns,targ,targsort,itargsort,
           itargladder,nt,maxnodes,zll,blength,
           ier,localonoff)

    tree.nboxes = nboxes[1]
    tree.nlev = nlev[1]

    if (ier[1] != 0)
        return tree, sorted_pts, ier[1]
    end

    ifixflag = [zero(Int32)]

    icolleagbox = tree.icolleagbox

    ccall( (:lrt2d_restrict_,LIBMBHFMM2D),
           Void, (Ref{Int32},Ref{Int32},Ref{Int32},
                  Ref{Int32},Ref{Int32},Ref{Int32},
                  Ref{Int32},Ref{Int32},Ref{Int32},
                  Ref{Int32},Ref{Int32},Ref{Int32}),
           levelbox,iparentbox,ichildbox,icolbox,
           irowbox,icolleagbox,nboxes,nlev,nblevel,
           iboxlev,istartlev,ifixflag)

    if (ifixflag[1] != 0)
        iflag = Array{Int32}(maxboxes)
        @printf "enforcing level restriction ...\n"
        ccall( (:lrt2d_fix_,LIBMBHFMM2D),
               Void, (Ref{Int32},Ref{Int32},Ref{Int32},
                      Ref{Int32},Ref{Int32},Ref{Int32},
                      Ref{Int32},Ref{Int32},Ref{Int32},
                      Ref{Int32},Ref{Int32},Ref{Int32},
                      Ref{Int32},Ref{Int32}),
               levelbox,iparentbox,ichildbox,icolbox,
               irowbox,icolleagbox,nboxes,nlev,nblevel,
               iboxlev,istartlev,iflag,maxboxes,
               itemparray)
    end

    tree.nboxes = nboxes[1]
    tree.nlev = nlev[1]

    @printf "testing tree ...\n"

    ccall( (:lrt2d_testtree_,LIBMBHFMM2D),
           Void, (Ref{Int32},Ref{Int32},Ref{Int32},
                  Ref{Int32},Ref{Int32},Ref{Int32},
                  Ref{Int32},Ref{Int32},Ref{Int32},
                  Ref{Int32}),
           levelbox,iparentbox,ichildbox,icolbox,
           irowbox,nboxes,nlev,nblevel,
           iboxlev,istartlev)
    
    ccall( (:lrt2d_mkcolls_,LIBMBHFMM2D),
           Void, (Ref{Int32},Ref{Int32},Ref{Int32},
                  Ref{Int32},Ref{Int32},Ref{Int32},
                  Ref{Int32},Ref{Int32},Ref{Int32},
                  Ref{Int32}),
           icolbox,irowbox,icolleagbox,nboxes,
           nlev,iparentbox,ichildbox,nblevel,
           iboxlev,istartlev)
    
    neighbors = tree.neighbors
    nnbrs = tree.nnbrs
    
    ccall( (:lrt2d_mknbrs_,LIBMBHFMM2D),
           Void, (Ref{Int32},Ref{Int32},Ref{Int32},
                  Ref{Int32},Ref{Int32},Ref{Int32},
                  Ref{Int32},Ref{Int32}),
           neighbors,nnbrs,nboxes,ichildbox,
           iparentbox,icolleagbox,icolbox,irowbox)

@printf "sorting points into tree ...\n"

ccall( (:lrt2d_ptsort_,LIBMBHFMM2D),
       Void, (Ref{Int32},Ref{Int32},Ref{Int32},
              Ref{Int32},Ref{Int32},Ref{Int32},
              Ref{Int32},Ref{Int32},Ref{Int32},              
              Ref{Int32},Ref{Float64},Ref{Float64},
              Ref{Int32},Ref{Int32},Ref{Int32},
              Ref{Float64},Ref{Float64},Ref{Int32}),
       levelbox,icolbox,irowbox,nboxes,nlev,
       iparentbox,ichildbox,nblevel,iboxlev,istartlev,
       src,srcsort,isrcsort,isrcladder,ns,
       zll,blength,ier)

ccall( (:lrt2d_ptsort_wc_,LIBMBHFMM2D),
       Void, (Ref{Int32},Ref{Int32},Ref{Int32},
              Ref{Int32},Ref{Int32},Ref{Int32},
              Ref{Int32},Ref{Int32},Ref{Int32},              
              Ref{Int32},Ref{Float64},Ref{Float64},
              Ref{Int32},Ref{Int32},Ref{Int32},
              Ref{Float64},Ref{Float64},Ref{Int32}),
       levelbox,icolbox,irowbox,nboxes,nlev,
       iparentbox,ichildbox,nblevel,iboxlev,istartlev,
       targ,targsort,itargsort,itargladder,nt,
       zll,blength,ier)



for i = 1:nboxes[1]
    localonoff[i] = 0
    if (itargladder[2,i]-itargladder[1,i]+1 > 0 || ifalltarg)
        localonoff[i] = 1
    end
end

return tree, sorted_pts, ier[1]

end
