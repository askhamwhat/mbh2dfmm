

include("MBHFMM2D.jl")

using PyPlot

## parameters

lambda = 1.0e-2

ifcharge = false
ifdipole = false
ifquad = false
ifoct = true
ifpot = true
ifgrad = true
ifhess = true

## define source and target points

ns = 10000
nt = 2000

# sources on an ellipse, targets inside

src = Array(Float64,2,ns)
rnorm = Array(Float64,2,ns)
dsdt = Array(Float64,ns)
targ = Array(Float64,2,nt)

charge = zeros(Float64,1)
dipstr = zeros(Float64,1)
dipvec = zeros(Float64,2,1)
quadstr = zeros(Float64,1)
quadvec = zeros(Float64,3,1)
octstr = Array(Float64,ns) # octopole strength
octvec = Array(Float64,4,ns) # octopole "direction"

# storage for solution

pottarg1 = zeros(Float64,nt)
gradtarg1 = zeros(Float64,2,nt)
hesstarg1 = zeros(Float64,3,nt)
pottarg2 = zeros(Float64,nt)
gradtarg2 = zeros(Float64,2,nt)
hesstarg2 = zeros(Float64,3,nt)

center = [6.1;-1.3]

h = 2.0*pi/ns
a = 1.123
b = 0.7

for i = 1:ns
    t = (i-1)*h
    src[:,i] = center + [a*cos(t);b*sin(t)]
    dx = -a*sin(t)
    dy = b*cos(t)
    dsdt[i] = sqrt(dx^2+dy^2)
    rnorm[1,i] = dy/dsdt[i]
    rnorm[2,i] = -dx/dsdt[i]
    nu = rnorm[:,i]
    tau = [-nu[2];nu[1]]
    octvec[:,i] = boxfmm2d_formoctvec(nu,nu,tau)
    octstr[i] = h*dsdt[i]*cos(5*t)
end

for i = 1:nt
    xt = randn(); yt = randn();
    elnorm = sqrt((xt/a)^2 + (yt/b)^2)
    sc = rand();
    xt = sc*xt/elnorm; yt = sc*yt/elnorm
    targ[:,i] = center + [xt;yt]
end


## make parameters structure

fmmpars = MBHFMM2DParams(lambda,src,targ,ifcharge, ifdipole,
                         ifquad, ifoct, charge, dipstr,
                         dipvec, quadstr, quadvec, octstr,
                         octvec, iprec = 3, ifalltarg = false)

## direct calculation

@printf "PERFORMING DIRECT CALCULATION ...\n"

@time mbhfmm2d_direct!(fmmpars,targ,ifpot,pottarg1,ifgrad,
                 gradtarg1,ifhess,hesstarg1)

## set up fmm

# maxboxes determines size of storage arrays, set it big
# for now
# maxnodes determines max number of sources/targets
# per box

@printf "FORMING FMM ...\n"

@time fmmstor = mbhfmm2d_form(fmmpars,maxnodes=30,maxboxes=100000)

@printf "EVALUATING FMM ...\n"

@time mbhfmm2d_targ!(fmmpars,fmmstor,targ,ifpot,pottarg2,ifgrad,
                 gradtarg2,ifhess,hesstarg2)

println("relative error, pot ",
        vecnorm(pottarg1-pottarg2)/vecnorm(pottarg1))
println("relative error, grad ",
        vecnorm(gradtarg1-gradtarg2)/vecnorm(gradtarg1))
println("relative error, hess ",
        vecnorm(hesstarg1-hesstarg2)/vecnorm(hesstarg1))
