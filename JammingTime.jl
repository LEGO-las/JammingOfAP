#returns jamming time with estimated error

len=200 #length
nc=400 #number of particles
k=2 #capacity of cell

a=1. #rate forward
b=0. #rate bacward
c=0.001 #switching rate

#surveying run:
T_it=50 #time of one iteration
imax=100 #number of iterations
n_t=30 #number of independent trajectories

ntd=200 #number of independent trajectories in the main run


function evolutionIII(gridC,gridK,positions,inside,directions,a,b,c,Tmax,k_max,Tin)
    #i=1
    len=length(gridC)
    nc=length(positions)
    pr_p=a/(a+b+c/2)
    pr_t=(a+b)/(a+b+c/2)
    tau=1/(nc*(2*a+2*b+c))
    
    t=Tin
    while t<Tmax
        t=t-tau*log(rand())
        particle=rand(1:nc)
        r=rand()
        if r>pr_t
            directions[particle]=~directions[particle]
        else

            if (r<pr_p)⊻(directions[particle])
                s=1
            else
                s=-1
            end
            
            if inside[particle]
               if (s>0)&&(rand()>0.5) 
                    pt=positions[particle]
                    if gridC[pt]<k_max-0.5
                        gridC[pt]=gridC[pt]+1
                        gridK[pt]=gridK[pt]-1
                        inside[particle]=false
                    end
                end 
            else
                if rand()<0.5
                    pt=positions[particle]
                    pn=mod1(pt+s,len)
                    if gridC[pn]<k_max-0.5
                        gridC[pn]=gridC[pn]+1
                        gridC[pt]=gridC[pt]-1
                        positions[particle]=pn
                    end
                else 
                    if s<0
                        pt=positions[particle]
                        if gridK[pt]<k_max-0.5
                            gridK[pt]=gridK[pt]+1
                            gridC[pt]=gridC[pt]-1
                            inside[particle]=true
                        end
                    end
                end
                        
                    
                
            end
            
            
        end
    end
    return t-Tmax
end


function initMS(len,nc)
    gridC=zeros(Int8,len)
    positions=zeros(Int64,nc)
    directions=falses(nc)
    gridK=zeros(Int8,len)
    inside=falses(nc)
    kam=1
    kapsa=false
    for i=1:nc
        CarryOn=true
        while CarryOn
            kam=rand(1:len)
            kapsa=rand(Bool)
            if kapsa
                if gridK[kam]<2
                    CarryOn=false
                end
            elseif gridC[kam]<2
                CarryOn=false
            end
        end
        inside[i]=kapsa
        positions[i]=kam
        directions[i]=kapsa
        if kapsa
            gridK[kam]+=1
        else
            gridC[kam]+=1
        end
    end
    
    return gridC,gridK,positions,inside,directions
end
  

function Current(gridC,positions,inside,directions,k,a,b)
    nc=length(positions)
    len=length(gridC)
    I=0
    for i=1:nc
        if ~inside[i]
            if ~directions[i]
                if gridC[mod1(positions[i]+1,len)]<k
                    I+=a
                end
                if gridC[mod1(positions[i]-1,len)]<k
                    I-=b
                end
            else
                if gridC[mod1(positions[i]-1,len)]<k
                    I-=a
                end
                if gridC[mod1(positions[i]+1,len)]<k
                    I+=b
                end
            end
        end
    end
    
    return I
end

function CurrentTraw(len,nc,a,b,c,k,T_it,imax,n_t)
    
    Is=zeros(n_t,imax+1)
    IT=zeros(imax+1)
    dIT=zeros(imax+1)
    

    for i_t=1:n_t
        gridCn,gridK,positions,inside,directions=initMS(len,nc)
        Tin=0
        Is[i_t,1]=Current(gridCn,positions,inside,directions,k,a,b)
        for i_s=2:imax+1
            Tin=evolutionIII(gridCn,gridK,positions,inside,directions,a,b,c,T_it,k,Tin)
            Is[i_t,i_s]=Current(gridCn,positions,inside,directions,k,a,b)
        end
    end
    
    
    return Is
end

ITss=CurrentTraw(len,nc,a,b,c,k,T_it,imax,n_t)


using Statistics

function pdfhist(list,from,to,n_bins)
    nd=length(list)
    pdf=zeros(n_bins)
    width=(to-from)/n_bins
    for i=1:nd
        t=list[i]
        if (t>=from) && (t<to)
            j=trunc(Int, (t-from)/width + 1)
            pdf[j]+=1
        end
    end
    
    return pdf ./ (nd*width), range(start=from+width/2,stop=to-width/2,step=width)
end

using LsqFit
function jamingtime(len,nc,a,b,c,T_it,imax,n_t,dI,nbin,ntd) #automatic version
    
    
    k=2

    #investigation run
    ITss=CurrentTraw(len,nc,a,b,c,k,T_it,imax,n_t)
    pI,XI=pdfhist(ITss[:],0,nc,nc÷dI)

    Id=0
    kolko=nc÷dI
    for i=3:kolko-2
        pt=pI[i]
        if  (pt>pI[i+1]) && (pt>pI[i-1]) && (0<pI[i+1]) && (0<pI[i-1]) && (pt>pI[i+2]) && (pt>pI[i-2])
            Id=XI[i] /2

            break
        end
    end

    #time estimation
    ntp=0.5
    for i=1:n_t
        if ITss[i,imax]<Id
            ntp+=1
        end
    end
    Te=T_it*imax*n_t/ntp*2

    #main run
    
    dT=Te/nbin - 1E-13
    Ts=range(start=dT,stop=2*Te-dT,step=2*dT)
    pdftj=zeros(nbin)

    nt=0
    while nt<ntd
        gridC,gridK,positions,inside,directions=initMS(len,nc)
        Tin=evolutionIII(gridC,gridK,positions,inside,directions,a,b,c,dT,k,0)
        i=0
        while (i<nbin+1) && (Current(gridC,positions,inside,directions,k,a,b)>Id)
            i+=1
            if i < nbin+1
            Tin=evolutionIII(gridC,gridK,positions,inside,directions,a,b,c,2*dT,k,Tin)
            end
        end
        if (i>0) && (i<nbin+1)
            pdftj[i]+=1
            nt+=1
        end
    end

    pdftj=pdftj ./(2*dT*ntd)

    itm=nbin
    for i=1:nbin
        if pdftj[i]==0
            itm=i-2
            break
        end
    end
    
    #fit
    @. model(x, p) = (p[1]) * exp(-x * p[1])

    p0=[1/dT]
    fit1 = LsqFit.curve_fit(model, Ts[1:itm],pdftj[1:itm], p0)
    cfde=coef(fit1)[1]
    odde=stderror(fit1)[1]

  
    return 1/cfde, odde/cfde/cfde
   
end


nbin=10
dI=5


jamingtime(len,nc,a,b,c,T_it,imax,n_t,dI,nbin,ntd)

