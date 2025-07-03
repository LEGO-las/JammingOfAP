#returns unjammed current with estimated error, then jammed current with error

len=200 #length
nc=300 #number of particles
k=2 #capacity of cell

a=1. #rate forward
b=0. #rate bacward
c=0.003 #switching rate

T_eq=1000 #time of equlibration
T_e=3000 #time of experiment
iMax=200 #number of measurments during single experiment
n_t=200 #number of independent experiments

using Statistics

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


function MeanCurrent(len,nc,a,b,c,k,T_e,T_eq,iMax,n_t,Itr)
    
    Ineq=zeros(n_t,iMax)
    gridC,gridK,positions,inside,directions=initMS(len,nc)
    Tin=0

    for i_t=1:n_t
        
        gridC,gridK,positions,inside,directions=initMS(len,nc)
        Tin=evolutionIII(gridC,gridK,positions,inside,directions,a,b,c,T_eq,k,0)
        for i_S=1:iMax
            Tin=evolutionIII(gridC,gridK,positions,inside,directions,a,b,c,T_e/iMax,k,Tin)
            Ineq[i_t,i_S]=Current(gridC,positions,inside,directions,k,a,b)
        end
    end
    
    IMS=zeros(Union{Missing, Float64},n_t,iMax)
    IMS[:,:]=Ineq[:,:]
    NMS=0
    
    IJ=zeros(Union{Missing, Float64},n_t,iMax)
    IJ[:,:]=Ineq[:,:]
    NJ=0
    for t=1:n_t
        for i=1:iMax
            if Ineq[t,i]<Itr
                IMS[t,i]=missing
                NJ+=1
            else
                IJ[t,i]=missing
                NMS+=1
            end
        end
    end
    iMS=mean(skipmissing(IMS))
    diMS=std(skipmissing(IMS))
    
    
    iJ=mean(skipmissing(IJ))
    diJ=std(skipmissing(IJ))
    
    
    return iMS,diMS/sqrt(NMS),iJ,diJ/sqrt(NJ)
end



Itr=nc*(4*len-nc)/(16*len)*(a-b)^2 /(a+b)


MeanCurrent(len,nc,a,b,c,k,T_e,T_eq,iMax,n_t,Itr) ./len

