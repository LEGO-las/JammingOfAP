#returns array of mean numbers of particles propelled in easy direction (with error), difficult direction (and error), free capacity (and error)
#at first for pocket, then for backbone


len=200 #length
nc=200 #number of particles
k=2 #capacity of cell
L1=20 #length of pockets

a=1.3 #rate forward
b=0.7 #rate bacward
c=0.01 #switching rate

T_eq=1000 #time of equlibration
T_it=100 #time between two measurements
imax=200 #number of measurments during single experiment
n_t=200 #number of independent experiments
function evolutionLongP(grid,gridK,positions,inside,directions,a,b,c,Tmax,k_max,Tin,L1)
    #i=1
    len=length(grid)
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
                if (rand()>0.5)
                    pt=positions[particle]
                    if (s>0)
                    
                        if mod1(pt,L1)==1 
                            if grid[pt]<k_max-0.5
                                grid[pt]=grid[pt]+1
                                gridK[pt]=gridK[pt]-1
                                inside[particle]=false
                            end
                        else
                            pn=mod1(pt+s,len)
                            if gridK[pn]<k_max-0.5
                                gridK[pn]=gridK[pn]+1
                                gridK[pt]=gridK[pt]-1
                                positions[particle]=pn
                            end
                        end
                    else 
                        pn=mod1(pt+s,len)
                        pnt=mod1(pt+s,L1)
                    
                        if pnt>1 
                            
                            if gridK[pn]<k_max-0.5
                                gridK[pn]=gridK[pn]+1
                                gridK[pt]=gridK[pt]-1
                                positions[particle]=pn
                            end
                        end
                    end
                end 
            else
                if rand()<0.5 
                    pt=positions[particle]
                    pn=mod1(pt+s,len)
                    if grid[pn]<k_max-0.5
                        grid[pn]=grid[pn]+1
                        grid[pt]=grid[pt]-1
                        positions[particle]=pn
                    end
                else 
                    pt=positions[particle]
                    if (s<0)&&(mod1(pt,L1)==1)
                        if gridK[pt]<k_max-0.5
                            gridK[pt]=gridK[pt]+1
                            grid[pt]=grid[pt]-1
                            inside[particle]=true
                           
                        end
                    end
                end
            end
            
            
        end
    end
    return t-Tmax
end

using Random


function initRK(len,nc,k)
    grid=zeros(Int8,len)
    positions=zeros(Int64,nc)
    directions=falses(nc)
    gridK=zeros(Int8,len)
    inside=falses(nc)
    kam=1
    kapsa=false
    for i=1:nc
        este=true
        while este
            kam=rand(1:len)
            kapsa=rand(Bool)
            if kapsa
                if gridK[kam]<k
                    este=false
                end
            elseif grid[kam]<k
                este=false
            end
        end
        inside[i]=kapsa
        positions[i]=kam
        directions[i]=rand(Bool)
        if kapsa
            gridK[kam]+=1
        else
            grid[kam]+=1
        end
    end
    
    return grid,gridK,positions,inside,directions
end

function DirectionsInPockets(len,L1,positions,inside,directionsy)
    nkaps=len÷L1
    outward=zeros(L1)
    inward=zeros(L1)
    nothing=zeros(L1)
    nc=length(positions)
    for ci=1:nc
        if inside[ci]
            if directionsy[ci]
                inward[mod1(positions[ci]-1,L1)]+=1
            else
                outward[mod1(positions[ci]-1,L1)]+=1
            end
        end
    end
    outward= outward/ (nkaps)
    inward= inward/ (nkaps)
    nothing=2 .- outward .- inward
    return outward, inward, nothing
end

function DirectionsInBB(len,L1,positions,inside,directionsy)
    nkaps=len÷L1
    outward=zeros(L1)
    inward=zeros(L1)
    nothing=zeros(L1)
    nc=length(positions)
    for ci=1:nc
        if inside[ci]==false
            if directionsy[ci]
                inward[mod1(positions[ci]-1,L1)]+=1
            else
                outward[mod1(positions[ci]-1,L1)]+=1
            end
        end
    end
    outward= outward/ (nkaps)
    inward= inward/ (nkaps)
    nothing=2 .- outward .- inward
    return outward, inward, nothing
end

function DirectionsInPocketsRun(len,nc,a,b,c,k,T_eq,T_it,imax,n_t,L1)
    outwards=zeros(n_t*imax,L1)
    inwards=zeros(n_t*imax,L1)
    nothings=zeros(n_t*imax,L1)
    outwardsB=zeros(n_t*imax,L1)
    inwardsB=zeros(n_t*imax,L1)
    nothingsB=zeros(n_t*imax,L1)

    for i_t=1:n_t
        gridn,gridK,positions,inside,directionsy=initRK(len,nc,k)
        Tin=evolutionLongP(gridn,gridK,positions,inside,directionsy,a,b,c,T_eq,k,0,L1)
        for i_s=1:imax
            Tin=evolutionLongP(gridn,gridK,positions,inside,directionsy,a,b,c,T_it,k,Tin,L1)
            outwards[(i_t-1)*imax+i_s,:],inwards[(i_t-1)*imax+i_s,:],nothings[(i_t-1)*imax+i_s,:]=DirectionsInPockets(len,L1,positions,inside,directionsy)
            outwardsB[(i_t-1)*imax+i_s,:],inwardsB[(i_t-1)*imax+i_s,:],nothingsB[(i_t-1)*imax+i_s,:]=DirectionsInBB(len,L1,positions,inside,directionsy)
        end
    end
    
    Mv,Dv,X=process(outwards,true,1)
    Md,Dd,X=process(inwards,true,1)
    Mn,Dn,X=process(nothings,true,1)
    MvB,DvB,X=process(outwardsB,true,1)
    MdB,DdB,X=process(inwardsB,true,1)
    MnB,DnB,X=process(nothingsB,true,1)
    
    return Mv,Dv,Md,Dd,Mn,Dn,MvB,DvB,MdB,DdB,MnB,DnB,X
end



using Statistics
function process(data,DevOrErr,x_step )
    n_sample=length(data[:,1])
    n_x=length(data[1,:])
    AV=zeros(n_x)
    errorbar=zeros(n_x)
    for x=1:n_x
        AV[x]=mean(data[:,x])
        errorbar[x]=std(data[:,x])
    end
    osx=range(start=0,stop=x_step*(n_x-1),step=x_step)
    
    if DevOrErr
        errorbar/= sqrt(n_sample)
    end
    
    return AV,errorbar,osx
end

DirectionsInPocketsRun(len,nc,a,b,c,k,T_eq,T_it,imax,n_t,L1)