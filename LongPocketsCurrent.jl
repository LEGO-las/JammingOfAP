#returns steady current in long pockets geometry with estimated error

len=200 #length
nc=200 #number of particles
k=2 #capacity of cell
L1=20 #length of pockets

a=1.3 #rate forward
b=0.7 #rate bacward
c=0.01 #switching rate

T_eq=1000 #time of equlibration
T_it=100 #time between two measurements
iMax=200 #number of measurments during single experiment
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



function Current(grid,positions,inside,directions,k,a,b)
    nc=length(positions)
    len=length(grid)
    I=0
    for i=1:nc
        if ~inside[i]
            if ~directions[i]
                if grid[mod1(positions[i]+1,len)]<k
                    I+=a
                end
                if grid[mod1(positions[i]-1,len)]<k
                    I-=b
                end
            else
                if grid[mod1(positions[i]-1,len)]<k
                    I-=a
                end
                if grid[mod1(positions[i]+1,len)]<k
                    I+=b
                end
            end
        end
    end
    
    return I
end

function MeanCurrent(len,nc,a,b,c,k,T_eq,T_it,imax,n_t,L1)
    
    Is=zeros(n_t,imax)

    for i_t=1:n_t
        gridn,gridK,positions,inside,directions=initRK(len,nc,k)
        Tin=evolutionLongP(gridn,gridK,positions,inside,directions,a,b,c,T_eq,k,0,L1)
        for i_s=1:imax
            Tin=evolutionLongP(gridn,gridK,positions,inside,directions,a,b,c,T_it,k,Tin,L1)
            Is[i_t,i_s]=Current(gridn,positions,inside,directions,k,a,b)
        end
    end

    
    return mean(Is),std(Is)/sqrt(n_t*imax)
end



using Statistics

MeanCurrent(len,nc,a,b,c,k,T_eq,T_it,iMax,n_t,L1) ./len