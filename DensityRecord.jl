#returns a record of evolution of local densties in the channel

len=200 #length
nc=400 #number of particles
k=2 #capacity of cell

a=1. #rate forward
b=0. #rate bacward
c=0.001 #switching rate

T_eq=5000 #time for equilibration
dT=50 #time of one iteration
Nmax=100 #number of iterations
n_t=30 #number of independent trajectories


function evolutionIII(grid,gridK,positions,inside,directions,a,b,c,Tmax,k_max,Tin)
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
               if (s>0)&&(rand()>0.5) 
                    pt=positions[particle]
                    if grid[pt]<k_max-0.5
                        grid[pt]=grid[pt]+1
                        gridK[pt]=gridK[pt]-1
                        inside[particle]=false
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
                    if s<0
                        pt=positions[particle]
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

function initRK(len,nc)
    grid=zeros(Int8,len)
    positions=zeros(Int64,nc)
    directions=falses(nc)
    gridK=zeros(Int8,len)
    inside=falses(nc)
    place=1
    kapsa=false
    for i=1:nc
        CarryOn=true
        while CarryOn
            place=rand(1:len)
            kapsa=rand(Bool)
            if kapsa
                if gridK[place]<2
                    CarryOn=false
                end
            elseif grid[place]<2
                CarryOn=false
            end
        end
        inside[i]=kapsa
        positions[i]=place
        directions[i]=rand(Bool)
        if kapsa
            gridK[place]+=1
        else
            grid[place]+=1
        end
    end
    
    return grid,gridK,positions,inside,directions
end


function pohybzber_vyvojgridu(grid,gridK,positions,inside,directions,a,b,c,dT,Nmax,k_max)
    i=1
    nc=length(positions)
    len=length(grid)
    Tin=0
    
    grids=zeros(len,Nmax)
    while i<=Nmax
        
        Tin=evolutionIII(grid,gridK,positions,inside,directions,a,b,c,dT,k_max,Tin)
        grids[:,i]=grid[:]
        i=i+1
        
        
    end
    
    return grids
end



gridP,gridK,positionsP,inside,directionsP=initRK(len,nc)

evolutionIII(gridP,gridK,positionsP,inside,directionsP,a,b,c,T_eq,k,0)

pohybzber_vyvojgridu(gridP,gridK,positionsP,inside,directionsP,a,b,c,dT,Nmax,k)