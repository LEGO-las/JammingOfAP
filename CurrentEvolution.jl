#returns n_t independent current evolutions

len=200 #length
nc=400 #number of particles
k=2 #capacity of cell

a=1. #rate forward
b=0. #rate bacward
c=0.001 #switching rate
T_it=50 #time of one iteration
imax=100 #number of iterations
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

            if xor((r<pr_p),(directions[particle]))
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


function initMS(len,nc)
    grid=zeros(Int8,len)
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
            elseif grid[kam]<2
                CarryOn=false
            end
        end
        inside[i]=kapsa
        positions[i]=kam
        directions[i]=kapsa
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

function CurrentTraw(len,nc,a,b,c,k,T_it,imax,n_t)
    
    Is=zeros(n_t,imax+1)
    IT=zeros(imax+1)
    dIT=zeros(imax+1)
    

    for i_t=1:n_t
        gridn,gridK,positions,inside,directions=initMS(len,nc)
        Tin=0
        Is[i_t,1]=Current(gridn,positions,inside,directions,k,a,b)
        for i_s=2:imax+1
            Tin=evolutionIII(gridn,gridK,positions,inside,directions,a,b,c,T_it,k,Tin)
            Is[i_t,i_s]=Current(gridn,positions,inside,directions,k,a,b)
        end
    end
    
    
    return Is
end

ITss=CurrentTraw(len,nc,a,b,c,k,T_it,imax,n_t) ./len
