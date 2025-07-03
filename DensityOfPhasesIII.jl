#returns densities of of "gas" phase with errors, then densities of "condensate" phase with errors
#in Geometry III

len=400 #length of the system
nc=800 #number of particles

Teq=1000 #equlibration time
dT=300 #time between measurements
nTs=2000 #number of measurements

cs=[0.003,0.006,0.01,0.04] #array of switching rates
as=[1.001,1.1,1.2,1.4,1.7] #array of forward jump rates

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
        este=true
        while este
            place=rand(1:len)
            kapsa=rand(Bool)
            if kapsa
                if gridK[place]<2
                    este=false
                end
            elseif grid[place]<2
                este=false
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


function rhoPhasesMaxInt(prho,Dr,window,k,int1,int2)
    
    rhosA=range(start=0,stop=k,step=1/window/2)
    
    max1=argmax(prho[int1[1]:int1[2]])+int1[1]-1
    max2=argmax(prho[int2[1]:int2[2]])+int2[1]-1
    if max1==int1[2]
        while prho[max1]<prho[max1+1]
            max1+=1
        end
    end
    if max2==int2[1]
        while prho[max2]<prho[max2-1]
            max2-=1
        end
    end
    
    rhoout1=rhosA[max1]
    rhoout2=rhosA[max2]
    
    int1=[maximum([max1-1,1]), max1+3]
    int2=[max2-3 ,minimum([max2+1,1+2*k*window])]

    return rhoout1, rhoout2 , int1,int2
end

function CGdensityOver(grid,gridK,window,k)
    rhoA=zeros(k*window*2+1)
    len=length(grid)
    for i=1:len-window
        nic=sum(grid[i:window+i-1])
        nik=sum(gridK[i:window+i-1])
        rhoA[nic+nik+1]+=1
    end
    return rhoA
end

function RhosPhasesRunInt(len,a,b,c,k,windows,nc,Dr,Teq,dT,nTs,ints1,ints2)
    nos=length(windows)
    rhoA = [Int[] for _ in 1:nos]
    for i=1:nos
        rhoA[i]=zeros(k*windows[i]*2+1)
    end


    grid,gridK,positions,inside,directions=initRK(len,nc)
    Tin=evolutionIII(grid,gridK,positions,inside,directions,a,b,c,Teq,k,0)


    for i=1:nTs
        Tin=evolutionIII(grid,gridK,positions,inside,directions,a,b,c,dT,k,Tin)
        for j=1:nos
            rhoAtemp=CGdensityOver(grid,gridK,windows[j],k)
            rhoA[j]+=rhoAtemp
        end
    end
    rhoVs = [Float64[] for _ in 1:nos]
    for i=1:nos
        rhoVs[i]=zeros(2)
    end

    for i=1:nos
        rtemp1,rtemp2,int1temp,int2temp=rhoPhasesMaxInt(rhoA[i],Dr,windows[i],k,ints1[i],ints2[i])
        rhoVs[i][1]=rtemp1
        rhoVs[i][2]=rtemp2
        
        ints1[i]=int1temp
        ints2[i]=int2temp
    end
    
    vys=(rhoVs[2] .+rhoVs[1])/2
    errs=(rhoVs[2] .-rhoVs[1])/2
    #println(rhoVs)
    
    
    return vys[1], abs(errs[1]), vys[2], abs(errs[2])
end

using StatsBase


k=2

windows=[10,20]


ncs=length(cs)
nas=length(as)

mr1=zeros(ncs,nas)
er1=zeros(ncs,nas)

mr2=zeros(ncs,nas)
er2=zeros(ncs,nas)

for ia=1:nas
    
    ints1=[[1,20],[1,40] ]

    ints2=[[21,41],[41,81] ]
    for ic=1:ncs
        mr1[ic,ia],er1[ic,ia],mr2[ic,ia],er2[ic,ia]= RhosPhasesRunInt(len,as[ia],as[ia]-1.,cs[ic],k,windows,nc,Dr,Teq,dT,nTs,ints1,ints2)
    end
end

mr1,er1,mr2,er2