


function evolutionI(grid,gridK,positions,inside,directions,loops,a,b,c,Tmax,k_max,Tin)
    #i=1
    len=length(grid)
    nc=length(positions)
    pr_p=a/(a+b+c)
    pr_t=(a+b)/(a+b+c)
    tau=1/(nc*(a+b+c))
    
    
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
            pt=positions[particle]
            pn=mod1(pt+s,len)
            if grid[pn]!=k_max
                grid[pn]=grid[pn]+1
                grid[pt]=grid[pt]-1
                positions[particle]=pn
                
                
                if (pt==1)&&(pn==len)
                    loops[particle]-=1
                elseif (pn==1)&&(pt==len)
                    loops[particle]+=1
                end
                
            end
        end
    end
    return t-Tmax
end

function evolutionII(grid,gridK,positions,inside,directions,loops,a,b,c,Tmax,k_max,Tin)
    #i=1
    len=length(grid)
    nc=length(positions)
    pr_p=a/(a+b+c)
    pr_t=(a+b)/(a+b+c)
    tau=1/(nc*(a+b+c))
    
    
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
            pt=positions[particle]
            pn=mod1(pt+s,len)
            if grid[pn]!=k_max*2
                grid[pn]=grid[pn]+1
                grid[pt]=grid[pt]-1
                positions[particle]=pn
                
                
                if (pt==1)&&(pn==len)
                    loops[particle]-=1
                elseif (pn==1)&&(pt==len)
                    loops[particle]+=1
                end
                
            end
        end
    end
    return t-Tmax
end

function evolutionIII(grid,gridK,positions,inside,directions,loops,a,b,c,Tmax,k_max,Tin)
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
                    if grid[pt]!=k_max
                        grid[pt]=grid[pt]+1
                        gridK[pt]=gridK[pt]-1
                        inside[particle]=false
                    end
                end 
            else
                if rand()<0.5
                    pt=positions[particle]
                    pn=mod1(pt+s,len)
                    if grid[pn]!=k_max
                        grid[pn]=grid[pn]+1
                        grid[pt]=grid[pt]-1
                        positions[particle]=pn
                        
                        if (pt==1)&&(pn==len)
                            loops[particle]-=1
                        elseif (pn==1)&&(pt==len)
                            loops[particle]+=1
                        end
                    end
                else 
                    if s<0
                        pt=positions[particle]
                        if gridK[pt]!=k_max
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


function evolutionIV(grid,gridK,positions,inside,directions,loops,a,b,c,Tmax,k_max,Tin)
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
            pt=positions[particle]
            
            if inside[particle]
               if ( (s>0)⊻(mod(pt,2)==1 ) )&&(rand()>0.5) 
                    
                    if grid[pt]<k_max-0.5
                        grid[pt]=grid[pt]+1
                        gridK[pt]=gridK[pt]-1
                        inside[particle]=false
                    end
                end 
            else
                if rand()<0.5 
                    
                    pn=mod1(pt+s,len)
                    if grid[pn]<k_max-0.5
                        grid[pn]=grid[pn]+1
                        grid[pt]=grid[pt]-1
                        positions[particle]=pn
                        
                        
                        if (pt==1)&&(pn==len)
                            loops[particle]-=1
                        elseif (pn==1)&&(pt==len)
                            loops[particle]+=1
                        end
                    end
                else 
                    if (s<0)⊻(mod(pt,2)==1 ) 
                        
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


function evolutionV(grid,gridK,positions,inside,directions,loops,a,b,c,Tmax,k_max,Tin)
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
                    pn=mod1(pt+s,len)
                    if gridK[pn]!=k_max
                        gridK[pn]=gridK[pn]+1
                        gridK[pt]=gridK[pt]-1
                        positions[particle]=pn
                        
                        if (pt==1)&&(pn==len)
                            loops[particle]-=1
                        elseif (pn==1)&&(pt==len)
                            loops[particle]+=1
                        end
                    end
                elseif s>0
                    pt=positions[particle]
                    if grid[pt]!=k_max
                        grid[pt]=grid[pt]+1
                        gridK[pt]=gridK[pt]-1
                        inside[particle]=false
                    end
                end 
            else
                if rand()<0.5 
                    pt=positions[particle]
                    pn=mod1(pt+s,len)
                    if grid[pn]!=k_max
                        grid[pn]=grid[pn]+1
                        grid[pt]=grid[pt]-1
                        positions[particle]=pn
                        
                        if (pt==1)&&(pn==len)
                            loops[particle]-=1
                        elseif (pn==1)&&(pt==len)
                            loops[particle]+=1
                        end
                    end
                else 
                    if s<0
                        pt=positions[particle]
                        if gridK[pt]!=k_max
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

function initRD(len,nc)
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
            if grid[place]<4
                CarryOn=false
            end
        end
        positions[i]=place
        directions[i]=rand(Bool)
        grid[place]+=1
    end
    
    return grid,gridK,positions,inside,directions
end

function initR(len,nc)
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
            if grid[place]<2
                CarryOn=false
            end
        end
        positions[i]=place
        directions[i]=rand(Bool)
        grid[place]+=1
    end
    
    return grid,gridK,positions,inside,directions
end


function d2OneTrajectory(len,nc,a,b,c,T_eq,T_it,imax,k_max,fevolution)


    positions_zaznam=zeros(nc,imax)
    loops=zeros(nc)
    if fevolution==evolutionI
        grid,gridK,positions,inside,directions=initR(len,nc)
    elseif fevolution==evolutionII
        grid,gridK,positions,inside,directions=initRD(len,nc)
    else
        grid,gridK,positions,inside,directions=initRK(len,nc)
    end
    Tin=fevolution(grid,gridK,positions,inside,directions,loops,a,b,c,T_eq,k_max,0)
    positions_zaznam[:,1]=positions[:]
    positions_zaznam[:,1]+=len*loops

    for i=2:imax
        Tin=Tin=fevolution(grid,gridK,positions,inside,directions,loops,a,b,c,T_it,k_max,Tin)
        positions_zaznam[:,i]=positions[:]
        positions_zaznam[:,i]+=len*loops
    end
    for j=1:nc
        positions_zaznam[j,:]=positions_zaznam[j,:] .- positions_zaznam[j,1]
    end


    positions_s=zeros(imax)
    for j=1:nc
        positions_s[:]=positions_s[:] .+ positions_zaznam[j,:]
    end
    positions_s[:]=positions_s[:]/nc

    d2=zeros(imax)
    for j=1:nc
        d2[:]=d2[:] .+ (positions_zaznam[j,:] .- positions_s[:]  ) .^2
    end
    d2[:]=d2[:]/(nc-1)
    
    return d2
end



using CurveFit

function DiffusionCoefficient(len,nc,a,b,c,k,T_eq,T_it0,rel_eps,n_t,fevolution)
    
    T_it=T_it0
    imax=80 


    d2s=zeros(n_t,imax)

    reOK=false
    f1OK=false
    f2OK=false

    Mb,Db=0.,0.

    while !(reOK & f1OK & f2OK) 

        for i_t=1:n_t
            d2s[i_t,:]=d2OneTrajectory(len,nc,a,b,c,T_eq,T_it,imax,k,fevolution)
        end

        M,D,X=process(d2s,true,T_it)

        a1,b1=linear_fit(X[40:80], M[40:80])

        a2,b2=linear_fit(X[60:80], M[60:80])

        Mb=(b1+b2)/2
        Db=abs(b1-b2)

        if Db/Mb<rel_eps
            reOK=true
        else
            reOK=false
            #println("eps")
        end

        Err_fit1=0
        Err_fit2=0

        Err_dat=0

        for i=40:80
            Err_dat+=D[i]^2
            Err_fit2+=(a2+b2*X[i]-M[i] )^2
            Err_fit1+=(a1+b1*X[i]-M[i] )^2
        end

        if Err_fit1<Err_dat
            f1OK=true
        else
            f1OK=false
            #println("1")
        end

        if Err_fit2<Err_dat
            f2OK=true
        else
            f2OK=false
            #println("2")
        end

        T_it*=2
    end
    
    return Mb,Db  #,T_it/2
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


len=200
nc=10
a=1.
b=0.
c=0.001
T_eq=3000.
T_it0=20.

fevolution=evolutionV  #pick geometry: evolutionI, evolutionII, evolutionIII, evolutionIV, evolutionV


n_t=500

rel_eps=0.05 #required relative precision
k=2 #capacity


DiffusionCoefficient(len,nc,a,b,c,k,T_eq,T_it0,rel_eps,n_t,fevolution)


