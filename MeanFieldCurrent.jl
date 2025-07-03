#returns meanfield estimation of current

len=200 #length
nc=300 #number of particles

a=1. #rate forward
b=0. #rate bacward
c=0.003 #switching rate


function init_guess2III(ni,rho)
    
    ps=zeros(3)
    ps[3]=rho/2 + (1-sqrt(1+2*rho-rho^2))/2
    ps[2]=rho - 2*ps[3]
    ps[1]=1 + ps[3] - rho
    
    
    n=length(ni)
    
    ns=zeros(3)
    for i=1:n
        Ri=Int(ni[i])
        if mod(Ri,2)==0
            ns[div(Ri,2)+1]+=1
        end
    end
    #println(ns)
    ps= ps ./ ns
    init_guess=zeros(n)
    for i=1:n
        Ri=Int(ni[i])
        if mod(Ri,2)==0
            init_guess[i]=ps[div(Ri,2)+1]
        end
    end
    return init_guess
end
  


function ni2III()
    n=36
    ni=zeros(n)
    # 1 - Pocket
    for p=1:6
        s1=2+(p-1)*6
        s2=3+(p-1)*6
        ni[s1]+=1
        ni[s2]+=1
    end
    # 1 - BackBone
    for k=1:6
        s1=6+k
        s2=12+k
        ni[s1]+=1
        ni[s2]+=1
    end
    # 2 - Pocket
    for p=1:6
        s1=4+(p-1)*6
        s2=5+(p-1)*6
        s3=6+(p-1)*6
        ni[s1]+=2
        ni[s2]+=2
        ni[s3]+=2
    end
    # 2 - BackBone
    for k=1:6
        s1=18+k
        s2=24+k
        s3=30+k
        ni[s1]+=2
        ni[s2]+=2
        ni[s3]+=2
    end
    return ni
end

function TransitionMatrix2III(state,dt,a,c)
    n=36
    TransitionMatrix=zeros(n,n)

    #SWITCHING DIRECTION
    #switching between 1 - Pocket
    for p=1:6
        s1=2+(p-1)*6
        s2=3+(p-1)*6
        TransitionMatrix[s1,s2]=c*dt
        TransitionMatrix[s2,s1]=c*dt
    end
    #switching between 1 - BackBone
    for k=1:6
        s1=6+k
        s2=12+k
        TransitionMatrix[s1,s2]=c*dt
        TransitionMatrix[s2,s1]=c*dt
    end
    #switching between 2 - Pocket
    for p=1:6
        s1=4+(p-1)*6
        s2=5+(p-1)*6
        s3=6+(p-1)*6
        TransitionMatrix[s2,s1]=2*c*dt
        TransitionMatrix[s2,s3]=2*c*dt
        TransitionMatrix[s1,s2]=c*dt
        TransitionMatrix[s3,s2]=c*dt
    end
    #switching between 2 - BackBone
    for k=1:6
        s1=18+k
        s2=24+k
        s3=30+k
        TransitionMatrix[s2,s1]=2*c*dt
        TransitionMatrix[s2,s3]=2*c*dt
        TransitionMatrix[s1,s2]=c*dt
        TransitionMatrix[s3,s2]=c*dt
    end

    #JUMPS
    #01-10
    TransitionMatrix[2,7]=a*dt
    TransitionMatrix[13,3]=a*dt
    #20-11
    TransitionMatrix[8,19]=2*a*dt
    TransitionMatrix[14,25]=a*dt
    TransitionMatrix[25,9]=a*dt
    TransitionMatrix[31,15]=a*dt
    #11-20
    TransitionMatrix[15,6]=2*a*dt
    TransitionMatrix[14,5]=a*dt
    TransitionMatrix[5,9]=a*dt
    TransitionMatrix[4,8]=a*dt
    #21-12
    TransitionMatrix[10,20]=2*a*dt
    TransitionMatrix[11,21]=2*a*dt
    TransitionMatrix[26,11]=a*dt
    TransitionMatrix[16,26]=a*dt
    TransitionMatrix[27,12]=2*a*dt
    TransitionMatrix[17,27]=a*dt
    TransitionMatrix[32,17]=a*dt
    TransitionMatrix[33,18]=2*a*dt


    #JUMPING AWAY
    p_v=sum(state[1:18])
    #1->0
    for k=1:6
        s0=k
        s11=6+k
        s12=12+k
        TransitionMatrix[s0,s11]=p_v*a*dt
        TransitionMatrix[s0,s12]=p_v*a*dt
    end
    #2->1
    for k=1:6
        s11=6+k
        s12=12+k
        s21=18+k
        s22=24+k
        s23=30+k
        TransitionMatrix[s11,s21]=p_v*a*dt*2
        TransitionMatrix[s12,s23]=p_v*a*dt*2
        TransitionMatrix[s11,s22]=p_v*a*dt
        TransitionMatrix[s12,s22]=p_v*a*dt
    end

    #JUMPING IN
    p_p=sum(state[7:12])+2*sum(state[19:24])+sum(state[25:30])
    p_l=sum(state[13:18])+2*sum(state[31:36])+sum(state[25:30])
    #0->1
    for k=1:6
        s0=k
        s11=6+k
        s12=12+k
        TransitionMatrix[s11,s0]=p_p*a*dt
        TransitionMatrix[s12,s0]=p_l*a*dt
    end
    #1->2
    for k=1:6
        s11=6+k
        s12=12+k
        s21=18+k
        s22=24+k
        s23=30+k
        TransitionMatrix[s21,s11]=p_p*a*dt
        TransitionMatrix[s23,s12]=p_l*a*dt
        TransitionMatrix[s22,s11]=p_l*a*dt
        TransitionMatrix[s22,s12]=p_p*a*dt
    end

    for i=1:n
        TransitionMatrix[i,i]=0
        TransitionMatrix[i,i]=1-sum(TransitionMatrix[:,i])
    end

    return TransitionMatrix
end



Nit=10000
ni=ni2III()
#SS=SSBPIII()
dt=0.01


state=init_guess2III(ni,nc/(2*len))
for _=1:Nit

    state= TransitionMatrix2III(state,dt,a,c) * state
end

p_p=sum(state[7:12])+2*sum(state[19:24])
p_l=sum(state[13:18])+2*sum(state[31:36])
p_v=sum(state[1:18])


a*p_v*(p_l-p_p)