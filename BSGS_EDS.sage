class Block:
    def __init__(self,value,subscript,mark):
        self.value=value
        self.subscript=subscript
        self.mark=mark

class Element:
    def __init__(self,value,subscript):
        self.value=value
        self.subscript=subscript

def W(Q,l):
    if(l==0):
        return 0
    if(l<0):
        sign=-1
        l=-l
    else:
        sign=1
    b4=2*A
    b6=4*B
    b8=-A^2
    h_0=0
    h_1=1
    h_2=2*Q[1]
    h_3=3*Q[0]^4+3*b4*Q[0]^2+3*b6*Q[0]+b8
    h_4=h_2*(2*Q[0]^6+5*b4*Q[0]^4+10*b6*Q[0]^3+10*b8*Q[0]^2-b4*b6*Q[0]+b4*b8-b6^2)
    h2ver=h_2^(-1)
    h_5=h_4*h_2**3-h_3**3
    h_6=h_3*(h_5*h_2**2-h_4**2)*h2ver
    h_7=h_5*h_3**3-h_2*h_4**3
    f=factor(l)
    def factor2(n):
        if (n%2==0):
            s=factor(n)[0][1]
            return(s,n/(2**s))
        else:
            return(0,n)
    s=factor2(l)[0]
    u_0=factor2(l)[1]
    u=[]
    u.append(u_0)
    for i in range(0,u_0):
        x=floor(u[i]/2)
        if (x==1):
            u.append(x)
            break
        else:
            u.append(x)
    u=u[::-1]
    def BL_db(BL):
        h_2ns3=BL[3]*BL[1]^3-BL[0]*BL[2]^3
        h_2ns2=BL[2]*(BL[4]*BL[1]**2-BL[0]*BL[3]**2)*h2ver
        h_2ns1=BL[4]*BL[2]^3-BL[1]*BL[3]^3
        h_2n=BL[3]*(BL[5]*BL[2]**2-BL[1]*BL[4]**2)*h2ver
        h_2na1=BL[5]*BL[3]^3-BL[2]*BL[4]^3
        h_2na2=BL[4]*(BL[6]*BL[3]**2-BL[2]*BL[5]**2)*h2ver
        h_2na3=BL[6]*BL[4]^3-BL[3]*BL[5]^3
        h_2na4=BL[5]*(BL[7]*BL[4]**2-BL[3]*BL[6]**2)*h2ver
        return [h_2ns3, h_2ns2, h_2ns1, h_2n, h_2na1, h_2na2, h_2na3, h_2na4]
    def BL_dba(BL):
        h_2ns2=BL[2]*(BL[4]*BL[1]**2-BL[0]*BL[3]**2)*h2ver
        h_2ns1=BL[4]*BL[2]^3-BL[1]*BL[3]^3
        h_2n=BL[3]*(BL[5]*BL[2]**2-BL[1]*BL[4]**2)*h2ver
        h_2na1=BL[5]*BL[3]^3-BL[2]*BL[4]^3
        h_2na2=BL[4]*(BL[6]*BL[3]**2-BL[2]*BL[5]**2)*h2ver
        h_2na3=BL[6]*BL[4]^3-BL[3]*BL[5]^3
        h_2na4=BL[5]*(BL[7]*BL[4]**2-BL[3]*BL[6]**2)*h2ver
        h_2na5=BL[7]*BL[5]^3-BL[4]*BL[6]^3
        return [h_2ns2, h_2ns1, h_2n, h_2na1, h_2na2, h_2na3, h_2na4, h_2na5]
    BL=[-h_2,-1,0,1,h_2,h_3,h_4,h_5]
    for i in range(1,len(u)):
        if(i==1):
            if u[i]==2:
                BL=[-1,0,1,h_2,h_3,h_4,h_5,h_6]
            elif u[i]==3:
                BL=[0,1,h_2,h_3,h_4,h_5,h_6,h_7]
        if(i==len(u)-1):
            break
        if(2*u[i]==u[i+1]):
            BL=BL_db(BL)
        else:
            BL=BL_dba(BL)
    for i in range(0,s):
        BL=BL_db(BL)
    return sign*BL[3]

def phi(Q):
    if(Q==0):
        return 0
    else:
        return (W(Q,p-1)/W(Q,p-1+d))^(x**2)/phiP

def phi2(n):
    if(n==0):
        return 0
    else:
        return phiP^(n^2-1)*W(P,n)

def phi3(Q):
    if(Q==0):
        return 0
    else:
        return (W(Q,p-1)/W(Q,p-1+d))^(x**2)

def ezhash(a,b):
    if a==0:
        return 0
    else:
        return lift(a)%b


def hashsearch(List,Num,m,SUBSCRIPT):
    hash_=ezhash(Num,m)
    arr=[]
    OBJ=[obj for obj in List[hash_] if obj.value == Num]
    for a in OBJ:
        arr.append((SUBSCRIPT+a.subscript*m)%d)
    return arr

def Precomputation1():
    m=ceil(sqrt(d))
    mP=m*P
    H=[]
    C=[]
    D=[-1]
    phimP=(W(mP,p-1)/W(mP,p-1+d))^(x**2)
    r=phimP/phiP
    for i in range(5):
        H.append(phi3(i*mP)/phimP)
        C.append(H[i]^2)
        if(i>1):
            D.append(H[i-2]*H[i])
    H2_ver=H[2]^(-1)
    half_m=ceil(m/2+3)
    for i in range(5,half_m):
        if(i%2==0):
            l=i/2
            H_i=(C[l-1]*D[l+1]-C[l+1]*D[l-1])*H2_ver
        else:
            l=(i-1)/2
            H_i=C[l]*D[l+1]-C[l+1]*D[l]
        H.append(H_i)
        C.append(H_i^2)
        D.append(H_i*H[i-2])
    for i in range(half_m,m+1):
        if(i%2==0):
            l=i/2
            H_i=(C[l-1]*D[l+1]-C[l+1]*D[l-1])*H2_ver
        else:
            l=(i-1)/2
            H_i=C[l]*D[l+1]-C[l+1]*D[l]
        H.append(H_i)
    for i in range(m+1):
        H[i]=H[i]*r
    return H

def Precomputation2():
    m=ceil(sqrt(d))
    A=[]
    C=[]
    D=[-1]
    for i in range(5):
        A.append(phi2(i))
        C.append(A[i]^2)
        if(i>1):
            D.append(A[i-2]*A[i])
    A2ver=A[2]^(-1)
    for i in range(5,m+2):
        if(i%2==0):
            l=i/2
            A_i=(C[l-1]*D[l+1]-C[l+1]*D[l-1])*A2ver
        else:
            l=(i-1)/2
            A_i=C[l]*D[l+1]-C[l+1]*D[l]
        A.append(A_i)
        C.append(A_i^2)
        D.append(A_i*A[i-2])
    C=C[0:m]
    return [D,C]

def bsgs_eds(Q, E, H, A, B):
    def dbl_ini(BLOCK):
        global ANS
        l=BLOCK.subscript
        C=BLOCK.value[1]^2
        D=BLOCK.value[0]*BLOCK.value[2]
        Hks2la2=hks2_inv*(D*B[l-2]-A[l-2]*C)
        Hks2la1=hks1_inv*(D*B[l-1]-A[l-1]*C)
        SUBSCRIPT=2*l-2
        for k in hashsearch(List,Hks2la2,m, SUBSCRIPT):
            if Q==k*P:
                ANS=k
        for k in hashsearch(List,Hks2la1,m, SUBSCRIPT+1):
            if Q==k*P:
                ANS=k
        return [Hks2la2,Hks2la1]

    def dbl(BLOCK):
        global ANS
        l=BLOCK.subscript
        C=BLOCK.value[1]^2
        D=BLOCK.value[0]*BLOCK.value[2]
        Hks2l=hk_inv*(D*B[l]-A[l]*C)
        Hks2ls1=hka1_inv*(D*B[l+1]-A[l+1]*C)
        SUBSCRIPT=2*l
        for k in hashsearch(List,Hks2l,m, SUBSCRIPT):
            if Q==k*P:
                ANS=k
        for k in hashsearch(List,Hks2ls1,m, SUBSCRIPT+1):
            if Q==k*P:
                ANS=k
        return [Hks2l,Hks2ls1]

    def child_ini(BLOCK_1,BLOCK_2):
        value=dbl_ini(BLOCK_1)+dbl_ini(BLOCK_2)
        subscript=2*BLOCK_1.subscript-1
        Bl1=Block(value[0:3],subscript,1)
        Bl2=Block(value[1:4],subscript+1,0)
        return [Bl1,Bl2]

    def child(BLOCK,LEVEL):
        value=dbl(BLOCK)
        subscript=2*BLOCK.subscript-1
        Bl1=Block(BL[LEVEL+1][1].value[1:3]+[value[0]],subscript,0)
        Bl2=Block(Bl1.value[1:3]+[value[1]],subscript+1,0)
        if(LEVEL==level-2):
            Bl1.mark=1
            Bl2.mark=1
        return [Bl1,Bl2]

    def showBL():
        print(BL[0].value,BL[0].subscript,BL[0].mark)
        for i in range(1,len(BL)):
            print(BL[i][0].value,BL[i][0].subscript,BL[i][0].mark,BL[i][1].value,BL[i][1].subscript,BL[i][1].mark)
        return 0

    def main():
        global ANS
        for i in range(1,level-1):
            if(ANS!=None):
                return ANS
            else:
                BL.append(child_ini(BL[i][0],BL[i][1]))
        BL[level-1][1].mark=1
        L=level-1
        while(ANS==None):
            if(L==0):
                break
            if(BL[L][0].mark==0):
                if(2*BL[L][0].subscript-1>m):
                    BL[L][0].mark=1
                else:
                    BL[L+1]=child(BL[L][0],L)
                    BL[L][0].mark=1
                L=L+1

            elif(BL[L][1].mark==0):
                if(2*BL[L][1].subscript-1>m):
                    BL[L][1].mark=1
                else:
                    BL[L+1]=child(BL[L][1],L)
                    BL[L][1].mark=1
                L=L+1
            else:
                L=L-1
        return ANS
    
    for i in range(-1,4):
        if Q==i*P:
            return i%d
    global ANS
    ANS=None
    m=ceil(sqrt(d))
    level=ceil(log(m, 2))
    table0=H
    table=[]
    for i in range(m):
        table.append(table0[i])
    List=[]
    for i in range(m):
        List.append([])
    for i in range(m):
        List[ezhash(table[i],m)].append(Element(table[i], i))
    hka1=phi(Q+P)
    hka1_inv=hka1^(-1)
    hk=phi(Q)
    hk_inv=hk^(-1)
    hks1= phi(Q-P)
    hks1_inv=hks1^(-1)
    hks2= phi(Q-2*P)
    hks2_inv=hks2^(-1)
    hks3=hka1_inv*(hk*hks2*B[2]-A[2]*hks1^2)
    INI=[hks1, hks2, hks3]
    for i in range(3):
        for k in hashsearch(List, INI[i],m, i):
            if Q==k*P:
                ANS=k
    BL=[Block([hks1, hks2, hks3],2,1)]
    [hks4,hks5]=dbl(BL[0])
    BL.append([Block([hks2, hks3, hks4],3,1),Block([hks3, hks4, hks5],4,0)])
    k=main()
    return k

p=536870951
A=474772770
B=446849777
F=GF(p)
E=EllipticCurve(GF(p),[A,B])
P=E(374494555, 41012407)
k0=21194194
d=P.order()
Q=k0*P
print('E:',E)
print('P:',P,'k:',k0,'Q:',Q)
x=1/d%(p-1)
phiP=phi3(P)
Pre_1=Precomputation1()
Pre_2=Precomputation2()
k=bsgs_eds(Q,E,Pre_1,Pre_2[0],Pre_2[1])
print('BSGS_EDS:',k)
