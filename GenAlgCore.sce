
//Desnormalization
//Interval [0, 1] -> [LMIN, LMAX]


function a=desnormalization(Vn,LMIN,LMAX)

a=Vn*LMAX+LMIN*(1-Vn);

if a>LMAX
    a=LMAX;
end
if a<LMIN
    a=LMIN;
end

endfunction




//Normalization. Interval [LMIN, LMAX] -> [0, 1]

function a=normalization(Van,LMIN,LMAX)

a=(Van-LMIN)/((LMAX-LMIN)+0.000001);

if a>1
    a=1;
end
if a<0
    a=0;
end

endfunction


function [posmax, valmax]=maxp(V)

temp=size(V);
NTB=max(temp);

posmax=0;
valmax=0;
for b=1:NTB
    if V(b)>valmax
        valmax=V(b);
        posmax=b;    
    end
end


endfunction

//Algoritmya Evolutionary computation 2004
//Proyect of simplification of algorithms
//Autor: Luis Torres T.
//All rights reserved
//May 2004
//Evolution Strategies and Genetic Algorithms
//Maxp


function [posmin, valmin]=minp(V)

temp=size(V);
NTB=max(temp);

posmin=0;
valmin=1000000000;
for b=1:NTB
    if V(b)<valmin
        valmin=V(b);
        posmin=b;    
    end
end

endfunction


function ne=NE()

ne = sum(rand(1,12))-6; //Estimation of the normal random variable
    
endfunction


function n=N(mu,sigma)
    
    n = mu + NE()*sigma;
endfunction


//Adjusment in the range
//Copyrights LMTT092010

function [P]=adjust(P,M)
[NTI NTPr] = size(P);
for k=1:NTI
    for pr=1:NTPr
        
        //Checks the inferior range
        if P(k,pr)<M(pr,1)
            P(k,pr)=M(pr,1);
        end
        //Checks the superior range
        if P(k,pr)>M(pr,2)
            P(k,pr)=M(pr,2);
        end
               
    end
end

endfunction

//*********************************************************************************
//******************************  Functions  **************************************
//*********************************************************************************
//*********************************************************************************

//Generalized Sphere model (minimization)
//  -100 <= X(i) <= 100
//Minimum f(0,0,0,0,...0) = 0 
function y = GenSphere(X)

k = max(size(X));

s=0;

for i=1:k
    s = s + X(i)^2;
end
y = s;


endfunction


//Generalized Schwefel's Problem (minimization)
//  -500 <= X(i) <= 500
//Minimum f(420.9687,...420.9687) = -414.9829*k 
function y = GenSchwefel(X)

k = max(size(X));
s = 0;
for i=1:k
    s=s+X(i)*sin(sqrt(abs(X(i))));
end

y=-s;

endfunction

//Generalized Rastrigins function (minimization)
//  -5.12 <= X(i) <= 5.12
//Minimum f(0,0,0,0,...0) = 0 
function y = GenRastrigins(X)

k = max(size(X));

s=0;

for i=1:k
    s = s + X(i)^2-10*cos(2*%pi*X(i))+10;
end
y = s;
endfunction

//Generalized Rosenbrock function (minimization)
//  -10 <= X(i) <= 10
//Minimum f(1,1,1, 1,...1) = 0 

function y = GenRosenbrock(X)

k = max(size(X));
s=0;

for i=1:k-1
    s = s + ((100*(X(i+1) - X(i)^2)^2) + (1 -  X(i))^2);
end
y = s;


endfunction

//Generalized Schwefel's Problem 1.2 (minimization)
//  -100 <= X(i) <= 100
//Minimum f(0,0,0,0,0,...0) = 0 
function y = GenSchwefel12(X)
k = max(size(X));
s=0;

for i=1:k
    s2=0;
    for j=1:i
        s2 = s2 + X(j);
    end
    s = s + s2^2;
end

y = s;

endfunction




//Generalized Grienwangk's function
//  -600 <= X(i) <= 600
//Minimum f(0,0,0,0,0,...0) = 0 
function y = GenGriewangk(X)
n = max(size(X));

    s1=0;
    for j=1:n
        s1 = s1 + (X(j)^2);
    end
    s2=1;
    for j=1:n
        s2 = s2 * (cos(X(j)/sqrt(j))); 
    end
    
    y = (1/4000)*s1 - s2 + 1;


endfunction





//Generalized Ackley's function
//  -32.768 <= X(i) <= 32.768
//Minimum f(0,0,0,0,0,...0) = 0 
function y = GenAckley(X)
    //Recpommended:
    a=20; b=0.2; c=2*%pi;
n = max(size(X));

    s1=0;
    for j=1:n
        s1 = s1 + (X(j)^2);
    end
    
    s2=0;
    for j=1:n
        s2 = s2 + cos(c*X(j))
    end
    
    y = -a*exp(-b*sqrt((1/n)*s1)) - exp((1/n)*s2) + a + exp(1);

endfunction


//Generalized Langermann's function
//  -1 <= X(i) <= 10
//Local minimum are unevenly distributed 
function y = GenLangermann(X)
    //Recpommended:
   
    n = max(size(X));
    a = [3 5 2 1 7]; b = [5 2 1 4 9]; c = [1 2 5 2 3];



s=0;
for ii=1:5


    s1=0;
    for j=1:n
        s1 = s1 + (X(j)-a(ii))^2;
    end
    
    s2=0;
    for j=1:n
        s2 = s2 + (X(j)-b(ii))^2
    end
    
    s = s + c(ii)*exp(-(1/%pi)*s1)*cos(%pi*s2);
    
    
    
end //for ii

    y=s;
    
endfunction


//Generalized Michalewics's function
//  0 <= X(i) <= pi
//global minimum n=5 is -4.687
//n=10 is -9.66
//other number is unknown
function y = GenMichalewicz(X)
    //Recommended: m=1; 
    //Larger m becomes like search  a needle in the haystack...
   m=1;
    n = max(size(X));
  
    s=0;
    for j=1:n
        s = s + sin(X(j))*(sin(j*X(j)^2/%pi))^(2*m);
    end

    y=-s;
    
endfunction



//Deceptive function


function y=g(x,i,alpha)
    
    if ((x>=0) & (x<((4/5)*alpha(i)))) then    
        y=(4/5)-x/alpha(i);
    end
    
    if (x>((4/5)*alpha(i)) & (x<=alpha(i))) then    
        y=5*x/alpha(i)-4;
    end
    
    if ((x>alpha(i)) & (x<=((1+4*alpha(i))/5))) then    
        y=(5*(x-alpha(i)))/(alpha(i)-1)+1;
    end
    
    if ((x>((1+4*alpha(i))/5)) & (x<=1)) then    
        y=((x-1)/(1 - alpha(i)))+4/5;
    end
    
    
endfunction


//FIrst case 
function y = GenDeceptive(X,alpha)
   
    vbeta=0.2;
    n = max(size(X));  
    s=0;
    for j=1:n
        s = s + g(X(j),j,alpha)
    end
    y=-(s)^vbeta/n;
    
endfunction


function [y, MR] = TestFunction(X,option)
    
MR=zeros(1,2);

select option

case 1 then
//GenSphere
MR(1)=-5.12;
MR(2)=5.12;

//FireFly

case 2 then
//GenSchwefel
MR(1)=-500;
MR(2)=500;

case 3 then
//GenRastrigins
MR(1)=-5.12;
MR(2)=5.12;

case 4 then
//GenRosenbrock
MR(1)=-10;
MR(2)=10;

case 5 then
//GenSchwefel12
MR(1)=-100;
MR(2)=100;

case 6 then
//GenGriewangk
MR(1)=-600;
MR(2)=600;

case 7 then
//GenAckley
MR(1)=-32.768;
MR(2)=32.768;

case 8 then
//GenLangermann
MR(1)=-1;
MR(2)=10;

case 9 then
//GenMichalewicz
MR(1)=0;
MR(2)=%pi;

case 10 then
//Deceptive functions
MR(1)=0;
MR(2)=1;

end
        
        select option
        
        case 1 then
        y = GenSphere(X);        
        
    case 11 then
        y = FireFly(X);
        
        case 2 then
        y = GenSchwefel(X);
        
        case 3 then
        y = GenRastrigins(X);
        
        case 4 then
        
        y = GenRosenbrock(X);
        
        case 5 then
       y = GenSchwefel12(X);
        
        case 6 then
        y = GenGriewangk(X);
        
        case 7 then
        y = GenAckley(X);
        
        case 8 then
        y = GenLangermann(X);
       
        case 9 then
        y = GenMichalewicz(X);
       
        case 10 then
       //Deceptive function
       alpha=[0.3 0.7];//Solution
       y = GenDeceptive(X,alpha);
       
       end
    
    
endfunction



//************************  END FUNCTIONS  ****************
//*********************************************************
//*********************************************************


function P=Create(TI,TR,MR)
    
    P=zeros(TI,TR);
    
    for k=1:TI
        for r=1:TR
           P(k,r) = desnormalization(rand(),MR(r,1),MR(r,2)); 
        end
        
    end
    
endfunction


function [FE,miny,maxy,Ix]=Evaluation(P,miny,maxy,Ix)
    
    [TI, TR] = size(P);FE=zeros(1,TI);
    
    for k=1:TI
        for r=1:TR
            X(r)=P(k,r);
         end
           [y, MR] = TestFunction(X,2);
         
         if y>maxy then
             maxy=y;
            //   Ix=X; //best solution found (maximization))
         end
         
         if y<miny then
             miny=y;
             Ix=X; ////best solution found (minimization)
         end
         
        FE(k)=normalization(y,miny,maxy);
        
        
    end
    
    
endfunction


function [Ix,Rep,miny]=RandomS(N)
    
    miny=100000000;
    maxy=-100000000;
    //GenSphere
    MR=zeros(10,2);
    MR(:,1)=-5.12;
    MR(:,2)=5.12;
    Rep=[];
    Ix=rand(1,10);
    
    for cycles=1:N
        
         P=Create(100,10,MR);
         
        [FE,miny,maxy,Ix]=Evaluation(P,miny,maxy,Ix);
        
        Rep=[Rep; miny];
        
    end
    
    
endfunction



//***********************************************************************
//*************** Reports  **********************************************
//*********************************************************************

function [Report,Table,Rprom] = ReportRandomS()
    
    Table=[];
    tprom=0;
    Report=[];
    Ixprom=zeros(1,10);
    
    for exper=1:10
        
        tic();[Ix,Rep,miny]=RandomS(100);a=toc();
        tprom = tprom + a;
        Ixprom=Ixprom+Ix';
        Report=[Report; Rep'];
    end
    
    tprom = tprom / 10;
    
    
    //pause,
    //Execution time
    
    //Mean and standard deviation
    Rprom=mean(Report,'r');
    Table = [tprom Rprom(100) Ixprom];
    
//    clf
//    plot(Report'); 
//    plot(Rprom','*k');
    

endfunction

