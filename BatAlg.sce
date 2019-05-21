/*
Bat Algorithm
    Functions Aviliable:
        ReportBat
        Bat
        PreEvaluateBat
        EvaluateBat
*/


function [Report,Table,RProm] = ReportBat(TotalOfElements)

    Table=[];
    tprom=0;
    Report=[];

    TC = TotalOfElements;

    TotalExper = 10;
    Dimensions = 10;

    Ixprom=zeros(1,TC);

    printf("Report Bat Initialized\n\n")
    printf("Total of Elements =%d \n",TC)
    printf("Dimensions size %d \n",Dimensions)

    for exper=1:TotalExper

        printf("Iteration #%d of %d ",exper,TotalExper)

        tic();
        [Ix,Rep,miny]=Bat(TC,Dimensions);

        a=toc();

        printf("- %f Seg\n",a)

        tprom = tprom + a;
        Ixprom=Ixprom+Ix';

        Report=[Report; Rep'];

    end

    tprom = tprom / 10;
    //Mean and standard deviation
    RProm=mean(Report,'r');
    Table = [tprom RProm(TC) Ixprom];

    printf("\n")

endfunction



function [Ix,Rep,miny]=Bat(TC,Dim)

    //Poblacion inicial
    Ti=TC;   // Total of Elements
    D = Dim; //Total of Dimensions

    //Variable Definitions
    miny=100000000;
    maxy=-100000000;
    Rep=[];
    FE=zeros(Ti,1);
    Ix=rand(Ti,1);

    //GenSphere
    MR=zeros(10,2);
    MR(:,1)=-5.12;
    MR(:,2)=5.12;
    MRR = MR
    
 
    
    for cycles=1:TC

        P=Create(TC,D,MR);

        [FE,f,r,A] = PreEvaluationBat(P,miny,maxy)

        [FE,miny,maxy,Ix,f,r,A]=EvaluationBat(P,miny,maxy,Ix,f,r,A,FE);  
       
        Rep=[Rep; miny];

    end


endfunction




function [FE,miny,maxy,Ix,f,r,A]=EvaluationBat(P,miny,maxy,Ix,f,r,A,FE)


    [TI, D] = size(P);
    X = zeros(1,D);
    fex = zeros(TI,D);
    feb = zeros(TI,D);

    Alfa =1
    Gamma = 1
    Beta = 1

    fmax = [1,D]
    fmin = [1,D]

    fmax(1:D) = 100;
    fmin(1:D) = 0;


    //{Block: Generate new solutions adjusting frequency}; 
    [FE,B,X,f,maxy,miny,P] = GenerateSolutionsByFreq(P,miny,maxy,fmax,fmin,f)
for i = 1:Ti 
    if rand(1,1,"uniform") > r(i) 
        //{Block: Select the best solution and generate a local one}; 
        [B,P] = SelectlAndGenerateLocal(FE,P)
    end
    //{Generate and evaluate a new solution};
    [y, MR] = TestFunction(X,1);
    fex(i)=normalization(y,miny,maxy);

    [y, MR] = TestFunction(B,1);
    feb(i)=normalization(y,miny,maxy);

    for i = 1:Ti 
        if rand(1,1,"uniform") < A(i) then 
            if  fex < feb  
                feb = FE(i); 
                B = X;
                A(i) = alfa * A(i); 
                r(i) = r(i)*(1 - exp(-Gamma * i))
            end
        end
    end
  
    //{Rank the bats and select the best one}; 
    [pos,val] = maxp(FE); 
    for d = 1:D 
        Ix(d) = P(pos,d)
    
    end
end
endfunction


function [FE,B,X,f,maxy,miny,P] = GenerateSolutionsByFreq(P,miny,maxy,fmax,fmin,f)
    [TI, D] = size(P)
    X = zeros(1,D)
    V = zeros(TI,D)
    B = zeros(1,D)

    for i = 1:Ti
        for d = 1:D 
            X(d) = P(i,d) 
        end 
        [y, MR] = TestFunction(X,1);
        FE(i) = normalization(y,miny,maxy);

        if y > maxy 
            maxy = y 
        end
        if y < miny 
                B = X
            miny = y; 
        end
    end 
    NN = Ti;
    //{Generate a new solution by ﬂying randomly}
    for i = 1:NN 
        for d = 1:D 
            f(i,d) = fmin(d) + rand(1,1,"uniform")*(fmax(d)-fmin(d));          
            V (i,d) = V (i,d) + (P(i,d)-B(d));
            P(i,d) = P(i,d) + V (i,d); 
        end
    end

endfunction


function [B,P] = SelectlAndGenerateLocal(FE,P)

    [pos,val] = maxp(FE)
    [TI, D] = size(P)
    for d = 1:D 
        B(d) = P(pos,d); 
    end 
    Ap = 0;
    for i = 1:Ti  
        Ap = Ap + A(i);
    end
     Ap = Ap/Ti; 
    for i = 1:Ti 
        for d = 1:D  
            P(i,d) = P(i,d) + Ap * rand(1,1,"uniform") 
        end
    end

endfunction


function [FE,f,r,A] = PreEvaluationBat(P,miny,maxy)
    
    [Ti,D] = size(P);
       //Bat Alg Variables
    f=zeros(Ti,D);
    r=zeros(Ti,1);
    A=zeros(Ti,1);   
    
    for i = 1:Ti
        for j = 1:D
            X(j) = P(i,j);
        end

        f(i) = 100; 
        r(i) = 0.5; 
        A(i) = 2; 
       
    end
endfunction

