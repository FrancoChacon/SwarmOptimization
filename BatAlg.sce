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

        [FE,miny,maxy,Ix]=EvaluationBat(P,miny,maxy,Ix,f,r,A);  
       
        Rep=[Rep; miny];

    end


endfunction




function [FE,miny,maxy,Ix]=EvaluationBat(P,miny,maxy,Ix,f,r,A); 

    Alfa =1
    Gamma = 1
    Beta = 1

    [TI, D] = size(P);
    X = zeros(1,D);


    //{Block: Generate new solutions adjusting frequency}; 
   // GenerateSolutionsByFreq(P);
    
    
    for i = 1:Ti 
        if U() > r(i) 
            //{Block: Select the best solution and generate a local one}; 
            //SelectlAndGenerateLocal();
        end
        //{Generate and evaluate a new solution};
        fex = evaluation(X); 
        feb = evaluation(B);
        for i = 1:Ti 
            if U() < A(i) and fex < feb 
                FEB = FE(i); 

                B = X;
                A(i) = α·A(i); 
                r(i) = r(i)(1−exp−γ·cycle)
            end
        end
        //{Rank the bats and select the best one}; 
        (pos,val) = maxp(FE); 
        for d = 1:D do

            B(d) = P(posmax,d)
        end




    endfunction



function [FE,B] = GenerateSolutionsByFreq(P)
    [TI, D] = size(P)
    X = zeros(1,D)

    for i = 1:Ti
        for d = 1:D 
            X(d) = P(i,d) 
        end 


        [y, MR] = TestFunction(X,1);
        FE(i)=normalization(y,miny,maxy);

        
        if FE(i) > maxy 
            maxy = FE(i) 
            B = X
        end
        if FE(i) < miny 
            miny = FE(i); 
        end
    end 
    
    //Duda NN
    NN = Ti;
    //{Generate a new solution by ﬂying randomly}
    for i = 1:NN 
        for d = 1:D 
              f(i,d) = fmin(d) + U()·(fmax(d)−fmin(d)); DUDA
            //f(i,d) = fmin(d) + U()·(fmax(d)−fmin(d)); DUDA
            V (i,d) = V (i,d) + (P(i,d)−B(d));
            P(i,d) = P(i,d) + V (i,d); 
        end
    end


endfunction


function SelectlAndGenerateLocal()
    (val,pos) = maxp(FE))
    for d = 1:D 
        B(d) = P(pos,d); 
    end 
    Ap = 0;
    for i = 1:Ti  
        Ap = Ap + A(i);
    end Ap = Ap/Ti; 
    for i = 1:Ti 
        for d = 1,D  
            P(i,d) = P(i,d) + Ap ·Un() 
        end
    end

endfunction


function [FE,f,r,A] = PreEvaluationBat(P,miny,maxy)
    
    [Ti,D] = size(P);
       //Bat Alg Variables
    f=zeros(Ti,1);
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

