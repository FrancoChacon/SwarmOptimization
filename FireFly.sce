/*
FireFly Algorithm

    Functions Aviliable:
    
        ReportFireFly
        
        FireFly
        PreEvaluateFireFly
        EvaluateFireFly
*/




function [Report,Table,RProm] = ReportFireFly(TotalOfElements)

    Table=[];
    tprom=0;
    Report=[];

    TC = TotalOfElements;

    TotalExper = 10;
    Dimensions = 10;

    Ixprom=zeros(1,TC);

    printf("Report FireFly Initialized\n\n")
    printf("Total of Elements =%d \n",TC)
    printf("Dimensions size %d \n",Dimensions)

    for exper=1:TotalExper

        printf("Iteration #%d of %d ",exper,TotalExper)

        tic();
        [Ix,Rep,miny]=FireFly(TC,Dimensions);

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



function [Ix,Rep,miny]=FireFly(TC,Dim)

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

        P=Create(TC,10,MR);
    
        [FE] = PreEvaluationFireFly(P,miny,maxy,FE);

        [FE,miny,maxy,Ix]=EvaluationFireFly(P,miny,maxy,Ix,FE);  
       
        Rep=[Rep; miny];

    end


endfunction




function [FE,miny,maxy,Ix,PL,FEL] = EvaluationFireFly(P,miny,maxy,Ix,FE)

    Alfa =1
    Gamma = 1
    Beta = 1

    [TI, D] = size(P);
    X = zeros(1,D);

    for i = 1:TI
        for j = 1:i

            for d = 1:D 
                X(d)=  P(i, d) 

            end

            if FE(i) > FE(j) then
                r = 0
                for d = 1:D 
                    r = r + (P(i, d) - P(j, d))*2
                end
                r = sqrt(r) 
                for d = 1:D 

                    X(d)= X(d) + Beta * exp((Gamma*r)^2)*(P(j, d) - P(i, d)) + Alfa*rand(1,"normal")

                end
            end

            [y, MR] = TestFunction(X,1);

            FE(i)=normalization(y,miny,maxy);

            if y>maxy then
                maxy=y;
            end

            if y<miny then           
                miny=y;
            end
            for d = 1:D 
                P(i, d) = X(d)   
            end
        end
    end

    [B ,k]=gsort(FE)
    IndexSortx= k

    for d = 1:D 

        Ix(d)= P(IndexSortx(1), d)

    end

endfunction


function [FE] = PreEvaluationFireFly(P,miny,maxy,FE)
    
    [Ti,D] = size(P);
    for i = 1:Ti
        for j = 1:D
            X(j) = P(i,j);
        end

        [y, MRR] = TestFunction(X,1);
        FE(i)=normalization(y,miny,maxy);
    end
endfunction

