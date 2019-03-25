/*
FireFly ALGORITHM

    Functions Aviliable:
        ReportFireFly
        FireFly
        EvaluateFireFly
        UpdateFireFly

*/




function [Report,Table,RProm] = ReportFireFly()
    
    Table=[];
    tprom=0;
    Report=[];
    Ixprom=zeros(1,10);
    for exper=1:10
        
        tic();[B,Rep,miny]=FireFly(100);
        
        a=toc();
        tprom = tprom + a;
        Bprom=Ixprom+B';
        Report=[Report; Rep'];
    end
    
    tprom = tprom / 10;
    
    
    //pause,
    //Execution time
    
    //Mean and standard deviation
    RProm=mean(Report,'r');
    Table = [tprom RProm(100) Ixprom];
    
//    clf
//    plot(Report'); 
//    plot(RProm','*k');
    

endfunction



function [B,Rep,miny]=FireFly(TC)
    
    miny=100000000;
    maxy=-100000000;
   
    MR=zeros(10,2);
    MR(:,1)=10//LUB(1);
    MR(:,2)=-1//LUB(2);
    Rep=[];
    
    
    //Poblacion inicial
     Ti=100; // Total de individuos
     P=Create(100,10, MR);
     B=P(1,:);  //Mejor solución global
     PL=P;   //Mejor posición local
     V=zeros(Ti,10); //Velocidad
     FEL=zeros(Ti,1); //Mejor evaluación local
     
     
    for cycles=1:TC
        
        
         //Evaluacion
        [FE,miny,maxy,B,PL,FEL]=EvaluationPSO(P,miny,maxy,B,PL,FEL);
        
        
        //Actualizacón de posición
        [P,V]=UpdatePSO(P,V);
      
        
        Rep=[Rep; miny];
    
    
        disp(miny)        
    end
    
    
endfunction


function [FE,miny,maxy,B,PL,FEL]=EvaluationPSO(P,miny,maxy,B,PL,FEL)
    
    [TI, D] = size(P);FE=zeros(1,TI);
    
    for k=1:TI
        for d=1:D
            X(d)=P(k,d);
         end
           [y, MR] = TestFunction(X,2);
         
         if y>maxy then
             maxy=y;
            //   Ix=X; //best solution found (maximization))
         end
         
         if y<miny then
             miny=y;
             B=X; //best solution found (minimization)
         end
         
         //Minimization
        FE(k) = 1 - normalization(y,miny,maxy);
        
        if FE(k)>FEL(k) then
            
            FEL(k)=FE(k);
            for d=1:D
                    PL(k,d) = X(d);
            end
        end
            
    end
    
    
endfunction


function [FE,miny,maxy,B,PL,FEL] = EvaluationFireFly(P,miny,maxy,B,PL,FEL)
    
    Alfa =1
    Gamma = 1
    Beta = 1
    
    [TI, D] = size(P);
    FE=zeros(1,TI);

  
    for i = 1:TI
        for j = 1:D
            if FE(i) > FE(j) then
                r= 0
                for d = 1:D 
                    r= r + (P(i, d) - P(j, d))*2
                end
                r=P(r)
                for d = 1:D 
                    X(d)= X(d) + Beta * exp((Gama*r)^2)*(P(j, d) - P(i, d)) + Alfa(o, 1) 
                end
            end
            FE(i) =  TestFunction(X,1);
            
            for d = 1:D 
                P(i, d)= X(d)
            end
        end
    end
    Ix= sort(FE)
    for d = 1:D 
        B(d)= P(Ix(1), d)
    end


endfunction




















function FireFlyAlgInitialization()
    Ti = 1;
    D = 1;

    P= [];
    FE = [];
    X= [];
    Ri = [];
    for i = 1:Ti 
        for d = 1:D 
            P(i, d) = Ri(d, UB(d), LB(d)); //Duda
            X(d) = P(i, d);
        end
        FE(i) = evaluation(X); //Duda

    end

endfunction

