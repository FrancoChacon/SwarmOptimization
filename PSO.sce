/*
PSO ALGORITHM

    Functions Aviliable:
        ReportPSO
        PSO
        EvaluatePSO
        UpdatePSO

*/

function [Report,Table,RProm] = ReportPSO()
    
    Table=[];
    tprom=0;
    Report=[];
    Ixprom=zeros(1,10);
    for exper=1:10
        
        tic();[B,Rep,miny]=PSO(100);a=toc();
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


function [B,Rep,miny]=PSO(TC)
    
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

function [P,V]=UpdatePSO(P,V)
    
    w=0.3;
    beta1=0.35;
    beta2=0.35;
     [Ti, D] = size(P);
     PL=zeros(Ti,10); //Velocidad
     B=zeros(Ti,1); //Mejor evaluación local
     
   
    for k=1:Ti
        for d=1:D
               V(k,d) = V(k,d)*w + beta1*rand()*(PL(k,d)-P(k,d))+beta2*rand()*(B(d)-P(k,d));
               
               P(k,d) = P(k,d) + V(k,d); 
            
        end
          
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

