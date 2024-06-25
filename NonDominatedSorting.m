function [pop,F]=NonDominatedSorting(pop)

    npop=numel(pop);

    for i=1:npop
        pop(i).DominationSet=[];
        pop(i).DominatedCount=0;
    end
    
    F{1}=[];
    
    for i=1:npop
        for j=i+1:npop
            p=pop(i);
            pp=pop(j);
            
            if Dominates(p.Objectives,pp.Objectives)
                p.DominationSet=[p.DominationSet j];
                pp.DominatedCount=pp.DominatedCount+1;
            end
            
            if Dominates(pp.Objectives,p.Objectives)
                pp.DominationSet=[pp.DominationSet i];
                p.DominatedCount=p.DominatedCount+1;
            end
            
            pop(i)=p;
            pop(j)=pp;
        end
        
        if pop(i).DominatedCount==0
            F{1}=[F{1} i];
            pop(i).Rank=1;
        end
    end
    
    k=1;
    
    while true
        
        Q=[];
        
        for i=F{k}
            p=pop(i);
            
            for j=p.DominationSet
                pp=pop(j);
                
                pp.DominatedCount=pp.DominatedCount-1;
                
                if pp.DominatedCount==0
                    Q=[Q ,j];
                    pp.Rank=k+1;
                end
                
                pop(j)=pp;
            end
        end
        
        if isempty(Q)
            break;
        end
        
        F{k+1}=Q;
        
        k=k+1;
        
    end
    

end