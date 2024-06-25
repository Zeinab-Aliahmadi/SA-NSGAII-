function pop=CalcCrowdingDistance(pop,F)

    nF=numel(F);
    
    for k=1:nF
        
        Objs=[pop(F{k}).Objectives];
        
        nObj=size(Objs,1);
        
        n=numel(F{k});
        
        d=zeros(n,nObj);
        
        for j=1:nObj
            
            [cj ,so]=sort(Objs(j,:));
            
            d(so(1),j)=inf;
            
            for i=2:n-1
                
                d(so(i),j)=abs(cj(i+1)-cj(i-1))/abs(cj(1)-cj(end));
                
            end
            
            d(so(end),j)=inf;
            
        end
        
        
        for i=1:n
            
            pop(F{k}(i)).CrowdingDistance=sum(d(i,:));
            
        end
        
    end


end

