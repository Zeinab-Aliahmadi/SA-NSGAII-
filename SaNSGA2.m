clc;
clear;
close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                           Copyright rules: To use the current code,                                                                   %                        %
%                                                 Cite the following paper:                                                                             %                       %
%                                   Seyed Zeinab Aliahmadi, Armin Jabbarzadeh, Lucas A. Hof (2024),                                                     %  
%       A Multi-objective Optimization Approach for Sustainable and Personalized Trip Planning:  A Self-adaptive Evolutionary Algorithm with Case Study %
%                                                   Expert Systems With Applications                                                                    %                                                                                                                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Problem Definition
filename='CaseStudy.xlsx';       

Nh=xlsread(filename,1,'B:B');    %The total number of available hotels (H)%
Ns=xlsread(filename,1,'D:D');    %The total number of available scenic spots(S)%
Nr=xlsread(filename,1,'F:F');    %The total number of available restaurants(R)%
Nsr=Ns+Nr;                       %The total number of available scenic spots+restaurants(R+S)%
N=Nh+Ns+Nr;                      %The total number of nodes(n)%
Nv=xlsread(filename,1,'H:H');    %The total number of days included in a tour(T) %
t=xlsread(filename,1,'J:J');     %The total number of days included in a tour(T) %
%Distance (N);
xr=xlsread(filename,1,'N:N');
yr=xlsread(filename,1,'O:O');

for i=1:numel(xr)
    for j=1:numel(yr)
        dis(i,j)=sqrt((yr(i)-yr(j))^2 + (xr(i)-xr(j))^2);
    end
end

%data has been given on: https://open.canada.ca/data/en/dataset/02ebdab9-cbf3-4f56-8c29-79fa0ed0ed2e
NL=zeros(N,N);
TL1=xlsread(filename,1,'R:R');
TL2=xlsread(filename,1,'S:S');
for i=1:N
    for j=1:N
    for k=1:numel(TL1)
      if TL1(k)>=xr(i) && TL1(k)<=xr(j)
          if TL2(k)>=yr(i) && TL2(k)<=yr(j)
          NL(i,j)=NL(i,j)+1;
          end
      end
    end
    end
end

%The Lower bound of time windows to visit a place (Nsr)%  
LTW=xlsread(filename,1,'W:W');

%The Uper bound of time windows to visit a  place (Nsr)%
UTW=xlsread(filename,1,'AA:AA');

%Operational time to visite a place (Nsr)%
OT1=xlsread(filename,1,'AE:AE');
OT2=xlsread(filename,1,'AI:AI');
OT3=xlsread(filename,1,'AM:AM');

%Global (utility) rate of each restaurant (Nr)%
UT=xlsread(filename,1,'AQ:AQ');

%Restaurant's food & service cost (Nr)%
RFC1=xlsread(filename,1,'AU:AU');
RFC2=xlsread(filename,1,'AY:AY');
RFC3=xlsread(filename,1,'BC:BC');

%Global (utility) rate of each hotel (Nh)%
UTH=xlsread(filename,1,'BG:BG');

%The booking cost of each available hotel  per day per person (Nh)%
RCH1=xlsread(filename,1,'BK:BK');
RCH2=xlsread(filename,1,'BO:BO');
RCH3=xlsread(filename,1,'BS:BS');
%The ticket price of each listed SS  per day per person (Ns)%
SCs=xlsread(filename,1,'DH:DH');

%The variable cost of each vehicle for 1 unit distance%
VC=xlsread(filename,1,'BW:BW');

%The variable cost (ticket price) of each vehicle for 1 unit distance%
FC=xlsread(filename,1,'CA:CA');

%average velocity of each vehicle (km/h)%
v_star1=xlsread(filename,1,'CE:CE');
v_star2=xlsread(filename,1,'CI:CI');
v_star3=xlsread(filename,1,'CM:CM');

%carbon emission of each vehicle for 1 unit distance per passenger
C=xlsread(filename,1,'CQ:CQ');

%Budget%
B=xlsread(filename,1,'CT:CT');

%Departure time of visitor from the hotel%
DT=xlsread(filename,1,'CW:CW');
%Average waiting time at a traffic light%
wt=xlsread(filename,1,'CZ:CZ');
%Maximum number of trips per day=maximum number of visited SS per day%
MT=xlsread(filename,1,'DC:DC');

%% FUZZY Probabilities
alfaOT=0.5;
alfaRCH=0.5;
alfaRFC=0.5;
alfav_star=0.5;

%% NSGA-II Parameters
npop=100;
MaxIt=100;    
%% Number of Objective Functions
Objectives=@(pop) objFunction(pop);      % objective Function
nObj=3;                                            

    sol.L=[];
    sol.Rank=[];
    sol.DominationSet=[];
    sol.DominatedCount=[];
    sol.CrowdingDistance=[];
    sol.At=zeros(t,MT+1);
    sol.TV=zeros(t,MT+1);
    pop=repmat(sol,npop,1);
       
tic;

    LOT=(2*alfaOT-1)*OT1+(2-2*alfaOT)*OT2;
    UOT=(2-2*alfaOT)*OT2+(2*alfaOT-1)*OT3;
    LRCH=(2*alfaRCH-1)*RCH1+(2-2*alfaRCH)*RCH2;
    URCH=(2-2*alfaRCH)*RCH2+(2*alfaRCH-1)*RCH3;
    LRFC=(2*alfaRFC-1)*RFC1+(2-2*alfaRFC)*RFC2;
    URFC=(2-2*alfaRFC)*RFC2+(2*alfaRFC-1)*RFC3;
    Lv_star=(2*alfav_star-1)*v_star1+(2-2*alfav_star)*v_star2;
    Uv_star=(2-2*alfav_star)*v_star2+(2*alfav_star-1)*v_star3;
     
    npop=2*npop;
for i=1:npop
    
     pop(i).H=randi(Nh,1);
     pop(i).S=randperm(Ns); 
     pop(i).R=randperm(Nr);
     pop(i).R=pop(i).R(1:t); 
     pop(i).SS=pop(i).S; 
for k=1:t
    if k<t
         r(k)=randi([2,MT]);                   
         for j=1:MT+1
         if j<r(k) ||  j>r(k) 
            pop(i).Lsr(k,j)=pop(i).SS(1);
            pop(i).SS(1)=[];
         end
         if j==r(k)   
            pop(i).Lsr(k,j)=Ns+pop(i).R(k);
         end                 
         end
    end
    
    if k==t
         r(t)=randi([2,numel(pop(i).SS)+1]);
         for j=1:numel(pop(i).SS)+1
         if j<r(t) ||  j>r(t) 
            pop(i).Lsr(t,j)=pop(i).SS(1);
            pop(i).SS(1)=[];
         end
         if j==r(t)   
            pop(i).Lsr(t,j)=Ns+pop(i).R(t);
         end                 
         end    
    end
end
     
for k=1:t
    m(k)=0;
    for j=1:numel(pop(i).Lsr(k,:))
    if pop(i).Lsr(k,j)==0
       m(k)=m(k)+1;
    end
    end
end 

for k=1:t 
    pop(i).L(k,1)=pop(i).H; 
for j=1:numel(pop(i).Lsr(k,:)) 
    if pop(i).Lsr(k,j)~=0
    pop(i).L(k,j+1)=Nh+pop(i).Lsr(k,j);
    else
    pop(i).L(k,j+1)=0;    
    end 
     pop(i).L(k,numel(pop(i).Lsr(k,:))-m(k)+2)=pop(i).H; 
end
end

if Nv==5
    pop(i).V=randi(Nv-1,t,numel(pop(i).L(1,:))-1);
for k=1:t
  for j=1:numel(pop(i).Lsr(k,:))-m(k)  
    if dis(pop(i).L(k,j),pop(i).L(k,j+1))/Uv_star(pop(i).V(k,j))<=0.005
    pop(i).V(k,j)=5;    
    end
  end
end

else
        pop(i).V=randi(Nv,t,numel(pop(i).L(1,:))-1);
end

for k=1:t
  for j=1:numel(pop(i).Lsr(k,:))-m(k)  
      if j==1 && pop(i).V(k,j)~=4
   pop(i).At(k,j)=DT+(dis(pop(i).L(k,j),pop(i).L(k,j+1))/Uv_star(pop(i).V(k,j))) + (wt*NL(pop(i).L(k,j),pop(i).L(k,j+1)));
   pop(i).TV(k,j)=pop(i).At(k,j)+UOT(pop(i).Lsr(k,j));  
      elseif j>1  && pop(i).V(k,j)~=4  
   pop(i).At(k,j)=pop(i).TV(k,j-1)+ (dis(pop(i).L(k,j),pop(i).L(k,j+1))/Uv_star(pop(i).V(k,j)))+ (wt*NL(pop(i).L(k,j),pop(i).L(k,j+1))); 
   pop(i).TV(k,j)=pop(i).At(k,j)+ UOT(pop(i).Lsr(k,j)); 
      elseif j==1 && pop(i).V(k,j)==4
   pop(i).At(k,j)=DT+(dis(pop(i).L(k,j),pop(i).L(k,j+1))/Uv_star(pop(i).V(k,j))) ;
   pop(i).TV(k,j)=pop(i).At(k,j)+UOT(pop(i).Lsr(k,j));  
      elseif j>1  && pop(i).V(k,j)==4  
   pop(i).At(k,j)=pop(i).TV(k,j-1)+ (dis(pop(i).L(k,j),pop(i).L(k,j+1))/Uv_star(pop(i).V(k,j))); 
   pop(i).TV(k,j)=pop(i).At(k,j)+ UOT(pop(i).Lsr(k,j)); 

      end
  end
end

    for k=1:t  
    for j=1:numel(pop(i).Lsr(k,:))-m(k) 
    if pop(i).At(k,j)<LTW(j) 
        if pop(i).TV(k,j)>UTW(j)|| pop(i).At(k,j)>UTW(j)  
           pop(i).L(:,:)=[];    
        end
    end
    end
    end
    
if ~isempty(pop(i).L(:,:))
for k=1:t
    for j=1:numel(pop(i).Lsr(k,:))-m(k)+1
    if pop(i).L(k,j)~=0 && pop(i).L(k,j+1)~=0     
       pop(i).TotalDD(k,j)= dis(pop(i).L(k,j),pop(i).L(k,j+1));
       pop(i).TotalTCC(k,j)=VC(pop(i).V(k,j))*pop(i).TotalDD(k,j)+FC(pop(i).V(k,j));
    end
    end
end
for k=1:t
   for j=1:numel(pop(i).L(k,:))
     pop(i).TotalD(k)=sum((pop(i).TotalDD(k,:))); 
     pop(i).TotalTC(k)=sum((pop(i).TotalTCC(k,:)));
   end
   pop(i).TotalD=sum(pop(i).TotalD(:));
   pop(i).TotalTC=sum(pop(i).TotalTC(:));
end

for k=1:t
pop(i).TotalRC(k)= URFC(pop(i).R(k));
pop(i).TotalUTR(k)=UT(pop(i).R(k));
end

pop(i).TotalRC=sum(pop(i).TotalRC);
pop(i).TotalUTR=sum(pop(i).TotalUTR);

pop(i).TotalC=t*URCH(pop(i).L(1,1))+ pop(i).TotalTC + pop(i).TotalRC+SCs;

if pop(i).TotalC<=B
pop(i).TotalC=pop(i).TotalC;
pop(i).TotalUT=t*UTH(pop(i).L(1))+pop(i).TotalUTR;
    for j=1:numel(pop(i).V)
        pop(i).TotalEm(j)=pop(i).TotalDD(j)*C(pop(i).V(j));
    end
pop(i).TotalEm=sum(pop(i).TotalEm);
else
pop(i).L=[]; 
end

if ~isempty(pop(i).L)
pop(i).Objectives=objFunction(pop(i));  
end     
end
end
for i=npop:-1:1
if isempty(pop(i).L)
pop(i)=[];
end
end
npop=numel(pop);

for i=1:numel(pop)-1
    for j=numel(pop):-1:i
        if pop(i).TotalC==pop(j).TotalC
           pop(j)=[];
        end
    end
end
npop=numel(pop);
%% Selfadaptive
f1min=pop(1).TotalC;
f2min=pop(1).TotalUT;
f3min=pop(1).TotalEm;
f1max=0;
f2max=0;
f3max=0;

for i=1:npop
    if pop(i).TotalC<f1min
   f1min=pop(i).TotalC;
    end
    if pop(i).TotalUT<f2min
   f2min=pop(i).TotalUT;
    end
    if pop(i).TotalEm<f3min
   f3min=pop(i).TotalEm;
    end
    if pop(i).TotalC>f1max
   f1max=pop(i).TotalC;
    end
    if pop(i).TotalUT>f2max
   f2max=pop(i).TotalUT;
    end
    if pop(i).TotalEm>f3max
   f3max=pop(i).TotalEm;
    end
end


for i=1:npop
   f1(i)= (f1max-pop(i).TotalC)/(f1max-f1min);
   f2(i)= (f2max-pop(i).TotalUT)/(f2max-f2min);
   f3(i)= (f3max-pop(i).TotalEm)/(f3max-f3min);
end

f1m=f1(1)-pop(1).TotalC;
f2m=f2(1)-pop(1).TotalUT;
f3m=f3(1)-pop(1).TotalEm;

for i=1:npop
  if  f1(i)-pop(i).TotalC<f1m
     f1m= f1(i)-pop(i).TotalC;
  end
  if  f2(i)-pop(i).TotalUT<f2m
     f2m= f2(i)-pop(i).TotalUT;
  end
  if  f3(i)-pop(i).TotalEm<f3m
     f3m= f3(i)-pop(i).TotalEm;
  end
end

sigma1=0;
sigma2=0;
sigma3=0;
for i=1:npop
   sigma1=sigma1+f1(i);
   sigma2=sigma2+f2(i); 
   sigma3=sigma3+f3(i); 
end
sigma1=(sigma1/npop);
sigma2=(sigma2/npop);
sigma3=(sigma3/npop);
k=0.6640625;
g=0.3203125;

pCrossover=k*cos((pi/2)*(1/(exp(sigma1+sigma2+sigma3))));
nCrossover=2*round(pCrossover*npop/2);      % Number of Offsprings

pMutation=1;
for i=1:npop
pMutation=pMutation+((f1(i)/f1m)*(f2(i)/f2m)*(f3(i)/f3m));
end
pMutation=g*pMutation;
nMutation=round(pMutation*npop); 
%%

% Non-Dominated Sorting
[pop, F]=NonDominatedSorting(pop);

% Calculate Crowding Distance
pop=CalcCrowdingDistance(pop,F);

% Sort Population
[pop, F]=SortPopulation(pop);

 
%% NSGA-II Main Loop

for it=1:MaxIt
%% Crossover
popc=repmat(sol,nCrossover,2);

 for i=1:nCrossover/2 
  
k=randi(npop);
f=randi(npop);

popc(i,1).H=pop(i).H;
popc(i,2).H=pop(i).H;

popc(i,1).S=pop(k).S;
popc(i,2).S=pop(f).S;

popc(i,1).R=randperm(Nr);
popc(i,2).R=pop(i).R(1:t);  

r1=randi([1,MT-1]);
r2=randi([1,MT-r1]);

popc(i,1).SS=popc(i,1).S;
popc(i,2).SS=popc(i,2).S;

if Ns-1>2
r=randi([2 Ns-1]);
else
r=2;    
end

popc(i,1).S(1:r)=pop(f).S(1:r);
popc(i,2).S(1:r)=pop(k).S(1:r);

for h=1:2
popc(i,h).S=unique(popc(i,h).S);
if numel(popc(i,h).S)<Ns
   popc(i,h).S=(1:Ns); 
end
end

for h=1:2 
    
for k=1:t
    if k<t
         r(k)=randi([2,MT]);                   
         for j=1:MT+1
         if j<r(k) ||  j>r(k) 
            popc(i,h).Lsr(k,j)=popc(i,h).SS(1);
            popc(i,h).SS(1)=[];
         end
         if j==r(k)   
            popc(i,h).Lsr(k,j)=Ns+popc(i,h).R(k);
         end                 
         end
    end
    
    if k==t
         r(t)=randi([2,numel(popc(i,h).SS)+1]);
         for j=1:numel(popc(i,h).SS)+1
         if j<r(t) ||  j>r(t) 
            popc(i,h).Lsr(t,j)=popc(i,h).SS(1);
            popc(i,h).SS(1)=[];
         end
         if j==r(t)   
            popc(i,h).Lsr(t,j)=Ns+popc(i,h).R(t);
         end                 
         end    
    end
     
end
for k=1:t
    m(k)=0;
    for j=1:numel(popc(i,h).Lsr(k,:))
    if popc(i,h).Lsr(k,j)==0
       m(k)=m(k)+1;
    end
    end
end 

for k=1:t 
    popc(i,h).L(k,1)=popc(i,h).H; 
for j=1:numel(popc(i,h).Lsr(k,:)) 
    if popc(i,h).Lsr(k,j)~=0
    popc(i,h).L(k,j+1)=Nh+popc(i,h).Lsr(k,j);
    else
    popc(i,h).L(k,j+1)=0;    
    end 
    popc(i,h).L(k,numel(popc(i,h).Lsr(k,:))-m(k)+2)=popc(i,h).H; 
end
end

if Nv==5
    popc(i,h).V=randi(Nv-1,t,numel(popc(i,h).L(1,:))-1);
for k=1:t
  for j=1:numel(popc(i,h).Lsr(k,:))-m(k)  
    if dis(popc(i,h).L(k,j),popc(i,h).L(k,j+1))/Uv_star(popc(i,h).V(k,j))<=0.005
    popc(i,h).V(k,j)=5;    
    end
  end
end

else
    popc(i,h).V=randi(Nv,t,numel(popc(i,h).L(1,:))-1);    
end

for k=1:t
  for j=1:numel(popc(i,h).Lsr(k,:))-m(k)  
      if j==1 && popc(i,h).V(k,j)~=4
   popc(i,h).At(k,j)=DT+(dis(popc(i,h).L(k,j),popc(i,h).L(k,j+1))/Uv_star(popc(i,h).V(k,j))) + (wt*NL(popc(i,h).L(k,j),popc(i,h).L(k,j+1)));
   popc(i,h).TV(k,j)=popc(i,h).At(k,j)+UOT(popc(i,h).Lsr(k,j));  
      elseif j>1  && popc(i,h).V(k,j)~=4  
   popc(i,h).At(k,j)=popc(i,h).TV(k,j-1)+ (dis(popc(i,h).L(k,j),popc(i,h).L(k,j+1))/Uv_star(popc(i,h).V(k,j)))+ (wt*NL(popc(i,h).L(k,j),popc(i,h).L(k,j+1))); 
   popc(i,h).TV(k,j)=popc(i,h).At(k,j)+ UOT(popc(i,h).Lsr(k,j)); 
      elseif j==1 && popc(i,h).V(k,j)==4
   popc(i,h).At(k,j)=DT+(dis(popc(i,h).L(k,j),popc(i,h).L(k,j+1))/Uv_star(popc(i,h).V(k,j))) ;
   popc(i,h).TV(k,j)=popc(i,h).At(k,j)+UOT(popc(i,h).Lsr(k,j));  
      elseif j>1  && popc(i,h).V(k,j)==4  
   popc(i,h).At(k,j)=popc(i,h).TV(k,j-1)+ (dis(popc(i,h).L(k,j),popc(i,h).L(k,j+1))/Uv_star(popc(i,h).V(k,j))); 
   popc(i,h).TV(k,j)=popc(i,h).At(k,j)+ UOT(popc(i,h).Lsr(k,j)); 

      end
  end
end





    for k=1:t  
    for j=1:numel(popc(i,h).Lsr(k,:))-m(k) 
    if popc(i,h).At(k,j)<LTW(j) 
        if popc(i,h).TV(k,j)>UTW(j)|| popc(i,h).At(k,j)>UTW(j)  
           popc(i,h).L(:,:)=[];    
        end
    end
    end
    end
    
if ~isempty(popc(i,h).L(:,:))
for k=1:t
    for j=1:numel(popc(i,h).Lsr(k,:))-m(k)+1
    if popc(i,h).L(k,j)~=0 && popc(i,h).L(k,j+1)~=0     
       popc(i,h).TotalDD(k,j)= dis(popc(i,h).L(k,j),popc(i,h).L(k,j+1));
       popc(i,h).TotalTCC(k,j)=VC(popc(i,h).V(k,j))*popc(i,h).TotalDD(k,j)+FC(popc(i,h).V(k,j));
    end
    end
end
for k=1:t
   for j=1:numel(popc(i,h).L(k,:))
     popc(i,h).TotalD(k)=sum((popc(i,h).TotalDD(k,:))); 
     popc(i,h).TotalTC(k)=sum((popc(i,h).TotalTCC(k,:)));
   end
   popc(i,h).TotalD=sum(popc(i,h).TotalD(:));
   popc(i,h).TotalTC=sum(popc(i,h).TotalTC(:));
end

for k=1:t
popc(i,h).TotalRC(k)= URFC(popc(i,h).R(k));
popc(i,h).TotalUTR(k)=UT(popc(i,h).R(k));
end

popc(i,h).TotalRC=sum(popc(i,h).TotalRC);
popc(i,h).TotalUTR=sum(popc(i,h).TotalUTR);

popc(i,h).TotalC=t*URCH(popc(i,h).L(1,1))+ popc(i,h).TotalTC + popc(i,h).TotalRC+SCs;

if popc(i,h).TotalC<=B
popc(i,h).TotalC=popc(i,h).TotalC;
popc(i,h).TotalUT=t*UTH(popc(i,h).L(1))+popc(i,h).TotalUTR;
    for j=1:numel(popc(i,h).V)
        popc(i,h).TotalEm(j)=popc(i,h).TotalDD(j)*C(popc(i,h).V(j));
    end
popc(i,h).TotalEm=sum(popc(i,h).TotalEm);
else
popc(i,h).L=[]; 
end

if ~isempty(popc(i,h).L)
popc(i,h).Objectives=objFunction(popc(i,h));  
end

      
end
end

end

    popc=popc(:);
    
for i=nCrossover :-1:1
if isempty(popc(i).L)
popc(i)=[];
end
end
nCrossover=numel(popc);

for i=nCrossover :-1:1
if isempty(popc(i).Objectives)
popc(i)=[];
end
end
nCrossover=numel(popc);


for i=1:numel(popc)-1
    for j=numel(popc):-1:i
        if popc(i).TotalC==popc(j).TotalC
           popc(j)=[];
        end
    end
end

nCrossover=numel(pop);

%%%%

     %% Mutation
     popm=repmat(sol,nMutation,1);

for i=1:nMutation 
    
r=randi(npop);
popm(i).H=pop(r).H;
popm(i).S=pop(r).S;
popm(i).R=pop(i).R;

a=randi([2 Ns]);
b=randi([2 Ns]);
c=popm(i).S(b);
d=popm(i).S(a);

popm(i).S(a)=c;
popm(i).S(b)=d;
     r1=randi([1,MT-1]);
     r2=randi([1,MT-r1]);
     popm(i).SS=popm(i).S;

for k=1:t
    if k<t
         r(k)=randi([2,MT]);                   
         for j=1:MT+1
         if j<r(k) ||  j>r(k) 
            popm(i).Lsr(k,j)=popm(i).SS(1);
            popm(i).SS(1)=[];
         end
         if j==r(k)   
            popm(i).Lsr(k,j)=Ns+popm(i).R(k);
         end                 
         end
    end
    
    if k==t
         r(t)=randi([2,numel(popm(i).SS)+1]);
         for j=1:numel(popm(i).SS)+1
         if j<r(t) ||  j>r(t) 
            popm(i).Lsr(t,j)=popm(i).SS(1);
            popm(i).SS(1)=[];
         end
         if j==r(t)   
            popm(i).Lsr(t,j)=Ns+popm(i).R(t);
         end                 
         end    
    end
end
     
for k=1:t
    m(k)=0;
    for j=1:numel(popm(i).Lsr(k,:))
    if popm(i).Lsr(k,j)==0
       m(k)=m(k)+1;
    end
    end
end

for k=1:t 
    popm(i).L(k,1)=popm(i).H; 
for j=1:numel(popm(i).Lsr(k,:)) 
    if popm(i).Lsr(k,j)~=0
    popm(i).L(k,j+1)=Nh+popm(i).Lsr(k,j);
    else
    popm(i).L(k,j+1)=0;    
    end 
     popm(i).L(k,numel(popm(i).Lsr(k,:))-m(k)+2)=popm(i).H; 
end
end

if Nv==5
    popm(i).V=randi(Nv-1,t,numel(popm(i).L(1,:))-1);
for k=1:t
  for j=1:numel(popm(i).Lsr(k,:))-m(k)  
    if dis(popm(i).L(k,j),popm(i).L(k,j+1))/Uv_star(popm(i).V(k,j))<=0.005
    popm(i).V(k,j)=5;    
    end
  end
end

else
    popm(i).V=randi(Nv,t,numel(popm(i).L(1,:))-1);    
end
for k=1:t
  for j=1:numel(popm(i).Lsr(k,:))-m(k)  
      if j==1 && popm(i).V(k,j)~=4
   popm(i).At(k,j)=DT+(dis(popm(i).L(k,j),popm(i).L(k,j+1))/Uv_star(popm(i).V(k,j))) + (wt*NL(popm(i).L(k,j),popm(i).L(k,j+1)));
   popm(i).TV(k,j)=popm(i).At(k,j)+UOT(popm(i).Lsr(k,j));  
      elseif j>1  && popm(i).V(k,j)~=4  
   popm(i).At(k,j)=popm(i).TV(k,j-1)+ (dis(popm(i).L(k,j),popm(i).L(k,j+1))/Uv_star(popm(i).V(k,j)))+ (wt*NL(popm(i).L(k,j),popm(i).L(k,j+1))); 
   popm(i).TV(k,j)=popm(i).At(k,j)+ UOT(popm(i).Lsr(k,j)); 
      elseif j==1 && popm(i).V(k,j)==4
   popm(i).At(k,j)=DT+(dis(popm(i).L(k,j),popm(i).L(k,j+1))/Uv_star(popm(i).V(k,j))) ;
   popm(i).TV(k,j)=popm(i).At(k,j)+UOT(popm(i).Lsr(k,j));  
      elseif j>1  && popm(i).V(k,j)==4  
   popm(i).At(k,j)=popm(i).TV(k,j-1)+ (dis(popm(i).L(k,j),popm(i).L(k,j+1))/Uv_star(popm(i).V(k,j))); 
   popm(i).TV(k,j)=popm(i).At(k,j)+ UOT(popm(i).Lsr(k,j)); 
      end
  end
end

    for k=1:t  
    for j=1:numel(popm(i).Lsr(k,:))-m(k) 
    if popm(i).At(k,j)<LTW(j) 
        if popm(i).TV(k,j)>UTW(j) || popm(i).At(k,j)>UTW(j)  
           popm(i).L(:,:)=[];    
        end
    end
    end
    end
    
if ~isempty(popm(i).L(:,:))
for k=1:t
    for j=1:numel(popm(i).Lsr(k,:))-m(k)+1
    if popm(i).L(k,j)~=0 && popm(i).L(k,j+1)~=0     
       popm(i).TotalDD(k,j)= dis(popm(i).L(k,j),popm(i).L(k,j+1));
       popm(i).TotalTCC(k,j)=VC(popm(i).V(k,j))*popm(i).TotalDD(k,j)+FC(popm(i).V(k,j));
    end
    end
end
for k=1:t
   for j=1:numel(popm(i).L(k,:))
     popm(i).TotalD(k)=sum((popm(i).TotalDD(k,:))); 
     popm(i).TotalTC(k)=sum((popm(i).TotalTCC(k,:)));
   end
   popm(i).TotalD=sum(popm(i).TotalD(:));
   popm(i).TotalTC=sum(popm(i).TotalTC(:));
end

for k=1:t
popm(i).TotalRC(k)= URFC(popm(i).R(k));
popm(i).TotalUTR(k)=UT(popm(i).R(k));
end

popm(i).TotalRC=sum(popm(i).TotalRC);
popm(i).TotalUTR=sum(popm(i).TotalUTR);

popm(i).TotalC=t*URCH(popm(i).L(1,1))+ popm(i).TotalTC + popm(i).TotalRC+SCs;

if popm(i).TotalC<=B
popm(i).TotalC=popm(i).TotalC;
popm(i).TotalUT=t*UTH(popm(i).L(1))+popm(i).TotalUTR;
    for j=1:numel(popm(i).V)
        popm(i).TotalEm(j)=popm(i).TotalDD(j)*C(popm(i).V(j));
    end
popm(i).TotalEm=sum(popm(i).TotalEm);
else
popm(i).L=[]; 
end

if ~isempty(popm(i).L)
popm(i).Objectives=objFunction(popm(i));  
end

      
end

end
 popm=popm(:);

 
for i=nMutation :-1:1
if isempty(popm(i).L)
popm(i)=[];
end
end
nMutation=numel(popm);

for i=1:numel(popm)-1
    for j=numel(popm):-1:i
        if popm(i).TotalC==popm(j).TotalC
           popm(j)=[];
        end
    end
end

    pop=[pop
        popc
         popm];
     

     
    toc;
     
    % Non-Dominated Sorting
    [pop ,F]=NonDominatedSorting(pop);

    % Calculate Crowding Distance
    pop=CalcCrowdingDistance(pop,F);

    % Sort Population
    [pop ,F]=SortPopulation(pop);
    
    % Truncate
    pop=pop(1:npop);
    
    % Non-Dominated Sorting
    [pop, F]=NonDominatedSorting(pop);

    % Calculate Crowding Distance
    pop=CalcCrowdingDistance(pop,F);

    % Sort Population
    [pop, F]=SortPopulation(pop);
     
    % Store F1
     F1=pop(F{1});
    
    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Number of F1 Members = ' num2str(numel(F1))]);
    
    % Plot F1 Costs
     figure(1);
    PlotObjs(F1);

end
%% Results
pareto=pop(F{1},:);


for i=1:numel(pareto)-1
    for j=numel(pareto):-1:i+1
        if pareto(i).TotalC==pareto(j).TotalC
           pareto(j)=[];
        end
    end
end

