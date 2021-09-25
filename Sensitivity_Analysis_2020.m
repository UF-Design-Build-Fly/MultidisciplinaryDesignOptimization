clc;clear;
format shortE

%List of constants
g=32.2; %gravitional accel (%ft/s/s)
rho=0.00235308; %Desnity air at Tucson, Az (slug/ft^3)
HSAR=4.5;%Horizontal Stab Aspect Ratio %DEBUG - get these from jasmin
VSAR=5;%Vertical Stab Aspect Ratio
Dlap=3000; %Lap distance in ft
Ds=1.5; %diameter sensor (in)
numV=51;%Number of "regular" variables
gmin=64; %ground mission minimum time. 64 + num containers

%main code
[wings]=wing();
Span=5;
wings; %DEBUG - dont need it
AR=(4:0.2:15);
[a,b]=size(wings(:,:,1));
wingbits=zeros(1,2); %stores info about the airplane
tailbits=zeros(1,2);
powerbits=zeros(a,15);
ARs=(4:0.5:10); %sensor aspect ratio
Ws=(0.25:0.25:1); %weight sensor
workss=zeros(a,numV); %a is the size of "wings". variable matrix
works=zeros(a,15); %accounting for power system
worksss=zeros(a,numV,length(AR),6,length(ARs),length(Ws)); %6d matrix of aircraft variables and information

[A,B,C,D,E,F]=size(worksss);


%Explaining the matrix worksss (I just kept add an s,idk, it just works...)
%Position 1=airfoil  %DEBUG - position is dimension
%position 2= regular variables (CLWing, weight, etc)
%google doc of all variables - vector variable list: https://docs.google.com/spreadsheets/d/1FmU5eVV3CRve9tXYGc6SgpZ6opuZdF7lGGSVIrlqlDo/edit#gid=0 
%position 3= Aspect ratio of wing
%position 4=propulsion system;
%position 5= aspect ratio sensor
%position6 = weight sensor


Winner=zeros(a,numV,length(AR),6,length(ARs),length(Ws)); %equivalent to the worksss matrix for the best aircraft

bestm=zeros(64,numV); %DEBUG -figure out what this does


%DEBUG - variable debug determines how iterations to run
for x=1:1 %length(AR)
        %just building an airplane from here on out
        Cwing=Span/AR(x);
        WingS=Span*Cwing;
        wingbits=[WingS Cwing];
        work=wings(:,:,x);
        
        [HS_area, HS_chord, HS_weight, HS_Cd, VS_area, VS_chord, VS_weight, VS_Cd] = empennage(WingS, HSAR, VSAR);
        
        tailbits=[HS_area, HS_chord, HS_weight, HS_Cd, VS_area, VS_chord, VS_weight, VS_Cd];
        
        for airfoil=1:a
            works(airfoil,:)=[work(airfoil,:) wingbits tailbits];
        end
    
        for Cell=1:6
            
            %power stuff...
            
            power=powerSelections(Cell);
            
                 for i=1:a
                   powerbits(i,:)=power;
                 end
                 
            workss=[works powerbits ];
            [aa,bb]=size(workss);
            
            filler=zeros(aa,B-bb);
            workplus=[workss filler];
            worksss(:,:,x,Cell,1)=workplus;
            size(worksss);

             for h=1:length(ARs)
                 worksss(:,:,x,Cell,h)=worksss(:,:,x,Cell,1);
             end

            for g=1:length(Ws)
                worksss(:,:,x,Cell,:,g)=worksss(:,:,x,Cell,:,1);
            end

            %M3 start
            
            for airfoil=1:a

             n=airfoil; %n=5 for testing, NACA 2410, else set equal to wing (n=wing)

                for i=1:13

                    ARs(i);

                    for j=1:length(Ws)

                    %M2 start

                        for cargo=1:64 %higher? 

                            [FSA,fh,Wf,Lf,W,SAfT] = SaestsNew(ARs(i),cargo);

                            Wwing=worksss(n,4,x,Cell,i,j);
                            Whs=worksss(n,10,x,Cell,i,j);
                            Wvs=worksss(n,14,x,Cell,i,j);
                            Wm2=worksss(n,29,x,Cell,i,j);

                            WeightT= Wwing + Whs + Wvs + Wm2+ W;%
                            worksss(n,31,x,Cell,i,j)=WeightT;%empty weight
                             WeightT=WeightT+Ws(j)*cargo;
                            
                            Thrust=worksss(n,25,x,Cell,i,j);
                            RPM=worksss(n,21,x,Cell,i,j);
                            pitch=worksss(n,19,x,Cell,i,j);
                            dp=worksss(n,18,x,Cell,i,j);
                            CLw=worksss(n,1,x,Cell,i,j);
                            WingS=worksss(n,6,x,Cell,i,j);
                            HStabS=worksss(n,8,x,Cell,i,j);
                            VStabS=worksss(n,12,x,Cell,i,j);
                            CDWing=worksss(n,3,x,Cell,i,j);
                            Cw=worksss(n,7,x,Cell,i,j);
                            CDHStab=worksss(n,11,x,Cell,i,j);
                            Chstab=worksss(n,9,x,Cell,i,j);
                            CDVStab=worksss(n,15,x,Cell,i,j);
                            Cvstab=worksss(n,13,x,Cell,i,j);

                          [V,Drag,Weight] = GenVelocityTest(Thrust,RPM,pitch,dp,CLw,WeightT,rho,AR(x),WingS,HStabS,VStabS,CDWing,Cw,CDHStab,Chstab,CDVStab,Cvstab,FSA,SAfT,Lf,Wf,0,ARs(i),Ws(j));
                          
                             if (0.5*rho*CLw*WingS*V^2)< Weight
                                 V=1;
                             end
                             
                             Velocity=V;

                            %if passes lift check, takeoff check

                            CLm= worksss(n,2,x,Cell,i,j);
                            Weight;
                            Dtakeoff=0;

                            if V~=1
                                Dtakeoff=TakeoffChecker(Thrust,Weight,rho,WingS,AR(x),Cwing,CLm,CDWing,dp,fh,RPM,pitch);
                            end

                            if Dtakeoff>95
                                V=1;
                            end

                            %find M2 Raw score
                            worksss(n,32,x,Cell,i,j)=cargo;
                            worksss(n,33,x,Cell,i,j)=Weight;%m2 weight
                            worksss(n,34,x,Cell,i,j)=V;
                            V;
                            worksss(n,35,x,Cell,i,j)=Dtakeoff;
                            disp(['  AR ','  Airfoil','  Cell ',' ARs ','  Ws ','  Cargo'])
                            disp([x n Cell i j cargo]) %Lets you keep track of where calcs are at

                            timem2=3*Dlap/V;
                            
                            if timem2>(5*60)
                                timem2=1000000000;
                            end
                            
                            M2Raw=cargo/timem2;
                            worksss(n,36,x,Cell,i,j)=timem2;
                            worksss(n,37,x,Cell,i,j)=M2Raw;
                            bestm(cargo,:)=worksss(n,:,x,Cell,i,j);
                        
                        end
                        
                     %M3 weight
                 %M3 Raw score
                 %find best M2 score (normalize all data to that)
                 %find best M3 score (normalize all data to that)
                 %find best sum M2+M3
                 %record best score for given cargo containers
                 
                      bestm;
                      m2max=0;
                      m3max=0;
                     
                      WeightM3=0;
                      VM3=0;
                      timeM3=0;
                      M3raw=0;
                      M3norm=0;
                      M=0;
                      
                      for k=1:64
                          
                     WeightM3=bestm(k,30)+bestm(k,31);
                     [V,Drag,Weight] = GenVelocityTest(Thrust,RPM,pitch,dp,CLw,WeightM3,rho,AR(x),WingS,HStabS,VStabS,CDWing,Cw,CDHStab,Chstab,CDVStab,Cvstab,FSA,SAfT,Lf,Wf,1,ARs(i),Ws(j));
                     
                     bestm(k,39)=Weight;
                     bestm(k,40)=V;
                     bestm(k,41)=Dlap/V;%time perlap M3
                     bestm(k,42)=floor(V*10*60/Dlap);%num laps
                     bestm(k,43)=bestm(k,42)*Ds*Ws(j);
                     bestm(k,45)=54+10.*bestm(k,32);
                      end
                      
                      m2max=max(bestm(:,37));
                      m3max=max(bestm(:,43));
                      
                      bestm(:,38)=bestm(:,37)./m2max; %m2 norm
                      bestm(:,44)=bestm(:,43)./m3max; %m3 norm
                      bestm(:,46)=gmin./bestm(:,45);
                      bestm(:,47)=bestm(:,38) + bestm(:,44) + bestm(:,46);
                      
                      [M,I]=max(bestm(:,47));
                      bestm;
                      best=bestm(I,:);
                      
                      Winner(n,:,x,Cell,i,j)=best;
                            
                    end
            %find best score for given weight, keep that
                end
            
             end
        %find best score for given aspect ratio, keep that 
        end
end

bestwin=zeros(1,numV);
bestwin2=bestwin;
bestwin3=bestwin;
bestwin4=bestwin;
bestwin5=bestwin;

bestint=zeros(4,numV);%intermediate matrix

bestint2=zeros(length(ARs),numV);

bestin3=zeros(6,numV);

bestint4=zeros(length(AR),numV);

bestint5=zeros(a,numV);

%Search winning solutions

Winners=zeros(a,numV,length(AR),6,length(ARs));
Winnerss=zeros(a,numV,length(AR),6);
Winnersss=zeros(a,numV,length(AR));
Winnerssss=zeros(a,numV);
WinnerWinner=zeros(1,numV);%chicken dinner....

for airfoil=1:a
    
    for x=1:length(AR)
        
        for c=1:6
            
            for i=1:length(ARs)

                for j=1:length(Ws)
                    
                    bestint(j,:)=Winner(airfoil,:,x,c,i,j);
                                        
                end
                
                maxm2=max(bestint(:,37));
                maxm3=max(bestint(:,43));
                
                
                bestint(:,38)=bestint(:,37)/maxm2;
                bestint(:,44)=bestint(:,43)/maxm3;
                bestint(:,46)=gmin./bestint(:,45);
                bestint(:,47)=bestint(:,38)+bestint(:,44)+bestint(:,46);
                [M,I]=max(bestint(:,47));
                
                bestwin=bestint(I,:);
                bestwin(51)=Ws(I);
                bestwin(50)=ARs(i);
                bestwin(49)=c;
                bestwin(48)=AR(x);
                bestwin;
                Winners(airfoil,:,x,c,i)=bestwin; 
                
            end
            
            for i=1:length(ARs)
                
                bestint2(i,:)=Winners(airfoil,:,x,c,i);
                
            end
            
            maxm2=max(bestint2(:,37));
            maxm3=max(bestint2(:,43));
            
            bestint2(:,38)=bestint2(:,37)./maxm2;
            bestint2(:,44)=bestint2(:,43)./maxm3;
            bestint2(:,46)=gmin./bestint2(:,45);
            bestint2(:,47)=bestint2(:,38)+bestint2(:,44)+bestint2(:,46);
            
            temp=bestint2(:,47);
            [M,I]=max(temp(~isinf(temp)));
            
            bestint2;
            bestwin2=bestint2(I,:);
            Winnerss(airfoil,:,x,c)=bestwin2; 
            
        end
        
        for c=1:6
            bestint3(c,:)=Winnerss(airfoil,:,x,c);
        end
                    
        maxm2=max(bestint3(:,37));
        maxm3=max(bestint3(:,43));

        
        bestint3(:,38)=bestint3(:,37)./maxm2;
        bestint3(:,44)=bestint3(:,43)./maxm3;
        bestint3(:,46)=gmin./bestint3(:,45);
        
        bestint3(:,47)=bestint3(:,38)+bestint3(:,44)+bestint3(:,46);
        
            temp=bestint3(:,47);
            [M,I]=max(temp(~isinf(temp)));
                
        bestwin3=bestint3(I,:);
        Winnersss(airfoil,:,x)=bestwin3; 
    end
    
    for x=1:length(AR)
        
        bestint4(x,:)=Winnersss(airfoil,:,x);
        
    end
    maxm2=max(bestint4(:,37));
    maxm3=max(bestint4(:,43));

                
    bestint4(:,38)=bestint4(:,37)./maxm2;
    bestint4(:,44)=bestint4(:,43)./maxm3;
    bestint4(:,46)=gmin./bestint4(:,45);
    bestint4(:,47)=bestint4(:,38)+bestint4(:,44)+bestint4(:,46) ;
            temp=bestint4(:,47);
            [M,I]=max(temp(~isinf(temp)));

                
    bestwin4=bestint4(I,:);
    Winnerssss(airfoil,:)=bestwin4; 
    
end

maxm2=max(Winnerssss(:,37));
maxm3=max(Winnerssss(:,43));

Winnerssss(:,38)=Winnerssss(:,37)./maxm2;
Winnerssss(:,44)=Winnerssss(:,43)./maxm3;
Winnerssss(:,46)=gmin./Winnerssss(:,45);

Winnerssss(:,47)=Winnerssss(:,38)+Winnerssss(:,44)+Winnerssss(:,46);
            temp=Winnerssss(:,47);
            [M,I]=max(temp(~isinf(temp)));
                
WinnerWinner=Winnerssss(I,:);

disp(WinnerWinner);
