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
[wings]=wing;
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
%%
%DEBUG - make this better from the charts I made
function power = powerSelections(cells)
%powerSelection: Given number of battery cells return ideal power system
%from a lookup table.
%There is an ideal combination of battery, motor, propeller size, etc.. for
%each given battery voltage. This function takes a battery voltage input
%and returns all the parameters necessary to construct the ideal system.
%Future work could make ideal selection more automated.
%powerSystem: systems from 3 to 8 cells.
%            Motor Name,  KV, Propeller diameter (inches), Propeller pitch (inches),
%            Voltage, RPM, Current Draw (A), Power Consumption (W),
%            Flight time (minutes), thrust (pound force), Pitch speed
%            (ft/s), Efficiency (thrust/watt)*100, battery available watt
%            hours, estimated system weight (pounds).
powerSystem = zeros(15, 15);
%                'name' kv   pd  pp    v    rpm   amps    watts time        lbf  p-speed    effi      battW miss1,       miss2             
powerSystem(1,:) = [2	760	13	6.5	 11.1	6449  23.42	   260	12.46		3.591	58.2	1.381153846	54	1.123571035	1.123571035];
powerSystem(2,:) = [3	470	17	7	14.8	5972	35.03	518.4	11.13		7.194	58.1	1.387731481	96.2	2.097488987	2.097488987];	
powerSystem(3,:) = [4	1080	9	4.5	18.5	16703	51.78	958	11.58663883		6.592	104.4	0.688100209	185	1.817439427	3.308628855];
powerSystem(4,:) = [5	300	20	8	22.2	5632	45.22	1003.9	11.94142843		12.578	62.6	1.252913637	199.8	2.6769163	4.273832599];	
powerSystem(5,:) = [6	300	18	8	25.9	6724	40.92	1059.7	9.531943003		12.172	74.7	1.14862697	168.35	2.315682819	3.551365639];
powerSystem(6,:) = [7	300	16	8	29.6	7203	35.91	1063	10.85983067		11.495	80.1	1.081373471	192.4	2.280613987	3.688102974];
power=powerSystem(cells,:);
end
%%
%%%velocity solver via false position root finding
function [V,Drag,WeightT] = GenVelocityTest(Thrust,RPM,pitch,dp,CLw,WeightT,rho,AR,WingS,HStabS,VStabS,CDWing,Cw,CDHStab,Chstab,CDVStab,Cvstab,FSA,SAfT,Lf,Wf,SensorCalc,ARs,Ws)
%A function brought to you by Bryce Moran and Matthew Maddalon.
%don't panic

%List of Inputs
    %Thrust=Thrust.... (lb)
    %dp=diameter propeller (in)
    %CLw= coeff lift, wing (no dim)
    %WeightT=Weight Total (lb)
    %rho=density air (slug/ft^3)
    %AR=Wing aspect ratio (no dim)
    %WingS=wing reference area (ft^2)
    %HStabS=horizontal stabilizer reference area (ft^2)
    %VStabS=Vertical stabilizer reference area (ft^2)
    %CDwing= coeff drag wing (no dim)
    %Cw= chord wing (ft)
    %CDHstab= coeff horizontal stabilizer (no dim)
    %Chstab= chord hoirzontal stab (ft)
    %CDVstab= coeff drag vertical stabilizer (no dim)
    %Cvstab=chord vertical stab (ft)
    %FSA=fuselage frontal surface area (ft^2)
    %SAfT= Surface Area fuselage Total (ft^2)
    %Lf=length fuselage (ft)
    %Wf= width fuselage (in)
    %SensorCalc = logical value, 1 or 0. 0 will not run sensor calculations
    %ARs=sensor aspect ratio
    %Ws = weight sensor

%Constants
    K=1/(pi*0.8*AR); %e=0.8. K is an aero value. e is efficiency rating for the wing (rectangular)
    mu=0.3776e-6; %for 19 Centigrade %viscous friction coefficient. calculated for a certain temperature
    muk=1.612e-4; %kinmatic viscosity (ish) related to mu
    Ds=1.5; %Diameter sensor (in)

    
    %landing gear should not be calculated in drag solver
    %below is landing gear solver %DEBUG
    
    Cd_Al=0.05;         % Let Cd = 0.05 for the flatbar
    Cd_wheel = 0.20;     % Let Cd = 0.245 for the wheels
    top = Wf;         %length of top flat part (change to fuse width)
    t = 1/8;            %thickness of flatbar
    theta = 38*pi/180;  %Angle of inside of bend (minus 90)
    gearwidth = 1;      %width of flatbar
    width_wheel = 0.5;    %width of wheels
    radius_wheel = 1.5; %wheel radius
    rho_Al = 0.001122368;         %density of Al in slug/in^3
    rho_air = rho/(12^3);  %density of air at sea level in slug/in^3
        
    %ldg gear information
    %Courtesy of Matt Lutton
    
    height = dp/2+0.5; %dp is prop diamter - height of prop centerline plus tolerance
    
    % Profile area of aluminum pieces
    
    area_Al= (top + 2*(height/cos(theta)+height))*t;
    volume_Al = area_Al*gearwidth;
    weight_Al = rho_Al*volume_Al;
    area_wheel = width_wheel*radius_wheel*3*2;    %profile area for 3 wheels
    WeightT=WeightT+weight_Al;
    %%%end landing gear drag
    
    
    
%Induced

    Induced= @(v)(2*K*WeightT^2)/(rho*WingS*v^2); %It's the induced Drag. doesn't account for lift from tail.

%ldg gear parasitic drag

    drag_Al =@(v) Cd_Al*rho_air*(v*12)^2*area_Al*.5; %aluminum bar
    drag_wheel = @(v) Cd_wheel*rho_air*(v*12)^2*2*area_wheel*.5;
    ldgGearP=@(v) drag_Al(v)+drag_wheel(v); %p for parasitic

    %Fuselage Parasitic 
    
    CDf=0.75; %good round estimate
    
%Parasitic drag= sum (wing + HStab + VStab + ldg gear + Fuselage

    Parasitic= @(v) (0.5*rho*(v^2)*WingS*CDWing) + (0.5*rho*(v^2)*HStabS*CDHStab)  + (0.5*rho*(v^2)*FSA*CDf)+ (0.5*rho*(v^2)*VStabS*CDVStab)+ldgGearP(v); %its the parasitic drag - from frontal cross section

%coeff skin friction. everything has a unique one

    Temp=273+19; %kelvin - average april weather at competition site.
    R=287; %gas consant - no need to change year to year
    gamma=1.4; %air constant - no change year-to-year
    a=3.28*sqrt(gamma*R*Temp);%with metric to imperial conversion. a is speed of sound. 3.28 is conversion factor to feet
    Cfl= @(v,L) 1.328/sqrt(v*L/muk);%coeff friction laminar
    Cft=@(v,L) (0.455)/((log10(v*L/muk)^2.58)*(1+0.144*(v/a)^2)^0.65); %coeff friction turbulent
    
    Cf=@(v,L,k) k*Cfl(v,L)+(1-k)*Cft(v,L); %trueish - some ratio of laminar and turbulent

%k for not fuselage is about 0.1
%k (a constant) for fuselage should be 0 (all turbulent)
%k, the constant you have no idea about, is the transition point along the
%length from laminar to tubulent flow
%
%-Bryce

    %ldg gear skin friction
    
    refAl=2*(height/cos(theta)+height)*gearwidth/144;
    WheelWet=(2*pi*(radius_wheel)^2+ pi*2*radius_wheel*width_wheel)/144;

%Skin friction Drag( sum: wing, HStab, VStab+...
%drag equation except using skin friction instead of parasitic
Skin= @(v) 2*(0.5*rho*(v^2)*Cf(v,Cw,1)*WingS) +2*(0.5*rho*(v^2)*Cf(v,Chstab,1)*HStabS) + 2*(0.5*rho*(v^2)*Cf(v,Cvstab,1)*VStabS) + (0.5*rho*(v^2)*Cf(v,Lf,0.9)*SAfT) + 2*(0.5*rho*(v^2)*Cf(v,(2/12) ,0)*2*(height/cos(theta)+height)*t/144) + 3*(0.5*rho*(v^2)*Cf(v,2*radius_wheel/12,0)*WheelWet/144);

%prop drag - accounts for accelerated flow across the rest of the aircraft
%all this for like 0.001 change :/  
%-Bryce

    Ap=pi*(dp/2)^2;
    Sl=SAfT + (dp-Wf)*12*Cw;%total wet fuselage Surface Area in prop wash
    q=@(v) 0.5*rho*v^2;% is q
    Dcl=@(v) 0.5*(Cw/dp)*CLw* Thrust/(q(v) * Ap);%change in CL
    dci=@(v) 2*Dcl(v)*CLw/(pi*AR);%change in Cdi
    Vs=@(v) (v+Thrust/(rho*Ap*v));%velocity of prop slipstream
        SAT=SAfT + 2*(dp-Wf)*12*Cw;
        kwing=(2*(dp-Wf)*12*Cw)/SAT;
        kfuse=SAfT/SAT;
    CFp=@(v) kwing*Cf(Vs(v),Cw,0.1) + kfuse*(Cf(Vs(v),Lf,0));  
    DC0=@(v) Cf(v,Lf,0)*(Sl/WingS)*( ((Vs(v))/v)^2 -1);%change in Cd0
    Prop=@(v) (DC0(v)+dci(v)) *0.5*rho*(Vs(v))^2;%chnage in drag (additional drag)

%Sensor drag
%DEBUG - delete this - we're not towing anythign anymore
if SensorCalc>0
    rs=Ds/2;
    rss=rs/12;
    SAs=4*pi*( ((rs^2)^1.6 + 2*(rs*Ds*ARs/2)^1.6)/3)^(1/1.6); %surface area sensor
    SAs=SAs/144; %ft^2
    CDs=@(ARs) 0.35-0.025*ARs;
    sensor= @(v) (0.5*rho*(v^2)*CDs(ARs)*pi*(rss^2)) + (0.5*rho*(v^2)*SAs*Cf(v,Ds*ARs/12,1));
    WeightT=WeightT+Ws;
else
    sensor=@(v) v*0;
end
    
%Total Drag
Drag= @(v)   Induced(v)+Parasitic(v) +Skin(v) + Prop(v) + sensor(v); %DEBUG

PS=pitch/12 * 5614 /60; %propeller speed

PD=pitch/dp; %this is dynamic thrust function. Article this is based off of is in resources google drive somewhere
FS=@(v) v;
Ts=Thrust;
s=0;
if PD<0.6
    T= @(v) Ts*(1 - FS(v)/(PS*(PD+0.2)/PD));
elseif (s/PS)*PD<(PD-0.6)
    T=@(v) Ts;
else
    T=@(v) Ts*(1- ((FS(v)*PD/PS)-(PD-0.6))/0.8);
end

func=@(v) (Drag(v)-T(v)) + FS(v) - FS(v);
%end of dynamic thrust

xl=-30; xu=1000; ea=1e-4; %false position root finder #1 to use initially, try the other when this breaks. V has to be between xl and xu
Error=1;
v=xl;
xrold=0;
while Error>ea/100
 xrold=v;
 v=(xl+xu)/2;
 check=(func(xl)*func(v));
 if check <0
 xu=v;
 elseif check >0
 xl=v;
 end
 s=FS(v);
 Error=abs((v-xrold)/v);
end
V=v;



%Emergency, if stuff breaks, its because of induced drag. Solver reverts to a
%different form of induced drag with the assumption that L=W. Without this
%it returns complex numbers :(. Happens when velocity becomes so small it's
%computationally zero. 
if V<0
    V=1; %DEBUG - make an error flag - functions look for exactly this value for now
elseif V>900 %no way we ever get this fast.
    V=1;
end

if V==1; %new induced drag if the other broke
    Induced=@(v) WeightT*(CLw)*K*WingS;
end

Drag= @(v)   Induced(v)+Parasitic(v) +Skin(v) + Prop(v);
func=@(v) (Drag(v)-T(v));

xl=-30; xu=1000; ea=1e-4; %this is root finder again but with different induced method to avoid problem from above if we get here
Error=1;
xr=xl;
xrold=0;
while Error>ea/100
 xrold=xr;
 xr=(xl+xu)/2;
 check=(func(xl)*func(xr));
 if check <0
 xu=xr;
 elseif check >0
 xl=xr;
 end
 Error=abs((xr-xrold)/xr);
end

V=xr;
if V<0
    V=1;
elseif V>900
    V=1; %DEBUG -- v == 1 is still the it's broke flag
end
%DEBUG - add broken down drag values here - return an array for sanity
%check purposes.
%DEBUG - look into describing drag in terms of parts instead of classes
%(wing instead of skin alone)
D=Drag(V);
% iter
end
%%
function [FSA,fh,fw,fl,Weightfuse,SAfT] = SaestsNew(AR,nc) %surface area estimate for fuselage
%DEBUG - need new version of this - just a look up table. how will fuselage
%grow
cw = 1.5; %width sensor
ch = 1.5; %height sensor
cl=cw*AR; %length sensor
fh = 0; %fuse height
fw = 0; %fuse width
fl = 0; %fuse length
rhoCarbon=120.486; %lb/ft^3
if nc<=4
    fh=cw;
    fw=cw;
    fl=nc*cw*AR;
elseif nc<=8
    fh=2*cw;
    fw=cw;
    fl=4*cl;
elseif nc<=16
    fh=2*cw;
    fw=2*cw;
    fl=4*cl;
elseif nc<=24
    fh=3*cw;
    fw=2*cw;
    fl=4*cl;
elseif nc<= 36
    fh=3*cw;
    fw=3*cw;
    fl=4*cl;
elseif nc<=48
   fh=4*cw;
   fw=3*cw;
   fl=4*cl;
elseif nc<=64
    fh=4*cw;
    fw=4*cw;
    fl=4*cl;
end

SAbrick=2*fl*fh+2*fl*fw;
Hf=sqrt( (0.5*fh)^2 +6^2);
SAfTri=2*Hf*fw + fh*6;
Hb=sqrt( fh^2 +8^2);
SAbTri=fw*8 + fw*Hb + fh*8;
SAfT=SAbrick+SAfTri+ SAbTri;
SAfT=SAfT/144;%Surface Area fuse Total,ft^2
FSA=fw*fh/144;%Frontal Surface Area,ft^2
Vol=SAfT*(1/16)/12;%The volume in ft^3
Weightfuse=Vol *rhoCarbon;
end
%%

%this funcion uses the monoplane equation (matrix method) for generating a
%wing

%the function creates a wing for each airfoil and aspect ratio to determine
%the approximate C_l, C_d, and weight for a wing with that configuration.
%The method assumes a rectangular wing (chord length is the same at root
%and tip).
%The code uses the Prandtl-Lifting-Line theory to estimate the general lift
%distribution for the wing.
% clc;clear;
function [wings]=wing %DEBUG - all the orange comes from this function

b=5;                        %span of the wing %DEBUG if wingspan is changing this needs to input
rho_spar=68.3;              %the density of the carbon fiber spar in the wing (lb/ft^3)
Ax_spar=pi*(.1875/12)^2;    %the cross-section of the carbon fiber spar in the wing
V_spar=b*Ax_spar;           %the volume of the wing spar
W_spar=V_spar*rho_spar;     %the eight of the wing spar

rho_foam=2.13;              %density of the foam used to make the wing (lb/ft^3)

%the next set of code creates a reference matrix for each airfoil, all data
%taken at alpha=3 degrees: %DEBUG - talk to jasmin about the alpha
%NACA - airfoil name
%cl - coefficient of lift for that airfoil
%cd - coefficient of drag for that airfoil
%alpha0 - zero lift angle of attack for that airfoil
%max_t - the maximum thickness of that airfoil 
%max_tx - the position of maximum thickness for that airfoil
%data - the complete matrix of airfoil information

%63210 - 63A210
%8036 - s8036
%8037 - s8037
%2091 - s2091
%3002 - s3002
%67 - e67
%66 - e66
%0 - clarkY
%7062 - sd7062
%3012 - hq 3.0/12
NACA=[1408 1410 1412 2408 2410 2412 2415 2418 2421 2424 4412 4415 4418 4421 6409 6412 63210 8036 8037 2091 3002 67 66 0 7062 3012];
cl=[.3377 .4082 .4472 .4445 .4688 .4528 .4065 .3356 .2726 .0542 .5279 .3020 .2943 .0590 .5997 .5632 .3963 .1968 .2186 .6095 .5698 .3545 .4302 .5275 .5115 .3807];
cd=[.02271 .02562 .02791 .02628 .02950 .03304 .03941 .04769 .06115 .08159 .04426 .05278 .06234 .07671 .04926 .05831 .03095 .04811 .05268 .03799 .0342 .04189 .03963 .03775 .04793 .04001];
alpha0=(pi/180)*[-.25 .125 .9 -.75 0 .25 -.25 -.125 .875 2.75 -1.25 -.25 -.25 1.375 -1.25 -2 -3 -1.625 -2.125 -1.25 -.5 .25 -.525 -.2805 -1.75 .625];
max_t=[.08 .1 .12 .08 .1 .12 .15 .18 .21 .24 .12 .15 .18 .21 .09 .12 .1 .16 .16 .101 .099 .116 .101 .117 .14 .12];
max_tx=[.3 .299 .299 .299 .299 .3 .3 .3 .3 .297 .3 .309 .3 .3 .293 .301 .349 .368 .335 .26 .306 .326 .287 .28 .255 .35];
data=[NACA' cl' cd' alpha0' max_t' max_tx'];


alphac=3*pi/180;        %the cruise angle of attack of the airfoil (3 degrees)
alphat=15*pi/180;       %the takeoff angle of attack (10 degrees)
alphae=alphac-alpha0;   %the effective angle of attack for the airfoil (cruise angle minus zero angle)
alphamax=alphat-alpha0; %the effective angle of attack at takeoff (takeoff angle minus 0 angle)

AR=8:.5:15;             %the wing aspect ratios being considered 


%DEBUG - make this an input instead of internal so AR's aren't hard-coded

theta=[pi/4 pi/2 3*pi/4 pi];    %the theta values for positions down the wing
n=[1 3 5 7];                    %n values for determining the wing's coefficient of drag

%psi and zeta are an arbitrary matrix and vector, respectively, to solve
%for the coefficients, the A vector, needed to determine the lift and drag
%coefficients for the wing
psi=zeros(4,4);
zeta=zeros(length(theta),1);
zetamax=zeta;

wings=zeros(length(NACA),3,length(AR));     %matrix to store wing data

%for loops to run the wing analysis for each aspect ratio and airfoil
for x=1:length(AR)
    for y=1:length(alphae)
    mu=pi/(2*AR(x));
        for z=1:length(theta)
            psi(z,:)=[sin(theta(z))*(mu+sin(theta(z))) , sin(3*theta(z))*(3*mu+sin(theta(z))) , sin(5*theta(z))*(5*mu+sin(theta(z))) , sin(7*theta(z))*(7*mu+sin(theta(z)))];
            zeta(z)=mu*alphae(y)*sin(theta(z));
            zetamax(z)=mu*alphamax(y)*sin(theta(z));
        end
        A=psi\zeta;                     %the coefficients used to determine C_l and C_di for the wing
        Amax=psi\zetamax;               %the coefficients used to determine C_lmax (takeoff)
        C_l=A(1)*pi*AR(x);              %C_l for the wing
        C_lmax=Amax(1)*pi*AR(x);      %C_lmax (C_l at takeoff) for the wing
        C_di=pi*AR(x)*dot(n,A.^2);      %C_di for the wing
        c=b/AR(x);                      %chord length of the wing (based on aspect ratio)
        t=c*max_t(y);                   %thickness of the wing
        
        %the next set of code determines the weight of the wing, by
        %approximating the cross-sectional area as 2 right triangles, with
        %the position of maximum thickness for the airfoil used to define
        %the bases of the triangles.
        ax=(c*max_tx(y)*t/2)+(c*(1-max_tx(y))*t/2)-Ax_spar;     %wing cross-sectional area, minus the area taken up by the spar 
        W_wing=(ax*b)*rho_foam;                                 %weight of the foam used for the wing
        W=W_spar+W_wing;                                        %weight of the wing and spar
        
        %the next 4 lines assign the wing data to positions in the wing
        %matrix. The first column is the wing C_l, the second is the wing
        %C_lmax, the third C_di, and fourth weight. Each row corresponds to 
        %an airfoil in the NACA vector, and each page corresponds to an 
        %aspect ratio. E.g., the first airfoil in the NACA vector is the 
        %1408, and the first AR is 8. So the first row of the first page in 
        %'wings' is the C_l, C_di, and weight for a wing with a NACA 1408 
        %airfoil, and an aspect ratio of 8.
        %The fifth column is the name of the airfoil (divided by 1000 so
        %that the other data is visible), for reference.
        wings(y,1,x)=C_l;                                       
        wings(y,2,x)=C_lmax;
        wings(y,3,x)=C_di;
        wings(y,4,x)=W;
        wings(:,5,x)=NACA'./1000; 
    end
end
end
%%
function D= TakeoffChecker(Thrust,W,rho,WingS,AR,Cwing,CLm,CD0,dp,fh,RPM,pitch)
%DEBG - this function might not actually work. gets an idea but not sure of
%accuracy. use at your own risk.
%Inputs
    %T=thrust
    %W=Mission 2 weight
    %rho=density air (slug/ft^3)
    %WingS=reference area wing (ft^2)
    %AR=aspect ratio, wing
    %Cwing=chord wing (ft)
    %CLm=coeff lift, max
    %CD0=wing cd0 @ chosen alpha
    %dp=prop diameter (in)
    %fh=fuselage height (in)
    %RPM =Motor RPM, max
    %Pitch= propeller pitch
%constants
    g=32.2; %acceleration due to gravity (ft/s/s)
    mu=0.008; %coeff rolling friction
    e=0.8; %efficiency of wing

K=1/(pi*AR*e);    %Wing K constant
h=(dp/2+0.5+fh)/12; %Height of wing above ground (ft)
phi= 1 -(2*e/pi)*log(1+ ( (pi*Cwing)/(8*h))^2); %scaling constant
Kg=phi*K;%K constant, ground effect
Clg=mu/(2*Kg); %coeff lift, ground effect
Kw=1.4817*(5.81e-5);%Intermediate constant
dCd0lg=(W/WingS)*(Kw)*(W/g); %change in CD0 due to ground effect
CD0lg=CD0-dCd0lg;%CD0 ground effect
CDg=CD0lg+Kg*Clg^2;

Ap=(pi*0.25*dp^2)/144; %Prop area(ft^2)
Vr=1.2*sqrt(2*W/(rho*CLm*WingS)); %Rotation speed
%dynamic thrust equation
PS=pitch/12 * 5614 /60;
PD=pitch/dp;
FS=@(v) v;
Ts=Thrust;
if PD<0.6
    T= @(v) Ts*(1 - FS(v)/(PS*(PD+0.2)/PD));
elseif (FS(v)/PS)*PD<(PD-0.6)
    T=@(v) Ts;
else
    T=@(v) Ts*(1- ((FS(v)*PD/PS)-(PD-0.6))/0.8);
end
A=g*(T(Vr)/W - mu); %A constant
%B constant. absolutely ignore that little a or nothing works
B=(g/W)*(0.5*rho*WingS*(CDg-mu*Clg) );% +a);
D=(1/(2*B))*log(A/(A-B*Vr^2)); %ground roll (ft)
end

%%
function [HS_area, HS_chord, HS_span, HS_weight, HS_Cd, VS_area, VS_chord, VS_span, VS_weight, VS_Cd] = empennage(wingarea, HSAR, VSAR)
    %DEBUG - make sure this is being output in correct order. This does not
    %solve via the monoplane equation. probably should change to monoplane
    %for higher accuracy even if it's a bit computationally slower
    rho_XPS = 2.3;   %lb/ft^3
    e = 0.7;  %efficiency factor
    alpha = 3.*pi/180;  %assume 3 degrees
    Cd0 = 0.0165;
    HS_area = wingarea/4;
    HS_chord = sqrt(HS_area/HSAR);
    HS_span = HS_area/HS_chord;
    airfoil = @(x) .6*(0.2969.*sqrt(x)-0.1260.*x-0.3516.*x.^2+0.2843.*x.^3-0.1015.*x.^4);  %NACA0012 airfoil is just assumed
    CSH_area = 2*integral(airfoil, 0, HS_chord); %cross sectional area of airfoil horizontal stab
    HS_weight = CSH_area*HS_span*rho_XPS;   %lb
    HS_Cdi = 4.*pi.*alpha.^2/(HSAR.*e);  %induced drag equation (NASA)
    HS_Cd = Cd0+HS_Cdi;
    VS_Cd = Cd0;     % no alpha, no Cdi
    VS_area = HS_area/2;
    VS_chord = sqrt(VS_area/VSAR);
    VS_span = VS_area/VS_chord;
    CSV_area = 2*integral(airfoil,0,VS_chord);   % cross sectional area of airfoil vertical stab
    VS_weight = CSV_area*VS_span*rho_XPS;    %lb
end