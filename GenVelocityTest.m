%%
%%%velocity solver via false position root finding
function [V,Drag] = GenVelocityTest(G, Thrust,RPM,pitch,dp,CLw,WeightT,rho,AR,WingS,HStabS,VStabS,CDWing,Cw,CDHStab,Chstab,CDVStab,Cvstab,FSA,SAfT,Lf,Wf,SensorCalc,ARs,Ws)
%A function brought to you by Bryce Moran and Matthew Maddalon.
%don't panic (!!!)

%List of Inputs
    %G is vector input from landingGear.m
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
%-------------------------------Constants---------------------------------%
    % NOT USED? mu=0.3776e-6; %for 19 Centigrade %viscous friction coefficient. calculated for a certain temperature
    muk=1.612e-4; %kinmatic viscosity (ish) related to mu                                  
%-------------------------------------------------------------------------%
%----------------------Induced Drag from the Wing-------------------------%
K=1/(pi*0.8*AR); %e=0.8. K is an aero value. e is efficiency rating for the wing (rectangular)
Induced= @(v)(2*K*WeightT^2)/(rho*WingS*v^2); %It's the induced Drag. doesn't account for lift from tail.
%-------------------------------------------------------------------------%
%---------------------------Parasitic Drag--------------------------------%
CDf=0.75; %good round estimate
wingParasitic = @(v) (0.5*rho*(v^2)*WingS*CDWing);
hStabParasitic = @(v) (0.5*rho*(v^2)*HStabS*CDHStab);
vStabParasitic = @(v) (0.5*rho*(v^2)*VStabS*CDVStab);
gearWheelParasitic = G(2);   
fuselageParasitic = @(v) (0.5*rho*(v^2)*FSA*CDf);
Parasitic = wingParasitic(v) + hStabParasitic(v) + vStabParasitic(v) + gearWheelParasitic(v) + fuselageParasitic(v);
%-------------------------------------------------------------------------%
%---------------------------Skin Friction---------------------------------%
%Skin friction Drag( sum: wing, HStab, VStab+...
Temp=273+21; %kelvin - average april weather at competition site.
R=287; %gas consant - no need to change year to year
gamma=1.4; %air constant - no change year-to-year
a=3.28*sqrt(gamma*R*Temp);%with metric to imperial conversion. a is speed of sound. 3.28 is conversion factor to feet
Cfl= @(v,L) 1.328/sqrt(v*L/muk);%coeff friction laminar
Cft=@(v,L) (0.455)/((log10(v*L/muk)^2.58)*(1+0.144*(v/a)^2)^0.65); %coeff friction turbulent
Cf=@(v,L,k) k*Cfl(v,L)+(1-k)*Cft(v,L); %trueish - some ratio of laminar and turbulent
skinFrictionWing = @(v) 2*(0.5*rho*(v^2)*Cf(v,Cw,1)*WingS);
skinFrictionHStab = @(v) 2*(0.5*rho*(v^2)*Cf(v,Chstab,1)*HStabS);
skinFrictionVStab = @(v) 2*(0.5*rho*(v^2)*Cf(v,Cvstab,1)*VStabS);
skinFrictionFuselage = @(v) (0.5*rho*(v^2)*Cf(v,Lf,0.9)*SAfT);
skinFrictionGear = @(v) 2*(0.5*rho*(v^2)*Cf(v,(2/12) ,0)*G(3));
skinFrictionWheels = @(v) 3*(0.5*rho*(v^2)*Cf(v,2*radius_wheel/12,0)*WheelSA);
Skin= @(v) skinFrictionWing(v) + skinFrictionHStab(v) + skinFrictionVStab(v) + skinFrictionFuselage(v) + skinFrictionGear(v) + skinFrictionWheels(v);
%coeff skin friction. everything has a unique one
%k for not fuselage is about 0.1
%k (a constant) for fuselage should be 0 (all turbulent)
%k, the constant you have no idea about, is the transition point along the
%length from laminar to tubulent flow
%-Bryce
%prop drag - accounts for accelerated flow across the rest of the aircraft
%all this for like 0.001 change :/  
%-Bryce

%DEBUG: Commment out
%-------------------------------------------------------------------------%
%-------------------------------Prop Drag---------------------------------%
%     Ap=pi*(dp/2)^2;
%     Sl=SAfT + (dp-Wf)*12*Cw;%total wet fuselage Surface Area in prop wash
%     q=@(v) 0.5*rho*v^2;% is q
%     Dcl=@(v) 0.5*(Cw/dp)*CLw* Thrust/(q(v) * Ap);%change in CL
%     dci=@(v) 2*Dcl(v)*CLw/(pi*AR);%change in Cdi
%     Vs=@(v) (v+Thrust/(rho*Ap*v));%velocity of prop slipstream
%         SAT=SAfT + 2*(dp-Wf)*12*Cw;
%         kwing=(2*(dp-Wf)*12*Cw)/SAT;
%         kfuse=SAfT/SAT;
%     CFp=@(v) kwing*Cf(Vs(v),Cw,0.1) + kfuse*(Cf(Vs(v),Lf,0));  
%     DC0=@(v) Cf(v,Lf,0)*(Sl/WingS)*( ((Vs(v))/v)^2 -1);%change in Cd0
%     Prop=@(v) (DC0(v)+dci(v)) *0.5*rho*(Vs(v))^2;%chnage in drag (additional drag)


%-------------------------------------------------------------------------%
%-----------------------------Total Drag----------------------------------%
Drag= @(v)   Induced(v)+Parasitic(v) +Skin(v); 
%-------------------------------------------------------------------------%
%---------------------------Dynamic Thrust--------------------------------%
% PS=pitch/12 * 5614 /60; %propeller speed
% PD=pitch/dp; %this is dynamic thrust function. Article this is based off of is in resources google drive somewhere. %DEBUG - does this match the physics of the motor?
% FS=@(v) v; %DEBUG: research this...
% 
% Ts=Thrust;
% s=0;
% if PD<0.6
%     T= @(v) Ts*(1 - FS(v)/(PS*(PD+0.2)/PD));
% elseif (s/PS)*PD<(PD-0.6)
%     T=@(v) Ts;
% else
%     T=@(v) Ts*(1- ((FS(v)*PD/PS)-(PD-0.6))/0.8);
% end

%-------------------------------------------------------------------------%
%--------------------------Velocity Solver--------------------------------%
%false position root finder #1 to use initially, try the other when this breaks. V has to be between xl and xu
%FS=@(v) v; needed?  + FS(v) - FS(v);
T = @(v) 
func=@(v) (Drag(v)-Thrust);
xl=-30; xu=1000; ea=1e-4;
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
%-------------------------------------------------------------------------%
%Emergency, if stuff breaks, its because of induced drag. Solver reverts to a
%different form of induced drag with the assumption that L=W. Without this
%it returns complex numbers :(. Happens when velocity becomes so small it's
%computationally zero. %DEBUG
if V<0
    V=-1; %DEBUG - make an error flag - functions look for exactly this value for now
elseif V>900 %no way we ever get this fast.
    V=-1;
end
if V==-1; %new induced drag if the other broke
    Induced=@(v) WeightT*(CLw)*K*WingS;
end
Drag= @(v)   Induced(v)+Parasitic(v) +Skin(v); 
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
    V=-1;
elseif V>900
    V=-1; %DEBUG -- v == 1 is still the it's broke flag
end
%DEBUG - add broken down drag values here - return an array for sanity
%check purposes. return drag for each part rather than overall.
%DEBUG - look into describing drag in terms of parts instead of classes
%(wing instead of skin alone)
D=Drag(V); %function handle that sums all drag
% iter
end