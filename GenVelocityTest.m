%%
%%%velocity solver via false position root finding
function [plane] = GenVelocityTest(plane,missionNumber,rho)
%G, Thrust,CLw,WeightT,rho,AR,WingS,HStabS,VStabS,CDWing,Cw,CDHStab,Chstab,CDVStab,Cvstab,FSA,SAfT,Lf
%RPM,pitch,dp,CLmax)
%A function brought to you by Bryce Moran and Matthew Maddalon.
%Don't panic (!!!)

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
%CLflaps = max coefficient of lift for wing with flaps down
%Outputs:
%V = cruise velocity
%landingV = landing speed

Lf = plane.fuselage.length; %ft
SAfT = plane.fuselage.totalSA; %ft^2
radius_wheel = plane.fuselage.wheelRadius; %in
CLw = plane.wing.clw; %dimensionless

%-------------------------------Constants---------------------------------%
% NOT USED? mu=0.3776e-6; %for 19 Centigrade %viscous friction coefficient. calculated for a certain temperature
muk=1.612e-4; %kinmatic viscosity (ish) related to mu
%-------------------------------------------------------------------------%
if missionNumber == 1
    WeightT = plane.performance.totalEmptyWeight;
elseif missionNumber == 2
    WeightT = plane.performance.totalWeight2;
else
    WeightT = plane.performance.totalWeight3;
end
AR = plane.wing.aspectRatio;
WingS = plane.wing.planformArea; %ft^2
%----------------------Induced Drag from the Wing-------------------------%
K=1/(pi*0.8*AR); %e=0.8. K is an aero value. e is efficiency rating for the wing (rectangular)
Induced= @(v)(2*K*WeightT^2)/(rho*WingS*v^2); %It's the induced Drag. doesn't account for lift from tail.
%-------------------------------------------------------------------------%
%---------------------------Parasitic Drag--------------------------------%
CDHStab = plane.empennage.HScd;
CDVStab = plane.empennage.VScd;
VStabS = plane.empennage.VSarea; %relative to wing planform area - see empennage.m
HStabS = plane.empennage.HSarea;
FSA = plane.fuselage.frontalSurfaceArea;
CDWing = plane.wing.cd;
CDf=0.75; %good round estimate
wingParasitic = @(v) (0.5*rho*(v^2)*WingS*CDWing);
hStabParasitic = @(v) (0.5*rho*(v^2)*HStabS*CDHStab);
vStabParasitic = @(v) (0.5*rho*(v^2)*VStabS*CDVStab);
fuselageParasitic = @(v) (0.5*rho*(v^2)*FSA*CDf);
Parasitic = @(v) wingParasitic(v) + hStabParasitic(v) + vStabParasitic(v) + fuselageParasitic(v) + plane.fuselage.gearParaDrag(v);
%-------------------------------------------------------------------------%
%---------------------------Skin Friction---------------------------------%
%Skin friction Drag( sum: wing, HStab, VStab+...
Temp=294; %kelvin - average april weather at competition site.
R=287; %gas consant - no need to change year to year
gamma=1.4; %air constant - no change year-to-year
a=3.28*sqrt(gamma*R*Temp);%with metric to imperial conversion. a is speed of sound. 3.28 is conversion factor to feet
Chstab = plane.empennage.HSchord;
Cvstab = plane.empennage.VSchord;
Cfl= @(v,L) 1.328/sqrt(v*L/muk);%coeff friction laminar
Cft=@(v,L) (0.455)/((log10(v*L/muk)^2.58)*(1+0.144*(v/a)^2)^0.65); %coeff friction turbulent
Cf=@(v,L,k) k*Cfl(v,L)+(1-k)*Cft(v,L); %trueish - some ratio of laminar and turbulent
skinFrictionWing = @(v) 2*(0.5*rho*(v^2)*Cf(v,CDWing,1)*WingS);
skinFrictionHStab = @(v) 2*(0.5*rho*(v^2)*Cf(v,Chstab,1)*HStabS);
skinFrictionVStab = @(v) 2*(0.5*rho*(v^2)*Cf(v,Cvstab,1)*VStabS);
skinFrictionFuselage = @(v) (0.5*rho*(v^2)*Cf(v,Lf,0.9)*SAfT);
skinFrictionGear = @(v) 2*(0.5*rho*(v^2)*Cf(v,(2/12) ,0)*plane.fuselage.gearSA);
skinFrictionWheels = @(v) 3*(0.5*rho*(v^2)*Cf(v,2*radius_wheel/12,0)*plane.fuselage.wheelSA);
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
FS=@(v) v; %DEBUG: research this...
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
Thrust = plane.power.thrust;
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
% plane.performance.drag1 = Drag(v);
% plane.performance.inducedDrag = Induced(v)
% plane.performance.parasiticDrag = Parasitic(v)
% plane.performance.skinDrag = Skin(v);
% plane.performance.wingPara = wingParasitic(v);
% plane.performance.hStabPara = hStabParasitic(v);
% plane.performance.vStabPara = vStabParasitic(v);
% plane.performance.fusePara = fuselageParasitic(v);
% plane.performance.gearPara = plane.fuselage.gearParaDrag(v);

%Parasitic = @(v)  + hStabParasitic(v) + vStabParasitic(v) + fuselageParasitic(v) + plane.fuselage.gearParaDrag(v);

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
if V==-1 %new induced drag if the other broke
    Induced=@(v) WeightT*(CLw)*K*WingS;
end
Drag= @(v)   Induced(v)+Parasitic(v) +Skin(v);
func=@(v) (Drag(v)-Thrust);
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
Vstall = sqrt(2*WeightT/(rho*WingS*plane.wing.clFlap));
landingV = 1.3*Vstall;
if missionNumber == 1
    plane.performance.velocity1 = V;
    plane.performance.landingSpeed1 = landingV;
    plane.performance.drag1 = D;
elseif missionNumber == 2
    plane.performance.velocity2 = V;
    plane.performance.landingSpeed2 = landingV;
    plane.performance.drag2 = D;
else
    plane.performance.velocity3 = V;
    plane.performance.landingSpeed3 = landingV;
    plane.performance.drag3 = D;
end
