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

PD=pitch/dp; %this is dynamic thrust function. Article this is based off of is in resources google drive somewhere. %DEBUG - does this match the physics of the motor?
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
%computationally zero. %DEBUG
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
%check purposes. return drag for each part rather than overall.
%DEBUG - look into describing drag in terms of parts instead of classes
%(wing instead of skin alone)
D=Drag(V); %function handle that sums all drag
% iter
end
