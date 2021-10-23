function G = landingGear(propDiam, fuselageHeight, fuselageWidth,rho,span)

%This function computes different parameters for the landing gear
%Inputs: PropDiam = diameter of propellor (inches)
%        fuselageHeight = height of the fuselage (inches)
%Outputs: A vector G with the following values:
%   G = [gearWeight gearParaDrag gearSA wheelSA]
%Height of landing gear is the difference between prop radius and fuselage height
%
%                                top
%         ^              -------------------
%         |            / (theta+90)          \
%         |           /                       \
%         |        c /                         \
%         h         /                           \ 
%         |        /                             \
%         |       |                               |
%         |       |                               |
%         |       |                               |
%
%
%         -        <---------------b--------------->
%
%-------------------------------------------------------------------------%
%----------------------Calculation of Gear Weight-------------------------%
height = .5*propDiam - fuselageHeight + 1; %1" clearence 
top = fuselageWidth;
theta = 38*pi/180;   %Angle of inside of bend (minus 90)
rho_Gear = 0.001122368;   %Density for landing gear material (Al) (slug/in^3)
t = 1/8;        %thickness of flatbar (in)
b = 0.25*span;
a = .58(b-top);
c = a/sin(theta);
d = c*cos(theta);
L = top+2*c+2*(height-d);
gearFrontArea= L*t;
gearwidth = 1;      %width of flatbar (in)
volume_Al = gearFrontArea*gearwidth;
gearWeight = rho_Gear*volume_Al;
%-------------------------------------------------------------------------%
%-------------------Calculation of Gear Surface Area----------------------%
width_wheel = 0.5;    %width of wheels (in)
radius_wheel = 1.5; %wheel radius
gearSA=L*gearwidth;
wheelSA=(2*pi*(radius_wheel)^2+ pi*2*radius_wheel*width_wheel)/144;
%-------------------------------------------------------------------------%
%------------------------Parasitic Drag Function--------------------------%
Cd_Al=0.05;         % Let Cd = 0.05 for the flatbar
Cd_wheel = 0.20;     % Let Cd = 0.245 for the wheels
rho_air = rho/(12^3);  %density of air at sea level in slug/in^3   CONVERSION
area_wheel = width_wheel*radius_wheel*3*2;    %profile area for 3 wheels
drag_Al =@(v) Cd_Al*rho_air*(v*12)^2*gearFrontArea*.5; %aluminum bar
drag_wheel = @(v) Cd_wheel*rho_air*(v*12)^2*2*area_wheel*.5;
gearParaDrag=@(v) drag_Al(v)+drag_wheel(v); %p for parasitic
%-------------------------------------------------------------------------%
%--------------------------------Output-----------------------------------%
G = [gearWeight gearParaDrag gearSA wheelSA];
%-------------------------------------------------------------------------%

