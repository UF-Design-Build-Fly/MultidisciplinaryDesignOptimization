function [plane] = landingGear(plane, rho)

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
%                  <--------------b-------------->
%
%-------------------------------------------------------------------------%
%----------------------Calculation of Gear Weight-------------------------%
propDiam = plane.power.propDiameter;
fuselageHeight = plane.fuselage.height;
height = .5*propDiam - fuselageHeight + 2; %2" clearence 
top = plane.fuselage.width;
span = plane.wing.span;
theta = 38*pi/180;   %Angle of inside of bend (minus 90)
rho_Gear = 0.001122368;   %Density for landing gear material (Al) (slug/in^3)
t = 1/8;        %thickness of flatbar (in)
b = 0.25*span;
a = .5*(b-top);
c = a/sin(theta);
d = c*cos(theta);
L = top+2*c+2*(height-d);
gearFrontArea= L*t;
gearwidth = 2;      %width of flatbar (in)
volume_Al = gearFrontArea*gearwidth;
plane.fuselage.gearWeight = rho_Gear*32.1740*volume_Al;
%-------------------------------------------------------------------------%
%-------------------Calculation of Gear Surface Area----------------------%
plane.fuselage.gearSA=L*gearwidth;
%-------------------------------------------------------------------------%
%------------------------Parasitic Drag Function--------------------------%
Cd_Al=0.05;         % Let Cd = 0.05 for the flatbar
Cd_wheel = 0.20;     % Let Cd = 0.245 for the wheels
rho_air = rho/(12^3);  %density of air at sea level in slug/in^3   CONVERSION

area_wheel = plane.fuselage.wheelWidth*plane.fuselage.wheelRadius*3*2;    %profile area for 3 wheels
%Bryce edit, taking out the 12's in the below equation %DEBUG
drag_Al =@(v) Cd_Al*rho_air*((v)^2)*gearFrontArea*.5; %aluminum bar
drag_wheel = @(v) Cd_wheel*rho_air*((v)^2)*area_wheel*.5;
plane.fuselage.gearParaDrag=@(v) drag_Al(v)+drag_wheel(v); %p for parasitic
%-------------------------------------------------------------------------%