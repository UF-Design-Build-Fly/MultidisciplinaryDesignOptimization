function G = landingGear(propDiam, fuselageHeight, fuselageWidth)

%This function computes different parameters for the landing gear
%Inputs: PropDiam = diameter of propellor (inches)
%        fuselageHeight = height of the fuselage (inches)
%Outputs: A vector G with the following values:
%  [gearWeight gearParaDrag gearSkinFriction

%Constants and geometry:
rho_Gear = 1;   %Density for landing gear material (slug/in^3)
top = fuselageWidth;
t = 1/8;            %thickness of flatbar (in)
theta = 38*pi/180;  %Angle of inside of bend (minus 90)
gearwidth = 1;      %width of flatbar (in)
width_wheel = 0.5;    %width of wheels (in)
radius_wheel = 1.5; %wheel radius


%Height of landing gear is the difference between prop radius and fuselage height
height = .5*propDiam - fuselageHeight + 1; %1" clearence 

area_Al= (top + 2*(height/cos(theta)+height))*t;
volume_Al = area_Al*gearwidth;
gearWeight = rho_Gear*volume_Al;
area_wheel = width_wheel*radius_wheel*3*2;    %profile area for 3 wheels
drag_Al =@(v) Cd_Al*rho_air*(v*12)^2*area_Al*.5; %aluminum bar
drag_wheel = @(v) Cd_wheel*rho_air*(v*12)^2*2*area_wheel*.5;
gearParaDrag=@(v) drag_Al(v)+drag_wheel(v); %p for parasitic
G = [gearWeight gearParaDrag

