function [plane] = TakeoffChecker(plane,M,rho)%,fh,RPM,pitch)
%DEBG - this function might not actually work. gets an idea but not sure of
%accuracy. use at your own risk.
%Inputs
    Ts=plane.power.thrust; %thrust
    %M = mission number
    if M==2
        W=plane.performance.totalWeight2; %Mission 2 weight
    elseif M==3
        W=plane.performance.totalWeight3; %Mission 3 weight
    end
    %rho=density air (slug/ft^3)
    WingS=plane.wing.planformArea; %reference area wing (ft^2)
    AR= plane.wing.aspectRatio;%aspect ratio, wing
    Cwing=plane.wing.chord; %chord wing (ft)
    CLm=plane.wing.clm; %coeff lift, max
    CD=plane.wing.cd; %wing cd @ chosen alpha
    dp=plane.power.propDiameter; %prop diameter (in)
    %fh=fuselage height (in)  %removed, originally used for wing height,
                              %wing height is now just prop radius
    %RPM =Motor RPM, max      %not sure why rpm was necessary
    %Pitch= propeller pitch   %removed because dynamic thrust equation is
                              %not used
%constants
    g=32.2; %acceleration due to gravity (ft/s/s)
    mu=0.008; %coeff rolling friction
    e=0.8; %efficiency of wing

K=1/(pi*AR*e);    %Wing K constant

h=((dp/2)+2)/12; %Height of wing above ground (ft)
phi= 1 -(2*e/pi)*log(1+ ( (pi*Cwing)/(8*h))^2); %scaling constant
Kg=phi*K;%K constant, ground effect
Clg=mu/(2*Kg); %coeff lift, ground effect
Kw=1.4817*(5.81e-5);%Intermediate constant
dCdlg=(W/WingS)*(Kw)*(W/g); %change in CD due to ground effect
CDlg=CD-dCdlg;%CD0 ground effect
CDg=CDlg+Kg*Clg^2;

%Ap was not used in anything so idk why it was calculated - Christian
%Ap=(pi*0.25*dp^2)/144; %Prop area(ft^2)
Vr=1.2*sqrt(2*W/(rho*CLm*WingS)); %Rotation speed; with 1.2 factor of safety

% We decided dynamic thrust was not necessary- Christian
% %dynamic thrust equation
% PS=pitch/12 * 5614 /60;
% PD=pitch/dp;
% FS=@(v) v;
% if PD<0.6
%     T= @(v) Ts*(1 - FS(v)/(PS*(PD+0.2)/PD));
% elseif (FS(v)/PS)*PD<(PD-0.6)
%     T=@(v) Ts;
% else
%     T=@(v) Ts*(1- ((FS(v)*PD/PS)-(PD-0.6))/0.8);
% end

%Groundroll takeoff distance equation
A=g*(Ts/W - mu); %A constant
%B constant. absolutely ignore that little a or nothing works
B=(g/W)*(0.5*rho*WingS*(CDg-mu*Clg) );% +a);

D=(1/(2*B))*log(A/(A-B*Vr^2)); %ground roll (ft)
if M==2
        plane.performance.takeoffDist2=D; %Mission 2 takeoff dist
    elseif M==3
        plane.performance.takeoffDist3=D; %Mission 3 takeoff dist
    end
end