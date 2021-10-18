function D= TakeoffChecker(Thrust,W,rho,WingS,AR,Cwing,CLm,CD0,dp,fh)%,RPM,pitch)
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
    %RPM =Motor RPM, max   %not sure why rpm was necessary
    %Pitch= propeller pitch   %removed because dynamic thrust equation is
                              %not used
%constants
    g=32.2; %acceleration due to gravity (ft/s/s)
    mu=0.008; %coeff rolling friction
    e=0.8; %efficiency of wing

K=1/(pi*AR*e);    %Wing K constant
%DEBG - this height equation is wrong, it should be fueselage height +
%landing gear height but idk what the landing gear height is
h=(dp/2+0.5+fh)/12; %Height of wing above ground (ft)
phi= 1 -(2*e/pi)*log(1+ ( (pi*Cwing)/(8*h))^2); %scaling constant
Kg=phi*K;%K constant, ground effect
Clg=mu/(2*Kg); %coeff lift, ground effect
Kw=1.4817*(5.81e-5);%Intermediate constant
dCd0lg=(W/WingS)*(Kw)*(W/g); %change in CD0 due to ground effect
CD0lg=CD0-dCd0lg;%CD0 ground effect
CDg=CD0lg+Kg*Clg^2;

%Ap was not used in anything so idk why it was calculated - Christian
%Ap=(pi*0.25*dp^2)/144; %Prop area(ft^2)
Vr=1.2*sqrt(2*W/(rho*CLm*WingS)); %Rotation speed; with 1.2 factor of safety

% We decided dynamic thrust was not necessary- Christian
% %dynamic thrust equation
% PS=pitch/12 * 5614 /60;
% PD=pitch/dp;
% FS=@(v) v;
Ts=Thrust;
% if PD<0.6
%     T= @(v) Ts*(1 - FS(v)/(PS*(PD+0.2)/PD));
% elseif (FS(v)/PS)*PD<(PD-0.6)
%     T=@(v) Ts;
% else
%     T=@(v) Ts*(1- ((FS(v)*PD/PS)-(PD-0.6))/0.8);
% end

%Groundroll takeoff distance equation
A=g*(Ts(Vr)/W - mu); %A constant
%B constant. absolutely ignore that little a or nothing works
B=(g/W)*(0.5*rho*WingS*(CDg-mu*Clg) );% +a);
D=(1/(2*B))*log(A/(A-B*Vr^2)); %ground roll (ft)
end