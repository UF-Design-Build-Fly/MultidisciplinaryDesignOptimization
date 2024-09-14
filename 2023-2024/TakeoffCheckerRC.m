function plane = TakeoffCheckerRC(plane,mission,rho)

    %DEBUG - this function might not actually work. gets an idea but not sure of
    %accuracy. use at your own risk.

    if (mission == 2)
        weight = plane.performance.totalWeight2; %Mission 2 weight
    elseif (mission == 3)
        weight = plane.performance.totalWeight3; %Mission 3 weight
    end

    %calculate velocity at end of takeoff distance
	%see if there is an aoa that generates enough lift @ said velocity to takeoff
	%include vertical thrust of motor @ aoa
	%max aoa limited to 20 degrees increments of 1 degree if lift-weight < ?2?lb increments of 0.1
    
    WingS = plane.wing.planformArea; %reference area wing (ft^2)
    %fh = fuselage height (in)  %removed, originally used for wing height,
                              %wing height is now just prop radius
    %RPM = Motor RPM, max      %not sure why rpm was necessary
    %Pitch = propeller pitch   %removed because dynamic thrust equation is
                              %not used
    %constants
    g = 32.2; %acceleration due to gravity (ft/s/s)
    mu = 0.008; %coeff rolling friction
    e = 0.8; %efficiency of wing

    K = 1/(pi*plane.wing.aspectRatio*e);    %Wing K constant
    
    h = plane.powerSystem.propDiameter/2 + 1/6; %Height of wing above ground (ft) 1/6 is clearence
    phi = 1 -(2*e/pi)*log(1+ ( (pi*plane.wing.chord)/(8*h))^2); %scaling constant
    Kg = phi*K;%K constant, ground effect
    Clg = mu/(2*Kg); %coeff lift, ground effect
    Kw = 1.4817*(5.81e-5);%Intermediate constant
    dCdlg = (weight/WingS)*(Kw)*(weight/g); %change in CD due to ground effect
    CDlg = plane.wing.cd-dCdlg;%CD0 ground effect
    CDg = CDlg+Kg*Clg^2;
    
    %Ap was not used in anything so idk why it was calculated - Christian
    %Ap = (pi*0.25*dp^2)/144; %Prop area(ft^2)
    Vr = 1.2*sqrt(2*weight/(rho*(plane.wing.clm+0.9)*WingS)); %Rotation speed; with 1.2 factor of safety
    
    % We decided dynamic thrust was not necessary- Christian
    % %dynamic thrust equation
    % PS = pitch/12 * 5614 /60;
    % PD = pitch/dp;
    % FS = @(v) v;
    % if PD<0.6
    %     T = @(v) Ts*(1 - FS(v)/(PS*(PD+0.2)/PD));
    % elseif (FS(v)/PS)*PD<(PD-0.6)
    %     T = @(v) Ts;
    % else
    %     T = @(v) Ts*(1- ((FS(v)*PD/PS)-(PD-0.6))/0.8);
    % end
    
    %Groundroll takeoff distance equation
    A = g*(plane.powerSystem.thrust/weight - mu); %A constant
    %B constant. absolutely ignore that little a or nothing works
    B = (g/weight)*(0.5*rho*WingS*(CDg-mu*Clg) );% +a);
    
    distance = (1/(2*B))*log(A/(A-B*Vr^2)); %ground roll (ft)
    
    if (mission == 2)
        plane.performance.takeoffDist2 = distance; %Mission 2 takeoff dist
    elseif (mission == 3)
        plane.performance.takeoffDist3 = distance; %Mission 3 takeoff dist
    end

end