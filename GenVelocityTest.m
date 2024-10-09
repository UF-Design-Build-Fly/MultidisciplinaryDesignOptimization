function plane = GenVelocityTest(plane, missionNumber, rho, temp)%, dThrustNeuralNet, dThrustStats)
%DEBUG: this sometimes finds complex numbers for drag and fails to find velocity. Is this a bug or just when it gets an impractical (underpowered) aircraft?

    %-------------------------------Constants---------------------------------%
    muk = 15.97e-6; %kinmatic viscosity (ish) related to mu (m2/s) @303K
    
    %--------------------------------Weight-----------------------------------%
    if missionNumber == 2
        weightTotal = plane.performance.totalWeight2;
    elseif missionNumber == 3
        weightTotal = plane.performance.totalWeight3;
    end
    
    %-----------------------------Induced Drag--------------------------------%
    % wingEfficiency = 0.8; %e is efficiency rating for the wing (rectangular)
    wingEfficiency = 1.78*(1-0.045*plane.wing.aspectRatio.^0.68)-0.64; % assuming rectangular wing
    %InducedDrag = @(v) (((2*weightTotal)/(rho*plane.wing.planformArea*(v.^2))).^2)/(pi*wingEfficiency*plane.wing.aspectRatio); %Doesn't account for lift from tail
    InducedDrag = @(v) (2*weightTotal^2)/(rho*(v.^2)*pi*plane.wing.span^2*wingEfficiency); %Doesn't account for lift from tail
    
    %---------------------------Parasitic Drag--------------------------------%
    CDf = 0.75; %good fuselage estimate
    wingParasitic = @(v) (0.5*rho*(v^2)*plane.wing.planformArea*plane.wing.cd);
    hStabParasitic = @(v) (0.5*rho*(v^2)*plane.empennage.HSarea*plane.empennage.HScd);
    vStabParasitic = @(v) (0.5*rho*(v^2)*plane.empennage.VSarea*plane.empennage.VScd);
    fuselageParasitic = @(v) (0.5*rho*(v^2)*plane.fuselage.frontalSurfaceArea*CDf);
    
    Cd_gear = 0.05;        % Let Cd = 0.05 for the flatbar (front part of tricycle)
    Cd_wheel = 0.245;     % Let Cd = 0.245 for the wheels
    dragBar = @(v) (0.5*rho*(v)^2*plane.fuselage.gearFrontalSA*Cd_gear); %bar
    dragWheels = @(v) (0.5*rho*(v)^2*plane.fuselage.wheelFrontalSA*Cd_wheel); %before Kevin's edits, v=40?, not sure why
    
    gearParaDrag = @(v) dragBar(v)+dragWheels(v); %p for parasitic. %DEBUG -- GenVelocityTest will call this function - but through a really convoluted series of function calls
                                                                                             %Maybe try to just set constants here and have GenVelocityTest compute drag directly?
    
    Parasitic = @(v) wingParasitic(v) + hStabParasitic(v) + vStabParasitic(v) + fuselageParasitic(v) + gearParaDrag(v);
    
    
    %-------------------------------------------------------------------------%
    
    
    %---------------------------Skin Friction drag----------------------------%
    %coeff skin friction. everything has a unique one
    %k for not fuselage is about 0.1
    %k (a constant) for fuselage should be 0 (all turbulent)
    %k, the constant you have no idea about, is the transition point along the
    %length from laminar to tubulent flow
    %-Bryce
    
	% Speed of sound in air (m/s)
    R = 287; %gas consant - no need to change year to year
    gamma = 1.4; %air constant - no change year-to-year
    a = sqrt(gamma*R*temp);%with metric to imperial conversion. a is speed of sound. 3.28 is conversion factor to feet
    %a = 331*sqrt(temp/273);
    
    Cfl = @(v,L) 1.328/sqrt(v*L/muk);%coeff friction laminar
    Cft = @(v,L) (0.455)/((log10(v*L/muk)^2.58)*(1+0.144*(v/a)^2)^0.65); %coeff friction turbulent %I'm not familiar with this equation. There are a dozen others to pick from. Source: https://www.cfd-online.com/Wiki/Skin_friction_coefficient 
    Cf = @(v,L,k) k*Cfl(v,L)+(1-k)*Cft(v,L); %trueish - some ratio of laminar and turbulent
    skinFrictionWing = @(v) 2*(0.5*rho*(v^2)*Cf(v,plane.wing.chord,1)*plane.wing.planformArea);
    skinFrictionHStab = @(v) 2*(0.5*rho*(v^2)*Cf(v,plane.empennage.HSchord,1)*plane.empennage.HSarea);
    skinFrictionVStab = @(v) 2*(0.5*rho*(v^2)*Cf(v,plane.empennage.VSchord,1)*plane.empennage.VSarea);
    skinFrictionFuselage = @(v) (0.5*rho*(v^2)*Cf(v,plane.fuselage.length,0.9)*plane.fuselage.totalSA);
    skinFrictionGear = @(v) (0.5*rho*(v^2)*Cf(v,(1.5*0.0254) ,0)*plane.fuselage.gearSA);
    skinFrictionWheels = @(v) (0.5*rho*(v^2)*Cf(v,2*plane.fuselage.wheelRadius,0)*plane.fuselage.wheelSA);

    Skin = @(v) skinFrictionWing(v) + skinFrictionHStab(v) + skinFrictionVStab(v) + skinFrictionFuselage(v) + skinFrictionGear(v) + skinFrictionWheels(v);
    %-------------------------------------------------------------------------%
    
    
    %-------------------------------Prop Drag---------------------------------%
    %DEBUG: Commment out. Used for aircraft configs where the prop is in front
    %of the fuselage. This year used dual tractor so no additional drag from
    %prop wash over fuselage.
    %prop drag - accounts for accelerated flow across the rest of the aircraft
    %all this for like 0.001 change :/
    %-Bryce
    %     Ap = pi*(dp/2)^2;
    %     Sl = SAfT + (dp-Wf)*12*Cw;%total wet fuselage Surface Area in prop wash
    %     q = @(v) 0.5*rho*v^2;% is q
    %     Dcl = @(v) 0.5*(Cw/dp)*CLw* Thrust/(q(v) * Ap);%change in CL
    %     dci = @(v) 2*Dcl(v)*CLw/(pi*AR);%change in Cdi
    %     Vs = @(v) (v+Thrust/(rho*Ap*v));%velocity of prop slipstream
    %         SAT = SAfT + 2*(dp-Wf)*12*Cw;
    %         kwing = (2*(dp-Wf)*12*Cw)/SAT;
    %         kfuse = SAfT/SAT;
    %     CFp = @(v) kwing*Cf(Vs(v),Cw,0.1) + kfuse*(Cf(Vs(v),Lf,0));
    %     DC0 = @(v) Cf(v,Lf,0)*(Sl/plane.wing.planformArea)*( ((Vs(v))/v)^2 -1);%change in Cd0
    %     Prop = @(v) (DC0(v)+dci(v)) *0.5*rho*(Vs(v))^2;%chnage in drag (additional drag)
    
    %-------------------------------------------------------------------------%
    
    
    
    %-----------------------------Total Drag----------------------------------%
    TotalDrag = @(v) InducedDrag(v) + Parasitic(v) + Skin(v);
    %-------------------------------------------------------------------------%
    
    
    %---------------------------Dynamic Thrust--------------------------------%
    %This function is only an estimate based on limited available testing data.
    %It should largely be considered an underestimator of performance. 
    %This largely determines cruise performance
    
    %Key assumptions are that the propeller spins at a constant rpm and that
    %air exit velocity is constant through the whole propeller radius. The RPM
    %changes suprisingly little with forward airspeed. The "propeller speed" is
    %by far the most important factor in determining max dynamic thrust.
    %See https://www.electricrcaircraftguy.com/2014/04/propeller-static-dynamic-thrust-equation-background.html
    
    % PS = (plane.powerSystem.propPitch/12) * (plane.powerSystem.rpm/60); %propeller speed
    % PD = plane.powerSystem.propPitch/plane.powerSystem.propDiamter;
    % 
    % Ts = plane.powerSystem.thrust;
    % s = 0;
    % if PD<0.6
    %     T = @(v) Ts*(1 - FS(v)/(PS*(PD+0.2)/PD));
    % elseif (s/PS)*PD<(PD-0.6)
    %     T = @(v) Ts;
    % else
    %     T = @(v) Ts*(1- ((FS(v)*PD/PS)-(PD-0.6))/0.8);
    % end
    
     %k1 = 4.39*10^-8; %semi-empirical model coefficients
     %k2 = 4.23*10^-4;
     %rpm = plane.powerSystem.rpm;
     %d = plane.powerSystem.propDiameter;
     %pitch = plane.powerSystem.propPitch;
     %N_lb = 0.224; %newtons to pounds conversion factor
    

    %T = @(v) min(plane.powerSystem.thrust, ((N_lb*k1)*rpm*((d^3.5)/sqrt(pitch)))*((k2)*rpm*pitch-(v)));
    
    
    %Empirical Data update: This function, for the cobra 4130/20 with a 20x13
    %prop, effectively finds the prop speed (61 mph) as the zero thrust
    %airspeed. In fact the motor is still producing the roughly 4 pounds of
    %thrust needed to maintain the measured maximum forward airspeed of ~55 mph
    %Motor zero thrust should be at 120 feet per second.
    %-------------------------------------------------------------------------%
    %Nathaniel's Very crude but somewhat accurate static thrust

    % Prop speed from static thrust is lower than at max speed
    % 3.7V assumed under load
    % propPitch in inches(to meters 0.0254)
    flightPropSpeed = 3.7 * plane.powerSystem.cells * plane.powerSystem.kv/60 * plane.powerSystem.propPitch * 0.0254;
    CrudeDynamicThrust = @(v) (9.80665*plane.powerSystem.thrust*(1 - v/flightPropSpeed));
    
    %--------------------------Velocity Solver--------------------------------%
    
    %func = @(v) (Drag(v)-T(v));
    %propDiameter = plane.powerSystem.propDiameter;
    %propPitch = plane.powerSystem.propPitch;
    %motorRPM = plane.powerSystem.rpm;
    %CalcNetForce = @(v) (TotalDrag(v) - CalcDynamicThrust(propDiameter, propPitch, motorRPM, v, dThrustNeuralNet, dThrustStats));
    CalcNetForce = @(v) (TotalDrag(v) - CrudeDynamicThrust(v));
    
    %i = 1; %debug code for visualizing thrust/drag curves
    %for v = 5:0.01:100
    %    d(i) = Drag(v);
    %    f(i) = func(v);
    %    w(i) = CalcDynamicThrust(d,p,rpm,v,net,stats);
    %    i = i+1;
    %    %I_d(i) = Induced(v);
    %    %P_d(i) = Parasitic(v);
    %    %S_d(i) = Skin(v);
    %end
    %v = 5:0.01:100;
    %hold on;
    %plot(v,d);
    %plot(v, f);
    %plot(v, w);
    %legend("Drag", "Net Force", "Thrust");
    
    velHigh = plane.powerSystem.propSpeed;
    velLow = 15;
    velocity = 0.85*velHigh;
    netForce = CalcNetForce(velocity);
    iterations = 0;
    while (netForce >= 0.05 || netForce <= -0.05)

        velocity = 0.5*(velLow+velHigh);
        netForce = CalcNetForce(velocity);
        if netForce > 0
            velHigh = velocity;
        elseif netForce < 0
            velLow = velocity;
        end

        iterations = iterations+1;
        if iterations > 20
            velocity = -1;
            break;
        end

    end
    
    %plane.performance.dynamicThrust = CalcDynamicThrust(propDiameter, propPitch, motorRPM, velocity, dThrustNeuralNet, dThrustStats); %max thrust available at cruise velocity
    wattFraction = plane.powerSystem.thrust/TotalDrag(velocity); %less power is consumed as motor is no longer able to give as much thrust
    plane.powerSystem.time = plane.powerSystem.time*(0.7/wattFraction); % 30% factor of safety in derating power consumption
    
    if missionNumber == 2
        plane.performance.drag2 = TotalDrag(velocity);
    elseif missionNumber == 3
        plane.performance.drag3 = TotalDrag(velocity);
    end
    %These values are saved so we can review them for reasonableness later
    plane.performance.inducedDrag = InducedDrag(velocity);
    plane.performance.parasiticDrag = Parasitic(velocity);
    plane.performance.skinDrag = Skin(velocity);
    plane.performance.wingPara = wingParasitic(velocity);
    plane.performance.hStabPara = hStabParasitic(velocity);
    plane.performance.vStabPara = vStabParasitic(velocity);
    plane.performance.fusePara = fuselageParasitic(velocity);
    plane.performance.gearPara = gearParaDrag(velocity);
    
    %Parasitic = @(v)  + hStabParasitic(v) + vStabParasitic(v) + fuselageParasitic(v) + gearParaDrag(v);
    
    if velocity<0
        velocity = -1; %DEBUG - make an error flag - functions look for exactly this value for now
    elseif velocity>100 %no way we ever get this fast.
        velocity = -1;
    elseif ((isnan(velocity)) || (velocity == inf)) %if solver doesn't converge
        disp("Found NAN or inf!") %DEBUG - can remove when debugging is finished.
        velocity = -1;
    end
    
    Vstall = sqrt(2*weightTotal/(rho*plane.wing.planformArea*plane.wing.clFlap));
    landingV = 1.3*Vstall; %From raymer textbook
    
    if missionNumber == 2
        plane.performance.velocity2 = velocity;
        plane.performance.landingSpeed2 = landingV;
    elseif missionNumber == 3
        plane.performance.velocity3 = velocity;
        plane.performance.landingSpeed3 = landingV;
    end
end
%-------------------Reference for old variable names from 2020---------------%
%G is vector input from landingGear.m
%Thrust = Thrust.... (lb)
%dp = diameter propeller (in)
%CLw = coeff lift, wing (no dim)
%WeightT = Weight Total (lb)
%rho = density air (slug/ft^3)
%AR = Wing aspect ratio (no dim)
%WingS = wing reference area (ft^2)
%HStabS = horizontal stabilizer reference area (ft^2)
%VStabS = Vertical stabilizer reference area (ft^2)
%CDwing = coeff drag wing (no dim)
%Cw = chord wing (ft)
%CDHstab = coeff horizontal stabilizer (no dim)
%Chstab = chord hoirzontal stab (ft)
%CDVstab = coeff drag vertical stabilizer (no dim)
%Cvstab = chord vertical stab (ft)
%FSA = fuselage frontal surface area (ft^2)
%SAfT = Surface Area fuselage Total (ft^2)
%Lf = length fuselage (ft)
%Wf = width fuselage (in)
%SensorCalc = logical value, 1 or 0. 0 will not run sensor calculations
%ARs = sensor aspect ratio
%Ws = weight sensor
%CLflaps = max coefficient of lift for wing with flaps down
%Outputs:
%V = cruise velocity
%landingV = landing speed

%Lf = plane.fuselage.length; %ft
%radius_wheel = plane.fuselage.wheelRadius; %in
%CLw = plane.wing.clw; %dimensionless
%----------------------------------------------------------------------------%