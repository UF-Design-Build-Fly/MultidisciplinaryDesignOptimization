function plane = GenVelocityTest(plane, missionNumber, rho, tempK)
    % Select mission mass
    if missionNumber == 2, Wkg = plane.performance.totalWeight2;
    else,                   Wkg = plane.performance.totalWeight3;
    end
    WN = Wkg * 9.80665;

    % Induced
    e = max(0.6, 1.78*(1-0.045*plane.wing.aspectRatio^0.68)-0.64);
    S = max(1e-6, plane.wing.planformArea);
    b = max(1e-6, plane.wing.span);
    Di = @(V) (2*WN^2) ./ (rho*(V.^2)*pi*e*b^2);

    % Parasitic
    Cd_wing = max(0, plane.wing.cd);
    wingPara = @(V) 0.5*rho*(V.^2)*S*Cd_wing;
    hPara    = @(V) 0.5*rho*(V.^2)*max(0,plane.empennage.HSarea) * max(0,plane.empennage.HScd);
    vPara    = @(V) 0.5*rho*(V.^2)*max(0,plane.empennage.VSarea) * max(0,plane.empennage.VScd);
    fusPara  = @(V) 0.5*rho*(V.^2)*max(0,plane.fuselage.frontalSurfaceArea) * 0.75;
    gearPara = @(V) 0.5*rho*(V.^2)*max(0,plane.fuselage.gearFrontalSA)  * 0.05;
    wheelPara= @(V) 0.5*rho*(V.^2)*max(0,plane.fuselage.wheelFrontalSA) * 0.245;

    % Skin friction
    mu = 15.97e-6; a = sqrt(1.4*287*max(150,tempK));
    Cfl = @(V,L) 1.328./sqrt(max(V.*L/mu,1e-6));
    Cft = @(V,L) 0.455./((log10(max(V.*L/mu,1.1)).^2.58) .* (1+0.144*(V./a).^2).^0.65);
    Cf  = @(V,L,k) k.*Cfl(V,L) + (1-k).*Cft(V,L);
    Skin = @(V) 2*(0.5*rho*(V.^2).*Cf(V,max(plane.wing.chord,1e-3),1).*S) + ...
                2*(0.5*rho*(V.^2).*Cf(V,max(plane.empennage.HSchord,1e-3),1).*max(0,plane.empennage.HSarea)) + ...
                2*(0.5*rho*(V.^2).*Cf(V,max(plane.empennage.VSchord,1e-3),1).*max(0,plane.empennage.VSarea)) + ...
                  (0.5*rho*(V.^2).*Cf(V,max(plane.fuselage.length,1e-3),0.9).*max(0,plane.fuselage.totalSA));

    % Banner drag (Mission 3 only) — vertical flat plate
    if missionNumber == 3
        Lin = plane.performance.bannerLength_in * 0.0254;
        Hin = plane.performance.bannerHeight_in * 0.0254;
        Hin = max(Hin, max(0.0508, Lin/5));      % H >= 2" and L/H <= 5
        A_banner = max(1e-4, Lin * Hin);
        Cd_banner = 1.2;
        BannerDrag = @(V) 0.5*rho*(V.^2)*Cd_banner*A_banner;
    else
        BannerDrag = @(V) 0;
    end

    TotalDrag = @(V) Di(V) + wingPara(V) + hPara(V) + vPara(V) + fusPara(V) + gearPara(V) + wheelPara(V) + Skin(V) + BannerDrag(V);

    % Thrust model
    T_static = max(0, plane.powerSystem.thrust) * 9.80665;  % kgf->N
    V_pitch  = max(15, plane.powerSystem.propSpeed);
    Thrust   = @(V) max(0, T_static .* max(0, 1 - V./V_pitch));

    % Solve Thrust(V) = Drag(V) by bisection
    V_lo=5; V_hi=max(20, min(120, V_pitch*0.99)); f=@(V) Thrust(V)-TotalDrag(V);
    if f(V_lo) < 0, Vstar=-1;
    else
        it=0; Vstar=-1;
        while it<40
            Vm=0.5*(V_lo+V_hi); val=f(Vm);
            if abs(val)<0.05, Vstar=Vm; break; end
            if val>0, V_lo=Vm; else, V_hi=Vm; end
            it=it+1;
        end
    end
    if Vstar<0 || Vstar>120 || ~isfinite(Vstar), Vstar=-1; end

    % Store
    if missionNumber==2, plane.performance.drag2 = TotalDrag(max(Vstar,1));
    else,                plane.performance.drag3 = TotalDrag(max(Vstar,1));
    end

    % Power/time scaling
    if Vstar>0
        D = TotalDrag(Vstar);
        frac = min(1.0, max(0.2, D / max(T_static, 1e-3)));
        plane.performance.dynamicThrust = T_static * frac;
        if isfinite(plane.powerSystem.time) && plane.powerSystem.time>0
            plane.powerSystem.time = plane.powerSystem.time * (0.7/frac);
        end
    end

    % Stall/landing
    CLmax = max(plane.wing.clFlap, plane.wing.clm);
    Vstall = sqrt(2*WN/(rho*S*max(1e-6,CLmax))); landingV = 1.3*Vstall;
    if missionNumber==2
        plane.performance.velocity2 = Vstar;
        plane.performance.landingSpeed2 = landingV;
    else
        plane.performance.velocity3 = Vstar;
        plane.performance.landingSpeed3 = landingV;
    end

    % Breakdown snapshot
    if Vstar>0
        plane.performance.inducedDrag=Di(Vstar);
        plane.performance.parasiticDrag=wingPara(Vstar)+hPara(Vstar)+vPara(Vstar)+fusPara(Vstar)+gearPara(Vstar)+wheelPara(Vstar);
        plane.performance.skinDrag=Skin(Vstar);
        plane.performance.wingPara=0.5*rho*(Vstar.^2)*S*plane.wing.cd;
        plane.performance.hStabPara=0.5*rho*(Vstar.^2)*max(0,plane.empennage.HSarea)*max(0,plane.empennage.HScd);
        plane.performance.vStabPara=0.5*rho*(Vstar.^2)*max(0,plane.empennage.VSarea)*max(0,plane.empennage.VScd);
        plane.performance.fusePara=0.5*rho*(Vstar.^2)*max(0,plane.fuselage.frontalSurfaceArea)*0.75;
        plane.performance.gearPara=0.5*rho*(Vstar.^2)*max(0,plane.fuselage.gearFrontalSA)*0.05 + ...
                                   0.5*rho*(Vstar.^2)*max(0,plane.fuselage.wheelFrontalSA)*0.245;
    end
end

