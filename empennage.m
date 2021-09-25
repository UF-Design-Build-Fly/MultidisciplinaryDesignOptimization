function [HS_area, HS_chord, HS_span, HS_weight, HS_Cd, VS_area, VS_chord, VS_span, VS_weight, VS_Cd] = empennage(wingarea, HSAR, VSAR)
    %DEBUG - make sure this is being output in correct order. This does not
    %solve via the monoplane equation. probably should change to monoplane
    %for higher accuracy even if it's a bit computationally slower
    rho_XPS = 2.3;   %lb/ft^3
    e = 0.7;  %efficiency factor
    alpha = 3.*pi/180;  %assume 3 degrees
    Cd0 = 0.0165;
    HS_area = wingarea/4;
    HS_chord = sqrt(HS_area/HSAR);
    HS_span = HS_area/HS_chord;
    airfoil = @(x) .6*(0.2969.*sqrt(x)-0.1260.*x-0.3516.*x.^2+0.2843.*x.^3-0.1015.*x.^4);  %NACA0012 airfoil is just assumed
    CSH_area = 2*integral(airfoil, 0, HS_chord); %cross sectional area of airfoil horizontal stab
    HS_weight = CSH_area*HS_span*rho_XPS;   %lb
    HS_Cdi = 4.*pi.*alpha.^2/(HSAR.*e);  %induced drag equation (NASA)
    HS_Cd = Cd0+HS_Cdi;
    VS_Cd = Cd0;     % no alpha, no Cdi
    VS_area = HS_area/2;
    VS_chord = sqrt(VS_area/VSAR);
    VS_span = VS_area/VS_chord;
    CSV_area = 2*integral(airfoil,0,VS_chord);   % cross sectional area of airfoil vertical stab
    VS_weight = CSV_area*VS_span*rho_XPS;    %lb
end