classdef EmpennageClass

  properties
    HSarea = -1;
    HSchord = -1;
    HSweight = -1;
    HScd = -1;
    VSarea = -1;
    VSchord = -1;
    VSweight = -1;
    VScd = -1;
  end

  methods (Static)

    function empennage = GenEmpennage(plane, HSAR, VSAR)
        empennage = plane.empennage;
        %solve via the monoplane equation. probably should change to monoplane
        %for higher accuracy even if it's a bit computationally slower
        %Define constants 
        rho_XPS = 2.3;      %density of foam (lb/ft^3)
        w_cf = 0.0450596;   %weight of carbon fiber sheets (lb/ft^2)
        e = 0.7;            %efficiency factor
        alpha = 3.*pi/180;  %assume 3 degrees
        Cd0 = 1.328*2.25/(sqrt(150000));  %zero-lift drag estimation (2.25 is estimate for Swet/Sref ratio)
        %Horizontal stab dimension calculations
        wingarea = plane.wing.planformArea;
        HS_area = wingarea*.25;              
        HS_span = sqrt(HS_area*HSAR);  
        HS_chord = HS_area/HS_span; 
        %Select NACA 0010 airfoil (t/c is 90% of predicted wing t/c)
        airfoil = @(x) 2.5*0.10*(0.2969.*sqrt(x)-0.1260.*x-0.3516.*x.^2+0.2843.*x.^3-0.1015.*x.^4);  %NACA0012 airfoil is just assumed
        %arcL = @(x) sqrt(1+ (2.5*0.10*(-0.5*0.2969.*x.^(-1.5)-0.1260-2*0.3516.*x+3*0.2843.*x.^2-4*0.1015.*x.^3)).^2); 
        CSH_area = 2*integral(airfoil, 0, HS_chord);      %cross sectional area of airfoil horizontal stab
        HS_wfoam = CSH_area*HS_span*rho_XPS;              %lb
        %HS_ArcL = 2*integral(arcL,0,HS_chord);            %arc length
        HS_weight = (2.25*HS_area*w_cf*2)+HS_wfoam;    %total weight of the carbon fiber sheets, (arc length*weight/ft*2 layers);
        
        HS_Cdi = 4.*pi.*alpha.^2/(HSAR.*e);          %induced drag equation (NASA)
        HS_Cd = Cd0+HS_Cdi;
        
        %Vertical stab dimension calculations           
        VS_area = HS_area*0.5;
        VS_span = sqrt(VS_area*VSAR);
        VS_chord = VS_area/VS_span;
        %Select NACA 0010 airfoil (t/c is 90% of predicted wing t/c)           
        CSV_area = 2*integral(airfoil,0,VS_chord);        % cross sectional area of airfoil vertical stab
        VS_wfoam = CSH_area*VS_span*rho_XPS;              %lb
        %VS_ArcL = 2*integral(arcL,0, VS_chord);       %arc length
        VS_weight = (2.25*VS_area*w_cf*2)+VS_wfoam;    %total weight of the carbon fiber sheets, (arc length*weight/ft*2 layers);
        VS_Cd = Cd0;       % no alpha, no Cdi
        
        %Determine Cdi using using Lifting Line theory
        AR = [HSAR, VSAR];
        alpha0=(pi/180)*.7;                  %for NACA 0010
        alphac=1*pi/180;
        alphat = 9*pi/180;
        alphae=alphac-alpha0;
        alphamax=alphat-alpha0;
        theta=[pi/4 pi/2 3*pi/4 pi];    %the theta values for positions down the wing
        n=[1 3 5 7];                    %n values for determining the wing's coefficient of drag
        %psi and zeta are an arbitrary matrix and vector, respectively, to solve
        %for the coefficients, the A vector, needed to determine the lift and drag
        %coefficients for the wing
        psi=zeros(4,4);
        zeta=zeros(length(theta),1);
        zetamax=zeta;
        for ii = 1:1:length(AR)
            mu=pi/(2*AR(ii));
            for z= 1:length(theta)
                psi(z,:)=[sin(theta(z))*(mu+sin(theta(z))) , sin(3*theta(z))*(3*mu+sin(theta(z))) , sin(5*theta(z))*(5*mu+sin(theta(z))) , sin(7*theta(z))*(7*mu+sin(theta(z)))];
                zeta(z)=mu*alphae*sin(theta(z));
                zetamax(z)=mu*alphamax*sin(theta(z));
            end
            A=psi\zeta;                     %the coefficients used to determine C_l and C_di for the wing
            Amax=psi\zetamax;               %the coefficients used to determine C_lmax (takeoff)
            C_di(ii)=pi*AR(ii)*dot(n,A.^2); %C_di for the wing
        end
        %HS_Cdi = 4.*pi.*alpha.^2/(HSAR.*e); %induced drag equation (NASA)
        %HS_Cd = Cd0+HS_Cdi;
        HS_Cd = Cd0+C_di(1);
        VS_Cd = Cd0;       % no alpha, no Cdi
        
        %convert output variables to object notation
        empennage.HSarea = HS_area;
        empennage.HSchord = HS_chord;
        empennage.HSweight = HS_weight;
        empennage.HScd = HS_Cd;
        empennage.VSarea = VS_area;
        empennage.VSchord = VS_chord;
        empennage.VSweight = VS_weight;
        empennage.VScd = VS_Cd;
        
    end

  end
end