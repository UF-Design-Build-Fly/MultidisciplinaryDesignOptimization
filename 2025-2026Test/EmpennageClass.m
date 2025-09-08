classdef EmpennageClass
  properties
    HSarea=-1; HSchord=-1; HSweight=-1; HScd=-1;
    VSarea=-1; VSchord=-1; VSweight=-1; VScd=-1;
  end
  methods (Static)
    function emp = GenEmpennage(plane, HSAR, VSAR)
        rho_XPS = 36.84; w_cf = 0.00190; Cd0 = 1.328*2.25/(sqrt(150000));
        Sw = max(1e-6, plane.wing.planformArea);
        HS_area = 0.25*Sw; HS_span = sqrt(HS_area*HSAR); HS_chord = HS_area/HS_span;
        t_c = 0.10; c = HS_chord; fy=@(x) 2.5*t_c*(0.2969*sqrt(x/c)-0.1260*(x/c)-0.3516*(x/c).^2+0.2843*(x/c).^3-0.1015*(x/c).^4);
        CSH = 2*integral(fy,0,c); HS_wfoam = CSH*HS_span*rho_XPS; HS_weight = (2.25*HS_area*w_cf*2)+HS_wfoam;
        HS_Cd = Cd0;

        VS_area = 0.5*HS_area; VS_span = sqrt(VS_area*VSAR); VS_chord = VS_area/VS_span;
        c2=VS_chord; fy2=@(x) 2.5*t_c*(0.2969*sqrt(x/c2)-0.1260*(x/c2)-0.3516*(x/c2).^2+0.2843*(x/c2).^3-0.1015*(x/c2).^4);
        CSV = 2*integral(fy2,0,c2); VS_wfoam = CSV*VS_span*rho_XPS; VS_weight = (2.25*VS_area*w_cf*2)+VS_wfoam;
        VS_Cd = Cd0;

        emp = plane.empennage;
        emp.HSarea=HS_area; emp.HSchord=HS_chord; emp.HSweight=HS_weight; emp.HScd=HS_Cd;
        emp.VSarea=VS_area; emp.VSchord=VS_chord; emp.VSweight=VS_weight; emp.VScd=VS_Cd;
    end
  end
end
