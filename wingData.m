%%

%this funcion uses the monoplane equation (matrix method) for generating a
%wing

%the function creates a wing for each airfoil and aspect ratio to determine
%the approximate C_l, C_d, and weight for a wing with that configuration.
%The method assumes a rectangular wing (chord length is the same at root
%and tip).
%The code uses the Prandtl-Lifting-Line theory to estimate the general lift
%distribution for the wing.

function [wings]=wingData(Aspect_Ratios, span) 
%Define our dimensional constants
AR=Aspect_Ratios;                 %the wing aspect ratios being considered 
%Define our material constants

%For CF wing
%rho_foam=2.13;              %density of the foam used to make the wing (lb/ft^3)
%w_cf = 0.0450596;           %weight of carbon fiber sheets (lb/ft^2)

%For balsa wing
rho_balsa=0.00579;               %density of aircraft balsa (lb/in^3)
w_spar=0.00958;                 %Weight (lb/in) of the CF spars, based on mcmaster part https://www.mcmaster.com/2153T85/

%no spar. approximation for Ibeam using carbon fiber sheet weight/area
%If we decide to do a spar
%rho_spar=68.3;              %the density of the carbon fiber spar in the wing (lb/ft^3)
%cs_spar=pi*(0.1875/12)^2;   %the cross-section of the carbon fiber spar in the wing
%V_spar=b*cs_spar;           %the volume of the wing spar
%W_spar=V_spar*rho_spar;     %the weight of the wing spar


%the next set of code creates a reference matrix for each airfoil, all data
%taken at alpha=3 degrees:
%af - airfoil name
%cl - coefficient of lift for that airfoil
%cd - coefficient of drag for that airfoil
%alpha0 - zero lift angle of attack for that airfoil
%max_t - the maximum thickness of that airfoil 
%max_tx - the position of maximum thickness for that airfoil
%data - the complete matrix of airfoil information

%Airfoils used, best ones selected from full wing analysis
%S4022 S3010 S2090 S4233 S7055, referred as: 4022 3010 2090 4233 7055
%FX 6-100 FX 6-100 10%, referred as: 6100 6110 
%Eppler 66 Eppler 393 Eppler 392, referred as: 66 393 392 
%goe433 goe332 goe404, referred as: 433 332 404 
%NACA 4415 NACA 1412 NACA 2414 NACA 2410 NACA M18, referred as: 4415 1412 2414 2410 M18 
%Clark Y Clark Z, referred as: 25 26

%Define our airfoil constants parameters, these are independant of span
af=[4022 3010 2091 4233 7055 6100 6110 66 392 393 433 332 404 1412 4415 2414 2410 18 25 26]; 
cl=[.904 .630 .743 .681 .710 .814 .809 .832 .768 .824 0.894 0.966 .810 0.521 0.839 0.646 0.580 0.68 0.735 0.765]; 
cd=[.015 .012 .013 .017 .013 .012 .012 .014 .015 .016 0.018 0.020 .016 0.012 0.015 0.013 0.011 0.017 0.013 0.014]; 
alpha0=(pi/180)*[-3.5 -3 -4 -3.5 -3.25 -3.25 -3.25 -3 -4.2 -4.5 -5 -6.9 -4.9 -0.8 -3.4 -1.8 -1.6 -2.4 -3.1 -3.8]; 
max_t=[.1126 .1032 .1011 .1364 .105 .0999 .0998 .1013 .1015 .1153 .1741 .1179 .1319 .120 .1499 .14 .10 .1202 .1171 .1175]; 
max_tx=[.3403 .2503 .2607 .3385 .3174 .2793 .2783 .3145 .3068 .3233 .2913 .2983 .2963 .2993 .3083 .2953 .2983 .3003 .2803 .3003]; 
data=[af' cl' cd' alpha0' max_t' max_tx'];

alphac=3*pi/180;        %the cruise angle of attack of the airfoil (3 degrees)
alphat=15*pi/180;       %the takeoff angle of attack (15 degrees)
alphae=alphac-alpha0;   %the effective angle of attack for the airfoil (cruise angle minus zero angle)
alphamax=alphat-alpha0; %the effective angle of attack at takeoff (takeoff angle minus 0 angle)

theta=[pi/4 pi/2 3*pi/4 pi];    %the theta values for positions down the wing
n=[1 3 5 7];                    %n values for determining the wing's coefficient of drag

%psi and zeta are an arbitrary matrix and vector, respectively, to solve
%for the coefficients, the A vector, needed to determine the lift and drag
%coefficients for the wing
psi=zeros(4,4);
zeta=zeros(length(theta),1);
zetamax=zeta;

wings=zeros(length(af),10,length(AR),length(span));     %matrix to store wing data

%for loops to run the wing analysis for each aspect ratio and airfoil
for x=1:length(AR)
    for y=1:length(alphae)
        mu=pi/(2*AR(x));
        for b=1:length(span) %new for loop to iterate over span values
        for z=1:length(theta)
            psi(z,:)=[sin(theta(z))*(mu+sin(theta(z))) , sin(3*theta(z))*(3*mu+sin(theta(z))) , sin(5*theta(z))*(5*mu+sin(theta(z))) , sin(7*theta(z))*(7*mu+sin(theta(z)))];
            zeta(z)=mu*alphae(y)*sin(theta(z));
            zetamax(z)=mu*alphamax(y)*sin(theta(z));
        end
        A=psi\zeta;                     %the coefficients used to determine C_l and C_di for the wing
        Amax=psi\zetamax;               %the coefficients used to determine C_lmax (takeoff)
        C_l=A(1)*pi*AR(x);              %C_l for the wing
        C_lmax=Amax(1)*pi*AR(x);        %C_lmax (C_l at takeoff) for the wing
        C_di=pi*AR(x)*dot(n,A.^2);      %C_di for the wing
        chord=span(b)/AR(x);               %chord length of the wing (based on aspect ratio)
        t=chord*max_t(y);               %thickness of the wing
        planArea=chord*span(b);
        
        %Flaps portion: pg. 8-13 https://www.fzt.haw-hamburg.de/pers/Scholz/HOOU/AircraftDesign_8_HighLift.pdf
        deltaC_lmax = (0.95)*(0.58)*(0.28)*1.15;   %I kept the ratios like this for future reference     
        C_lflaps = C_lmax + deltaC_lmax;
        
        %the next set of code determines the weight of the wing, by
        %approximating the cross-sectional area as 2 right triangles, with
        %the position of maximum thickness for the airfoil used to define
        %the bases of the triangles.
        
        %for CF wing
        %ax=(chord*max_tx(y)*t/2)+(chord*(1-max_tx(y))*t/2);             %wing cross-sectional area 
        %W_wing=(ax*span(b))*rho_foam;                                 %weight of the foam used for the wing
        planArea = chord*span(b);
        surfArea = 2*planArea;
        %%W_cf = (2*surfArea*w_cf)+(2*t*span(b)*w_cf);                        %total weight of the carbon fiber sheets, assumes 2 total layers of carbon and I-beam method
        %W_cf = (2*surfArea*w_cf)+(2*t*span(b)*w_cf);
        %W=W_cf+W_wing;                                          %weight of the wing and spar
        %W=W*1.4;                                                %weighting factor based on innacuracy in last year's CF weight estimates
        
        %for balsa wing
        ax=((chord*max_tx(y)*t/2)+(chord*(1-max_tx(y))*t/2))*(12^2);            %wing cross-sectional area (in^2)
        v_balsa=(((span(b)*12)/3)*ax*0.125)+(12*(span(b)*12)*(.125^2));         %volume of balsa in wing including ribs and stringers (in^3)
        W=(rho_balsa*v_balsa)+(2*w_spar*(span(b)*12));                        %weight of balsa+spar (lb)
        W=W*1.3;                                                       %Weighting factor to account for epoxy and monocoat 
        
        %the next 4 lines assign the wing data to positions in the wing
        %matrix. The first column is the wing C_l, the second is the wing
        %C_lmax, the third C_di, and fourth weight. Each row corresponds to 
        %an airfoil in the NACA vector, and each page corresponds to an 
        %aspect ratio. E.g., the first airfoil in the NACA vector is the 
        %1408, and the first AR is 8. So the first row of the first page in 
        %'wings' is the C_l, C_di, and weight for a wing with a NACA 1408 
        %airfoil, and an aspect ratio of 8.
        %The fifth column is the chord
        %the 6th column is the wing planform area.

        %The fifth column is the name of the airfoil (divided by 1000 so
        %that the other data is visible), for reference.
        
        %define # of airfoils. new wing = these values, wing 
        %syntax wing(row - airfoil, column - parameters , page - aspect ratio, dimension - span)
        wings(y,1,x,b)=C_l;    %s denotes the dimension associated with span (new)                                  
        wings(y,2,x,b)=C_lmax;
        wings(y,3,x,b)=C_di;
        wings(y,4,x,b)=C_lflaps;
        wings(y,5,x,b)=W;
        wings(y,6,x,b)=chord;
        wings(y,7,x,b)=planArea;
        wings(y,8,x,b)=surfArea;
        wings(:,9,x,b)=af'./1000;
        wings(y,10,x,b)=t;
        end
    end
end

end