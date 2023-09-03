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
%SD7062-14% SD8040-10%, referred as: 7062 8040
%Eppler 66 Eppler 374, referred as: 66 374 
%goe412 goe433, referred as: 412 433
%NACA 2414, referred as: 2414
%HQ 3.0/12, referred as: 3

%Define our airfoil constants parameters, these are independant of span
af=[7062 8040 66 374 412 433 2414 3]; 
cl=[.77485 .64702 .86236 .63729 .86640 .91441 .64901 .71923]; 
cd=[.01249 .00997 .01106 .01185 .01340 .01548 .01163 .71923]; 
alpha0=(pi/180)*[-4.235 -2.714 -3.589 -2.012 -5.805 -5.184 -2.123 -3.591]; 
max_t=[.1398 .1000 .1013 .1091 .1313 .1741 .1400 .1198]; 
max_tx=[.2715 .2933 .3145 .3434 .2973 .2913 .2953 .3504]; 
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