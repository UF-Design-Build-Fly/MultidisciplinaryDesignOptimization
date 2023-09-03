%%

%this funcion uses the monoplane equation (matrix method) for generating a
%wing

%the function creates a wing for each airfoil and aspect ratio to determine
%the approximate C_l, C_d, and weight for a wing with that configuration.
%The method assumes a rectangular wing (chord length is the same at root
%and tip).
%The code uses the Prandtl-Lifting-Line theory to estimate the general lift
%distribution for the wing.
% clc;clear;

function [wings]=wing(AR) %DEBUG - all the orange comes from this function

b=5;                        %span of the wing %DEBUG if wingspan is changing this needs to input
rho_spar=68.3;              %the density of the carbon fiber spar in the wing (lb/ft^3)
Ax_spar=pi*(.1875/12)^2;    %the cross-section of the carbon fiber spar in the wing
V_spar=b*Ax_spar;           %the volume of the wing spar
W_spar=V_spar*rho_spar;     %the eight of the wing spar

rho_foam=2.13;              %density of the foam used to make the wing (lb/ft^3)

%the next set of code creates a reference matrix for each airfoil, all data
%taken at alpha=3 degrees: %DEBUG - talk to jasmin about the alpha
%NACA - airfoil name
%cl - coefficient of lift for that airfoil
%cd - coefficient of drag for that airfoil
%alpha0 - zero lift angle of attack for that airfoil
%max_t - the maximum thickness of that airfoil 
%max_tx - the position of maximum thickness for that airfoil
%data - the complete matrix of airfoil information

%63210 - 63A210
%8036 - s8036
%8037 - s8037
%2091 - s2091
%3002 - s3002
%67 - e67
%66 - e66
%0 - clarkY
%7062 - sd7062
%3012 - hq 3.0/12
NACA=[1408 1410 1412 2408 2410 2412 2415 2418 2421 2424 4412 4415 4418 4421 6409 6412 63210 8036 8037 2091 3002 67 66 0 7062 3012];
cl=[.3377 .4082 .4472 .4445 .4688 .4528 .4065 .3356 .2726 .0542 .5279 .3020 .2943 .0590 .5997 .5632 .3963 .1968 .2186 .6095 .5698 .3545 .4302 .5275 .5115 .3807];
cd=[.02271 .02562 .02791 .02628 .02950 .03304 .03941 .04769 .06115 .08159 .04426 .05278 .06234 .07671 .04926 .05831 .03095 .04811 .05268 .03799 .0342 .04189 .03963 .03775 .04793 .04001];
alpha0=(pi/180)*[-.25 .125 .9 -.75 0 .25 -.25 -.125 .875 2.75 -1.25 -.25 -.25 1.375 -1.25 -2 -3 -1.625 -2.125 -1.25 -.5 .25 -.525 -.2805 -1.75 .625];
max_t=[.08 .1 .12 .08 .1 .12 .15 .18 .21 .24 .12 .15 .18 .21 .09 .12 .1 .16 .16 .101 .099 .116 .101 .117 .14 .12];
max_tx=[.3 .299 .299 .299 .299 .3 .3 .3 .3 .297 .3 .309 .3 .3 .293 .301 .349 .368 .335 .26 .306 .326 .287 .28 .255 .35];
data=[NACA' cl' cd' alpha0' max_t' max_tx'];


alphac=3*pi/180;        %the cruise angle of attack of the airfoil (3 degrees)
alphat=15*pi/180;       %the takeoff angle of attack (10 degrees)
alphae=alphac-alpha0;   %the effective angle of attack for the airfoil (cruise angle minus zero angle)
alphamax=alphat-alpha0; %the effective angle of attack at takeoff (takeoff angle minus 0 angle)

%DEBUG - make this an input instead of internal so AR's aren't hard-coded

theta=[pi/4 pi/2 3*pi/4 pi];    %the theta values for positions down the wing
n=[1 3 5 7];                    %n values for determining the wing's coefficient of drag

%psi and zeta are an arbitrary matrix and vector, respectively, to solve
%for the coefficients, the A vector, needed to determine the lift and drag
%coefficients for the wing
psi=zeros(4,4);
zeta=zeros(length(theta),1);
zetamax=zeta;

wings=zeros(length(NACA),3,length(AR));     %matrix to store wing data

%for loops to run the wing analysis for each aspect ratio and airfoil
for x=1:length(AR)
    for y=1:length(alphae)
    mu=pi/(2*AR(x));
        for z=1:length(theta)
            psi(z,:)=[sin(theta(z))*(mu+sin(theta(z))) , sin(3*theta(z))*(3*mu+sin(theta(z))) , sin(5*theta(z))*(5*mu+sin(theta(z))) , sin(7*theta(z))*(7*mu+sin(theta(z)))];
            zeta(z)=mu*alphae(y)*sin(theta(z));
            zetamax(z)=mu*alphamax(y)*sin(theta(z));
        end
        A=psi\zeta;                     %the coefficients used to determine C_l and C_di for the wing
        Amax=psi\zetamax;               %the coefficients used to determine C_lmax (takeoff)
        C_l=A(1)*pi*AR(x);              %C_l for the wing
        C_lmax=Amax(1)*pi*AR(x);      %C_lmax (C_l at takeoff) for the wing
        C_di=pi*AR(x)*dot(n,A.^2);      %C_di for the wing
        c=b/AR(x);                      %chord length of the wing (based on aspect ratio)
        t=c*max_t(y);                   %thickness of the wing
        
        %the next set of code determines the weight of the wing, by
        %approximating the cross-sectional area as 2 right triangles, with
        %the position of maximum thickness for the airfoil used to define
        %the bases of the triangles.
        ax=(c*max_tx(y)*t/2)+(c*(1-max_tx(y))*t/2)-Ax_spar;     %wing cross-sectional area, minus the area taken up by the spar 
        W_wing=(ax*b)*rho_foam;                                 %weight of the foam used for the wing
        W=W_spar+W_wing;                                        %weight of the wing and spar
        
        %the next 4 lines assign the wing data to positions in the wing
        %matrix. The first column is the wing C_l, the second is the wing
        %C_lmax, the third C_di, and fourth weight. Each row corresponds to 
        %an airfoil in the NACA vector, and each page corresponds to an 
        %aspect ratio. E.g., the first airfoil in the NACA vector is the 
        %1408, and the first AR is 8. So the first row of the first page in 
        %'wings' is the C_l, C_di, and weight for a wing with a NACA 1408 
        %airfoil, and an aspect ratio of 8.
        %The fifth column is the name of the airfoil (divided by 1000 so
        %that the other data is visible), for reference.
        wings(y,1,x)=C_l;                                       
        wings(y,2,x)=C_lmax;
        wings(y,3,x)=C_di+cd(y);
        wings(y,4,x)=W;
        wings(:,5,x)=NACA'./1000; 
    end
end
end
