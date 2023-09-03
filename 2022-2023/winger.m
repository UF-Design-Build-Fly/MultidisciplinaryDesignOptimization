function [T] =  winger(af, AR)
%wing-indexer (winger) function
%Inputs:
%airfoil (af) index (1-20)
%aspect ratio (AR) index (1-19)
%Output:
%Table with Cl, Clm, Cd0, W, airfoil name
%Note: airfoil names referred to as follows:
%S4022 S3010 S2090 S4233 S7055, referred as: 4022 3010 2090 4233 7055
%FX 6-100 FX 6-100 10%, referred as: 6100 6110 
%Eppler 66 Eppler 393 Eppler 392, referred as: 66 393 392 
%goe433 goe332 goe404, referred as: 433 332 404 
%NACA 4415 NACA 1412 NACA 2414 NACA 2410 NACA M18, referred as: 4415 1412 2414 2410 M18 
%Clark Y Clark Z, referred as: 25 26
[wings]=wing;
cl = wings(af, 1, AR);      %cl
clm = wings(af, 2, AR);     %cl max
cd0 = wings(af, 3, AR);     %cd i zero velocity coefficient of drag
w = wings(af, 4, AR);       %weight
airfoil = wings(af, 5, AR);    %airfoil name
T = table(cl,clm,cd0,w,airfoil,'VariableNames',{'cl','clm','cd0','w','airfoil'});
end


        