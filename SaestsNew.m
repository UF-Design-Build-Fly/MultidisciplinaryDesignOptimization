function [FSA,fh,fw,fl,Weightfuse,SAfT] = SaestsNew(AR,nc) %surface area estimate for fuselage
%DEBUG - need new version of this - just a look up table. how will fuselage
%grow
cw = 1.5; %width sensor
ch = 1.5; %height sensor
cl=cw*AR; %length sensor
fh = 0; %fuse height
fw = 0; %fuse width
fl = 0; %fuse length
rhoCarbon=120.486; %lb/ft^3
if nc<=4
    fh=cw;
    fw=cw;
    fl=nc*cw*AR;
elseif nc<=8
    fh=2*cw;
    fw=cw;
    fl=4*cl;
elseif nc<=16
    fh=2*cw;
    fw=2*cw;
    fl=4*cl;
elseif nc<=24
    fh=3*cw;
    fw=2*cw;
    fl=4*cl;
elseif nc<= 36
    fh=3*cw;
    fw=3*cw;
    fl=4*cl;
elseif nc<=48
   fh=4*cw;
   fw=3*cw;
   fl=4*cl;
elseif nc<=64
    fh=4*cw;
    fw=4*cw;
    fl=4*cl;
end

SAbrick=2*fl*fh+2*fl*fw;
Hf=sqrt( (0.5*fh)^2 +6^2);
SAfTri=2*Hf*fw + fh*6;
Hb=sqrt( fh^2 +8^2);
SAbTri=fw*8 + fw*Hb + fh*8;
SAfT=SAbrick+SAfTri+ SAbTri;
SAfT=SAfT/144;%Surface Area fuse Total,ft^2
FSA=fw*fh/144;%Frontal Surface Area,ft^2
Vol=SAfT*(1/16)/12;%The volume in ft^3
Weightfuse=Vol *rhoCarbon;
end