function [FSA,fh,fw,fl,Weightfuse,SAfT] = SaestsNew(ns) %surface area estimate for fuselage
%Inputs: sensor aspect ratio (2020), number of containers
%Akaash  - For 2021, input will just change to number of syringes, which will be ns
%DEBUG - need new version of this - just a look up table. how will fuselage grow
sw = ; %width syringe
sh = ; %height syringe
sl= ; %length syringe
nc= ; %number of crates
cbl= ; %conveyor belt length
fh = 0; %fuse height
fw = 0; %fuse width
fl = 0; %fuse length
rhoCarbon=120.486; %lb/ft^3
if ns<=9
    nc=0;
    cbl=0;
    fh=0;
    fw=0;
    fl=0; %make all of these 0 because then 
elseif nc<=19
    nc= 
    cbl=
    fh=
    fw=
    fl=
elseif nc<=29
    nc= 
    cbl=
    fh=
    fw=
    fl=
elseif nc<=39
    nc=
    cbl=
    fh=
    fw=
    fl=
elseif nc<=49
    nc=
    cbl=
    fh=
    fw=
    fl=
elseif nc<=59
    nc=
    cbl=
    fh=
    fw=
    fl=
end

%need to make sure that the length of the fueslage is less than 8 feet, so the if condition below will check to see if it is below 8

if fl>=96 
   fl=0;
   fh=0;
   fw=0; %make all of these 0 bc we can't use it 
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