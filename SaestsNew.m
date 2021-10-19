function [FSA,fh,fw,fl,Weightfuse,SAfT] = SaestsNew(ns) %surface area estimate for fuselage
%Inputs: sensor aspect ratio (2020), number of containers
%Akaash  - For 2021, input will just change to number of syringes, which will be ns
%DEBUG - need new version of this - just a look up table. how will fuselage grow

%Akaash: How I'm thinking about putting the syringes is in groups of five
%on the fuselage and then just stacking these groups of five on top of each other

sw=1.25; %width syringe (in)
sh=1.75; %height syringe (in)
sl=6; %length syringe (in)
nc=0; %number of crates
cl=3.5; %crate length (in) 
ch=cl; %crate height
rml=5; %release mechanism length (in)
rmar=2; %release mechanism aspect ratio <---- Just made it 2 for now, can change later
fh=0; %fuse height (in)
fw=0; %fuse width (in)
fl=0; %fuse length (in)
rhoCarbon=120.486; %lb/ft^3

if ns<=5 %this defines the function for the number of syringes we will carry 
    nc=0; %can only carry 1 crate if we carry 10 syringes
    rml=nc*cl+rml; %this will just be the length of the release mechanism so that we take it into account since we are dropping crates out the back of the fuselage
    fh=ch*rmar; %this will just be crate height x a ratio we want to ensure they are easily able to fit and also take into account mechanim height. 
    fw=sw*5; %I'm thinking of laying out the syringes in rows of five, so this is why this is sw*5
    if sh<fh %this if will be used to extend the length of the fuselage accordingly depending on how many syringes we will carry
        fl=rml+sl; 
    elseif sh>=fh %if the syringe height after stacking them exceeds the fuselage height, we have to extend the fueslage so we have more space to put them
        a = sh/fh;
        a = ceil(a); %this a variable will just 
        fl=a*sl+rml;
    end
elseif ns<=10
    nc=1;
    rml=nc*cl+rml;  
    fh=ch*rmar; 
    fw=sw*5; 
    sh=2*sh;
    if sh<fh
        fl=rml+sl; 
    elseif sh>=fh
        a = sh/fh;
        a = ceil(a);
        fl=a*sl+rml;
    end
elseif ns<=15
    nc=1;
    rml=nc*cl+rml;  
    fh=ch*rmar; 
    fw=sw*5; 
    sh=3*sh;
    if sh<fh
        fl=rml+sl; 
    elseif sh>=fh
        a = sh/fh;
        a = ceil(a);
        fl=a*sl+rml;
    end
elseif ns<=20
    nc=2;
    rml=nc*cl+rml+1.25;  
    fh=ch*rmar; 
    fw=sw*5; 
    sh=4*sh;
    if sh<fh
        fl=rml+sl; 
    elseif sh>=fh
        a = sh/fh;
        a = ceil(a);
        fl=a*sl+rml;
    end
elseif ns<=25
    nc=2;
    rml=nc*cl+rml+1.25;  
    fh=ch*rmar; 
    fw=sw*5; 
    sh=5*sh;
    if sh<fh
        fl=rml+sl; 
    elseif sh>=fh
        a = sh/fh;
        a = ceil(a);
        fl=a*sl+rml;
    end
elseif ns<=30
    nc=3;
    rml=nc*cl+rml+2.5;  
    fh=ch*rmar; 
    fw=sw*5; 
    sh=6*sh;
    if sh<fh
        fl=rml+sl; 
    elseif sh>=fh
        a = sh/fh;
        a = ceil(a);
        fl=a*sl+rml;
    end
elseif ns<=35
    nc=3;
    rml=nc*cl+rml+2.5;  
    fh=ch*rmar; 
    fw=sw*5; 
    sh=7*sh;
    if sh<fh
        fl=rml+sl; 
    elseif sh>=fh
        a = sh/fh;
        a = ceil(a);
        fl=a*sl+rml;
    end
elseif ns<=40
    nc=4;
    rml=nc*cl+rml+3.75;  
    fh=ch*rmar; 
    fw=sw*5; 
    sh=8*sh;
    if sh<fh
        fl=rml+sl; 
    elseif sh>=fh
        a = sh/fh;
        a = ceil(a);
        fl=a*sl+rml;
    end
elseif ns<=45
    nc=4;
    rml=nc*cl+rml+3.75;  
    fh=ch*rmar; 
    fw=sw*5; 
    sh=9*sh;
    if sh<fh
        fl=rml+sl; 
    elseif sh>=fh
        a = sh/fh;
        a = ceil(a);
        fl=a*sl+rml;
    end
elseif ns<=50
    nc=5;
    rml=nc*cl+rml+5;  
    fh=ch*rmar; 
    fw=sw*5; 
    sh=10*sh;
    if sh<fh
        fl=rml+sl; 
    elseif sh>=fh
        a = sh/fh;
        a = ceil(a);
        fl=a*sl+rml;
    end
elseif ns<=55
    nc=5;
    rml=nc*cl+rml+5;  
    fh=ch*rmar; 
    fw=sw*5; 
    sh=11*sh;
    if sh<fh
        fl=rml+sl; 
    elseif sh>=fh
        a = sh/fh;
        a = ceil(a);
        fl=a*sl+rml;
    end
end

%Akaash: I still do not understand these equations that were made last year
%and do not know if I should make new functions or not because I am not
%sure if these functions account for the extra space that will be needed
%in the fuselage for systems things. Honestly just need clarification on
%what these equations are really trying to say
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