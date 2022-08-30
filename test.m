%File for testing snippets of code without having to run the whole analysis
%Right now this just contains the code necessary to find the total number of airplanes checked by the analysis with
%the given number of parameters to iterate over.

clear; clc;
warning('off','all')
%rho=0.00235308; %Desnity air at Tuscon, Az (slug/ft^3)
rho = 0.002391; %Density air at Whichita, Ks with average climate data from April 2021
%Aspect_Ratios = 8:.5:15;%the wing aspect ratios being considered 
Aspect_Ratios = [7, 8, 9, 10, 11, 12, 13]; %other runs will have 11, 12, and 13 for four total analyses.
%syringes = 10:100;
syringes = 50:20:300;
load("MotorSpreadsheet.mat");
Num_Power_Systems = height(MotorSpreadsheet);

width_wheel = 0.5;    %width of wheels (in)
radius_wheel = 1.5; %wheel radius
VSAR = 2;
HSAR = 4.5;

span = 8;

[wings] = wingData(Aspect_Ratios, span); %call wing function to make airfoil data lookup table
[wingrow, wingcol, wingpg] = size(wings);

%only analyze naca4415 for new run of analysis.
%max_index = wingpg*wingrow*(Num_Power_Systems)*length(syringes)*length(syringes);
%plane(1:max_index) = struct(airplaneClass);
index = 1;
iterNum = 1;
for AR = 1:wingpg
    for airfoil = [11,15] %DEBUG --just naca and goe airfoils
        for powerIndex = 1:Num_Power_Systems
            for syringe_index = 1:length(syringes)
                for num_vials = 1:10 
                    index = index+1;
                end
            end
        end
    end
end
disp(index)