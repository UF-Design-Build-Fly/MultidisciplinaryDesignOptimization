clear; clc;
warning('off','all')
tic
%rho=0.00235308; %Desnity air at Tuscon, Az (slug/ft^3)
rho = 0.002391; %Density air at Whichita, Ks with average climate data from April 2021
%Aspect_Ratios = 8:.5:15;%the wing aspect ratios being considered 
Aspect_Ratios = [7, 8, 9, 10, 11, 12, 13]; %other runs will have 11, 12, and 13 for four total analyses.
%syringes = 10:100;
syringes = 50:10:300;
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
max_index = wingpg*wingrow*(Num_Power_Systems)*length(syringes)*length(syringes);
plane(1:max_index) = struct(airplaneClass);
index = 1;
iterNum = 1;
for AR = 1:wingpg
    for airfoil = [11,15] %DEBUG --just naca and goe airfoils
        for powerIndex = 1:Num_Power_Systems
            for syringe_index = 1:length(syringes)
                for num_vials = 1:floor(syringes(syringe_index)/10)
                    
                    %plane(index) = struct(airplaneClass);
                    plane(index).fuselage.wheelSA = (2*pi*(radius_wheel)^2+ pi*2*radius_wheel*width_wheel)/144;
                    %read values from wings matrix into aircraft properties
                    plane(index).wing.clw = wings(airfoil, 1, AR);
                    plane(index).wing.clm = wings(airfoil, 2, AR);     %cl max
                    plane(index).wing.cd = wings(airfoil, 3, AR);     %cd i zero velocity coefficient of drag
                    plane(index).wing.clFlap = wings(airfoil, 4, AR);  %weight
                    plane(index).wing.weight = wings(airfoil, 5, AR);    %airfoil name
                    plane(index).wing.chord = wings(airfoil, 6, AR);    %airfoil name
                    plane(index).wing.planformArea = wings(airfoil, 7, AR);    %airfoil name
                    plane(index).wing.surfaceArea = wings(airfoil, 8, AR);    %airfoil name
                    plane(index).wing.name = wings(airfoil, 9, AR);    %airfoil name
                    plane(index).wing.aspectRatio = Aspect_Ratios(AR); %ratio between length and width of wing.
                    %load values from power system table into aircraft
                    %properties
                            
                    plane(index) = powerSelections(plane(index), MotorSpreadsheet, powerIndex);
                    
                    %set payload and fuselage configuration
                    plane(index).fuselage.numSyringes = syringes(syringe_index);
                    plane(index).fuselage.numVials = num_vials;
                    plane(index) = fuselage(plane(index)); %calculate fuselage size based on number of vials and syringes
                    plane(index) = landingGear(plane(index), rho);
                    plane(index) = empennage(plane(index), HSAR, VSAR);
                    
                    plane(index) = findTotalWeight(plane(index));
                    
                    %simulate mission 2
                    plane(index) = GenVelocityTest(plane(index), 2, rho); %2 signifies mission 2 configuration
                    plane(index) = TakeoffChecker(plane(index), 2, rho);
                    plane(index) = mission2score(plane(index));
                    
                    %simulate mission 3
                    plane(index) = GenVelocityTest(plane(index), 3, rho); %2 signifies mission 2 configuration
                    plane(index) = TakeoffChecker(plane(index), 3, rho);
                    plane(index) = mission3score(plane(index));

                    plane(index) = sanityCheck(plane(index)); %make sure all the calculated values make sense and meet 
                    %competition requirements. In post-processing the only
                    %planes to be considered will be ones with a true flag
                    if plane(index).sanityFlag
                        index = index + 1; %only increment iteration when the airplane is reasonable
                    end
                    iterNum = iterNum+1;
                end
                
            end
        end
        disp("Analyzed to airfoil " + airfoil + " with preallocation")
        toc
    end
end

%Code below is from post.m file to run post-processing on the data. We should probably move this into a function call later.

score2 = zeros(1, length(plane)); %mission 2 score
score3 = score2; %mission 3 score
scoreg = score2; %ground score
vials = score2; %number of vials in given airplane
syringes = score2; %number of syringes in given airplane
for i = 1:length(plane) %load data from structure into arrays that are easier to work with
    if(plane(i).sanityFlag == 1) %all values remain zero unless the airplane passes sanity check. If the last plane analyzed is not sane the original sanity check can't throw it out.
        score2(i) = (plane(i).performance.score2);
        score3(i) = floor((plane(i).performance.score3));
        vials(i) = plane(i).fuselage.numVials;
        syringes(i) = plane(i).fuselage.numSyringes;
        scoreg(i) = 10 + 2*(3*syringes(i)/5) + 5*vials(i);%run, load and unload syringes, load vial
    end
end
[M2, I2] = max(score2); %find the airplanes with the best individual mission scores
[M3, I3] = max(score3);
[Mg, Ig] = max(scoreg);
score2 = 1 + score2/M2; %normalize scores against best performers
score3 = 2 + score3/M3;
scoreg = scoreg/Mg;
score = score2 + score3 + scoreg + 1;
[M, I] = maxk(score, 100); %find the top 100 airplanes
winners = plane(I);
scores = score(I);
save("winnersAR" + Aspect_Ratios(1) + ".mat", "winners"); %save top 100 airplanes. It is impractical to save all airplanes checked.
save("winnersAR" + Aspect_Ratios(1) + "_scores.mat", "scores");
%clear; %once finished don't keep hogging ram. Previous experience with running multiple instances on analysis on one computer shows that ram usage is the first limiting factor in enabling
       %the analysis to run, so clearing it as often as possible is important to avoid hogging ram from other programs. This ram hogging is the leading cause of crashing for this code.
%scatter(vials,score);