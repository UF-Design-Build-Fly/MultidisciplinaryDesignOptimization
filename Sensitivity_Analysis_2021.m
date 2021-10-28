clear; clc;
warning('off','all')
%rho=0.00235308; %Desnity air at Tuscon, Az (slug/ft^3)
rho = 0.002391; %Density air at Whichita, Ks with average climate data from April 2021
Aspect_Ratios = 8:.5:15;%the wing aspect ratios being considered 
syringes = 10:100;
load("MotorSpreadsheet.mat");
Num_Power_Systems = height(MotorSpreadsheet);

width_wheel = 0.5;    %width of wheels (in)
radius_wheel = 1.5; %wheel radius
VSAR = 2;
HSAR = 4.5;

[wings] = wingData(Aspect_Ratios, 8); %call wing function to make airfoil data lookup table
[wingrow, wingcol, wingpg] = size(wings);

index = 1;
max_index = wingpg*wingrow*Num_Power_Systems*length(syringes)*length(syringes);
for AR = 1:wingpg
    for airfoil = 1:wingrow
        for powerIndex = 1:Num_Power_Systems
            for syringe_index = 1:length(syringes)
                for num_vials = 1:floor(syringes(syringe_index)/10)
                    plane(index) = struct(airplaneClass);
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
                       
                    index = index + 1; %keep track of loop iteration number
                end
            end
        end
    end
end