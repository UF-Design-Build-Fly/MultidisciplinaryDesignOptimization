Aspect_Ratios = 4:0.2:15;
syringes = 10:100;
Num_Airfoils = 1; %DEBUG
Num_Power_Systems = 658; %DEBUG - match this to spreadheet size when loaded below

plane = struct(airplaneClass);
plane.fuselage.wheelSA = 5; %DEBUG - update with actual calculated value

load("MotorSpreadsheet.mat");

%call the wings function here, then iterate through by index in the loops below

index = 1;
for AR = 1:length(Aspect_Ratios)
    for airfoilIndex = 1:Num_Airfoils
        for powerIndex = 1:Num_Power_Systems
            for syringe_index = 1:length(syringes)
                for num_vials = 1:floor(num_syringes/10)
                    %set wing properties
    
                    plane(index) = powerSelections(plane(index), MotorSpreadsheet, powerIndex);
                    plane(index).fuselage.numSyringes = syringes(syringe_index);
                    plane(index).fuselage.numVials = num_vials;
                    plane(index).fuselage = SaestsNew(plane(index).fuselage);

                    plane(index) = genVelocitySolver(plane(index), 2); %2 signifies mission 2 configuration
                    plane(index) = takeoffCheck(plane(index));
                    if ~sanityCheck(plane(index))
                        plane(index).sanityFlag = false;
                    end
                    plane(index).performanceClass = mission2score(plane(index));

                    plane(index) = genVelocitySolver(plane(index), 3); %2 signifies mission 2 configuration
                    %plane(index) = takeoffCheck(plane(index));
                    if ~sanityCheck(plane(index))
                        plane(index).sanityFlag = false;
                        break; %stop wasting time scoring this airplane if it can't complete a mission or otherwise doesn't make sense as an aircraft configuration
                    end
                    %plane(index).performanceClass = mission3score(plane(index));



                    index = index + 1; %keep track of loop iteration number
                end
            end
        end
    end
end