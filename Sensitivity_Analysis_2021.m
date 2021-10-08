Num_Airfoils = 1; %constants use Underscore_And_Caps. Variables use camelCase
Aspect_Ratios = 1:2;%4:0.2:15;
Num_Power_Systems = 1;

plane = airplaneClass;
plane.wingClass.cl = 5;
plane.empennageClass.HSarea = 5;
plane.powerSystemClass.cells = 6;
plane.fuselageClass.weight = 10;
plane.performanceClass.weight2 = 5;

for AR = 1:length(Aspect_Ratios)
    for airfoilIndex = 1:Num_Airfoils
        for powerIndex = 1:Num_Power_Systems
            for num_syringes = 10:10:50
                for num_containers = 1:floor(num_syringes/10)
                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   %Jack Code
                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   %Build aircraft configuration %(new call to airplane constructor)
                        %wing
                        %empennage
                        %powerSystem
                        %fuselage
                    
                    %calculate m2 performance with estimated syringe capacity
                        %genVelocitySolver
                        %liftCalc/Check
                        %takeoffCalc/Check
                    %calculate m2 score

                    %calculate m3 performance with estimated vial capacity
                        %genVelocitySolver
                        %liftCalc/Check
                        %takeoffCalc/Check
                    %calculate m3 score


                    %Scoring and ranking:

                    %Option A: easy and fastest
                        %keep just the best performing airplane with a global best variable
                        %use global increment to track index

                    %Option B: easy and slow
                        %Keep vector of all aircraft tested in position by loop order (first tested at top)

                    %Option C: harder but faster
                        %Keep top N aircraft in ranked score order

                    %either way there ought to be a good method to export the results to an excel file for sorting/analysis later
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
                   



                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Outline of old analysis
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %estimate fuselage surface area
                    %calculate aircraft weight
                    %run velocity/drag solver
                    %make sure there's enough lift for the weight
                    %run takeoff distance checker
                    %display solver progress for debugging
                    %calculate m2 score
                    %see Jack's paper copy for more complete outline, or better yet just read the old code ;)
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
            end
        end
    end
end