Aspect_Ratios = 4:0.2:15;
Num_Airfoils = 1; %constants use Underscore_And_Caps. Variables use camelCase
Num_Power_Systems = 1;

plane = struct(airplaneClass);
index = 0;

%call the wings function here, then iterate through by index in the loops below


for AR = 1:length(Aspect_Ratios)
    for airfoilIndex = 1:Num_Airfoils
        for powerIndex = 1:Num_Power_Systems
            for batterySize = 10:10:100 %watt hours of given battery
                for num_syringes = 10:50
                    for num_vials = 1:floor(num_syringes/10)
                        index = index + 1;
                        
                        %set wing properties
        
                        plane(index) = powerSelections(powerIndex); %TODO: Update power selections to work with plane struct
                        plane(index).fuselage.numSyringes = num_syringes;
                        plane(index).fuselage.numVials = num_vials;
                        plane(index).fuselage = SaestsNew(plane(index).fuselage);

                        %plane(index) = genVelocitySolver(plane(index), 2); %2 signifies mission 2 configuration
                        %plane(index) = takeoffCheck(plane(index));
                        if ~sanityCheck(plane(index))
                            plane(index).sanityFlag = false;
                            break; %stop wasting time scoring this airplane if it can't complete a mission or otherwise doesn't make sense as an aircraft configuration
                        end
                        %plane(index).performanceClass = mission2score(plane(index));

                        %plane(index) = genVelocitySolver(plane(index), 3); %2 signifies mission 2 configuration
                        %plane(index) = takeoffCheck(plane(index));
                        if ~sanityCheck(plane(index))
                            plane(index).sanityFlag = false;
                            break; %stop wasting time scoring this airplane if it can't complete a mission or otherwise doesn't make sense as an aircraft configuration
                        end
                        %plane(index).performanceClass = mission3score(plane(index));
                    
                        
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
end