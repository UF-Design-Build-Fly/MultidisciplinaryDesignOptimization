function [newPlane] = powerSelections(plane, table, row)
%Given the motor spreadsheet and index, function will return power
%configuration for corresponding index/row.
%Motor Spreadsheet MUST be loaded before hand at the start of the analysis
%and passed through the function each time it is called.

%There is an ideal combination of battery, motor, propeller size, etc.. for
%each given battery voltage. This function returns all the parameters 
%necessary to construct the ideal system.
%Future work could make ideal selection more automated.
%powerSystem: systems from 3 to 6 cells.
%            Motor Name,  KV, Propeller diameter (inches), Propeller pitch (inches),
%            Voltage, RPM, Current Draw (A), Power Consumption (W),
%            Flight time (minutes), thrust (pound force), Pitch speed
%            (ft/s), Efficiency (thrust/watt)*100, battery available watt
%            hours, estimated system weight (pounds).

plane.powerSystem.motorName = table(row,:).motorName;
plane.powerSystem.cells = table(row,:).cells;
plane.powerSystem.kv = table(row,:).kv;
plane.powerSystem.propDiameter = table(row,:).propDiameter;
plane.powerSystem.propPitch = table(row,:).propPitch;
plane.powerSystem.voltage = table(row,:).voltage;
plane.powerSystem.rpm = table(row,:).rpm;
plane.powerSystem.amps = table(row,:).current;
plane.powerSystem.watts = table(row,:).watts;
plane.powerSystem.time = table(row,:).time;
plane.powerSystem.thrust = table(row,:).thrust;
plane.powerSystem.propSpeed = table(row,:).propSpeed;
plane.powerSystem.efficiency = table(row,:).efficiency;
plane.powerSystem.batteryCapacity = table(row,:).batteryCapacity;
plane.powerSystem.weight = table(row,:).estimatedWeight;

newPlane = plane;
end
