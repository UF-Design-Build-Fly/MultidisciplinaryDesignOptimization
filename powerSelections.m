function plane = powerSelections(plane, table, row)
%Given the motor spreadsheet and index, function will return power
%configuration for corresponding index/row. Thus function essentially
%converts spreadsheet entries into parameters for the "plane" object.
%Motor Spreadsheet must be loaded before hand at the start of the analysis
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

    plane.power.motorName = table(row,:).motorName;
    plane.power.cells = table(row,:).cells;
    plane.power.kv = table(row,:).kv;
    plane.power.propDiameter = table(row,:).propDiameter;
    plane.power.propPitch = table(row,:).propPitch;
    plane.power.voltage = table(row,:).voltage;
    plane.power.rpm = table(row,:).rpm;
    plane.power.amps = table(row,:).current;
    plane.power.watts = table(row,:).watts;
    plane.power.time = table(row,:).time;
    plane.power.thrust = table(row,:).thrust;
    plane.power.propSpeed = table(row,:).propSpeed;
    plane.power.efficiency = table(row,:).efficiency;
    plane.power.batteryCapacity = table(row,:).batteryCapacity;
    plane.power.weight = table(row,:).estimatedWeight/16;
    
end
