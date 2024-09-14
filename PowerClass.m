classdef PowerClass

  properties

	motorName = ""; 
	cells = -1;
	kv = -1;
	propDiameter = -1;
	propPitch = -1;
	voltage = -1;
	rpm = -1;
	amps = -1;
	watts = -1;
	time = -1;
	thrust = -1;
	propSpeed = -1;
	efficiency = -1;
	batteryCapacity = -1; %watt hours
	weight = -1;

	%propDrag1 = -1;
	%progDrag2 = -1;
	%propDrag3 = -1;
	
	%dThrust1 = -1;
	%dThrust2 = -1;
	%dThrust3 = -1;
	
  end

  methods (Static)

		function powerSystem = SetPowerSystemData(powerSystem, motorTable, tableIndex)
			
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
			%			Motor Name,  KV, Propeller diameter (inches), Propeller pitch (inches),
			%			Voltage, RPM, Current Draw (A), Power Consumption (W),
			%			Flight time (minutes), thrust (pound force), Pitch speed
			%			(ft/s), Efficiency (thrust/watt)*100, battery available watt
			%			hours, estimated system weight (pounds).
			
			powerSystem.motorName = char(motorTable(tableIndex,:).Motor);
			powerSystem.cells = motorTable(tableIndex,:).Cells;
			powerSystem.kv = motorTable(tableIndex,:).kV;
			powerSystem.propDiameter = motorTable(tableIndex,:).PropellerDiameterinches;
			powerSystem.propPitch = motorTable(tableIndex,:).PropellerPitchinches;
			powerSystem.voltage = motorTable(tableIndex,:).VoltageV;
			powerSystem.rpm = motorTable(tableIndex,:).RPMs;
			powerSystem.amps = motorTable(tableIndex,:).CurrentA;
			powerSystem.watts = motorTable(tableIndex,:).PowerConsumptionW;
			powerSystem.time = motorTable(tableIndex,:).WhFlightTimeSeconds;
			powerSystem.thrust = motorTable(tableIndex,:).Thrustlbs;
			powerSystem.propSpeed = motorTable(tableIndex,:).PitchSpeedfts*1.467; %Tabel values are in mph not ft/s
			powerSystem.efficiency = motorTable(tableIndex,:).Efficiencythrustwatt100;
			powerSystem.weight = motorTable(tableIndex,:).EstimatedSystemWeightlbs;

		end

  end

end