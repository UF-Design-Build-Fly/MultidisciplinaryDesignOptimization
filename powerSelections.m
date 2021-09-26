function power = powerSelections(cells)
%powerSelection: Given number of battery cells return ideal power system
%from a lookup table.
%There is an ideal combination of battery, motor, propeller size, etc.. for
%each given battery voltage. This function takes a battery voltage input
%and returns all the parameters necessary to construct the ideal system.
%Future work could make ideal selection more automated.
%powerSystem: systems from 3 to 8 cells.
%            Motor Name,  KV, Propeller diameter (inches), Propeller pitch (inches),
%            Voltage, RPM, Current Draw (A), Power Consumption (W),
%            Flight time (minutes), thrust (pound force), Pitch speed
%            (ft/s), Efficiency (thrust/watt)*100, battery available watt
%            hours, estimated system weight (pounds).
powerSystem = zeros(15, 15);
%                'name' kv   pd  pp    v    rpm   amps    watts time        lbf  p-speed    effi      battW miss1,       miss2             
powerSystem(1,:) = [2	760	13	6.5	 11.1	6449  23.42	   260	12.46		3.591	58.2	1.381153846	54	1.123571035	1.123571035];
powerSystem(2,:) = [3	470	17	7	14.8	5972	35.03	518.4	11.13		7.194	58.1	1.387731481	96.2	2.097488987	2.097488987];	
powerSystem(3,:) = [4	1080	9	4.5	18.5	16703	51.78	958	11.58663883		6.592	104.4	0.688100209	185	1.817439427	3.308628855];
powerSystem(4,:) = [5	300	20	8	22.2	5632	45.22	1003.9	11.94142843		12.578	62.6	1.252913637	199.8	2.6769163	4.273832599];	
powerSystem(5,:) = [6	300	18	8	25.9	6724	40.92	1059.7	9.531943003		12.172	74.7	1.14862697	168.35	2.315682819	3.551365639];
powerSystem(6,:) = [7	300	16	8	29.6	7203	35.91	1063	10.85983067		11.495	80.1	1.081373471	192.4	2.280613987	3.688102974];
power=powerSystem(cells,:);
end