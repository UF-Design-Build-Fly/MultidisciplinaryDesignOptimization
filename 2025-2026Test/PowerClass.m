classdef PowerClass
  properties
    motorName string = ""
    cells = -1; kv = -1; propDiameter = -1; propPitch = -1
    voltage = -1; rpm = -1; amps = -1; watts = -1; time = -1
    thrust = -1;           % [kg] from table (kgf)
    propSpeed = -1;        % [m/s]
    efficiency = -1; batteryCapacity = -1; weight = -1 % [Wh], [kg]
  end

  methods (Static)
    function ps = SetPowerSystemData(ps, motorTable, idx)
        gv = @(n,d) iff(n,motorTable.(n)(idx), d);
        names = motorTable.Properties.VariableNames;
        has  = @(n) any(strcmpi(n,names));
        ps.motorName     = string(gv('Motor',""));
        ps.cells         = gv('Cells',-1);
        ps.kv            = gv('kV',-1);
        ps.propDiameter  = gv('PropellerDiameterinches',-1);
        ps.propPitch     = gv('PropellerPitchinches',-1);
        ps.voltage       = gv('VoltageV',-1);
        ps.rpm           = gv('RPMs',-1);
        ps.amps          = gv('CurrentA',-1);
        ps.watts         = gv('PowerConsumptionW',-1);
        ps.time          = gv('WhFlightTimeSeconds',-1);
        ps.thrust        = gv('ThrustKg',-1);
        ps.efficiency    = gv('Efficiencythrustwatt100',-1);
        ps.weight        = gv('EstimatedSystemWeightKg',-1);

        if has('BatteryAvailableWh')
            ps.batteryCapacity = gv('BatteryAvailableWh',-1);
        elseif isfinite(ps.watts) && isfinite(ps.time) && ps.watts>0 && ps.time>0
            ps.batteryCapacity = ps.watts * (ps.time/3600); % Wh
        else
            ps.batteryCapacity = 100; % stay rule-compliant by default
        end

        % pitch speed
        if isfinite(ps.rpm) && isfinite(ps.propPitch) && ps.rpm>0 && ps.propPitch>0
            ps.propSpeed = (ps.rpm/60) * (ps.propPitch*0.0254);
        else
            V = max(ps.voltage, max(1,ps.cells)*3.7);
            rpm_est = max(ps.rpm, ps.kv*V);
            if isfinite(rpm_est) && rpm_est>0 && isfinite(ps.propPitch) && ps.propPitch>0
                ps.propSpeed = (rpm_est/60) * (ps.propPitch*0.0254);
            else
                ps.propSpeed = 60;
            end
        end

        function y=iff(cond,a,b), if cond, y=a; else, y=b; end; end
    end
  end
end
