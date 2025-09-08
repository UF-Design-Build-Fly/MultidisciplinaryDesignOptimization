function plane = FindTotalWeight(plane)
    % Defensive read of banner mechanism mass (default 0.20 kg if property missing)
    if isprop(plane.fuselage,'bannerMechanismMass') && ~isempty(plane.fuselage.bannerMechanismMass) ...
            && isnumeric(plane.fuselage.bannerMechanismMass) && isfinite(plane.fuselage.bannerMechanismMass)
        bannerMech = max(0, plane.fuselage.bannerMechanismMass);
    else
        bannerMech = 0.20;  % default mass if class doesn't define the property yet
    end

    plane.performance.totalEmptyWeight = ...
          max(0, plane.empennage.HSweight) ...
        + max(0, plane.empennage.VSweight) ...
        + max(0, plane.fuselage.weight) ...
        + max(0, plane.fuselage.gearWeight) ...
        + max(0, plane.fuselage.wheelWeight) ...
        + max(0, plane.powerSystem.weight) ...
        + max(0, plane.wing.weight) ...
        + bannerMech;

    % M2 payload (ducks/pucks)
    mp = plane.performance.massPerPassenger;
    mc = plane.performance.massPerCargo;
    p  = plane.performance.passengers;
    c  = plane.performance.cargo;
    plane.performance.m2Weight = max(0, p*mp + c*mc);

    x1 = max(0, plane.performance.x1_addon_kg);

    plane.performance.totalWeight2 = plane.performance.totalEmptyWeight + plane.performance.m2Weight + x1;
    plane.performance.totalWeight3 = plane.performance.totalEmptyWeight + x1; % banner hardware is in empty
end
