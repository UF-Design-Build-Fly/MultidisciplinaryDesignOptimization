classdef WingClass
    properties
        span = -1;            % [m]
        clw = -1;             % wing CL (cruise)
        clm = -1;             % CLmax
        cd  = -1;             % wing Cd (cruise)
        weight = -1;          % [kg]
        name = -1;
        surfaceArea = -1;     % [m^2] ~ 2*planform
        chord = -1;           % [m]
        aspectRatio = -1;     % AR
        planformArea = -1;    % [m^2]
        clFlap = -1;          % CLmax with flaps
        thickness = -1;       % [m]
    end

    methods (Static)
        function wing = SetWingData(wing, wings, airfoilIndex, aspectRatioIndex, spanIndex)
            sz = size(wings);
            airfoilIndex     = min(max(airfoilIndex,1), sz(1));
            aspectRatioIndex = min(max(aspectRatioIndex,1), sz(3));
            spanIndex        = min(max(spanIndex,1), sz(4));

            wing.clw          = wings(airfoilIndex, 1, aspectRatioIndex, spanIndex);
            wing.clm          = wings(airfoilIndex, 2, aspectRatioIndex, spanIndex);
            wing.cd           = max(0, wings(airfoilIndex, 3, aspectRatioIndex, spanIndex));
            wing.clFlap       = max(wing.clm, wings(airfoilIndex, 4, aspectRatioIndex, spanIndex));
            wing.weight       = max(0, wings(airfoilIndex, 5, aspectRatioIndex, spanIndex));
            wing.chord        = max(1e-6, wings(airfoilIndex, 6, aspectRatioIndex, spanIndex));
            wing.planformArea = max(1e-6, wings(airfoilIndex, 7, aspectRatioIndex, spanIndex));
            wing.surfaceArea  = max(1e-6, wings(airfoilIndex, 8, aspectRatioIndex, spanIndex));
            wing.name         = wings(airfoilIndex, 9, aspectRatioIndex, spanIndex);
            wing.thickness    = max(0, wings(airfoilIndex,10, aspectRatioIndex, spanIndex));
        end

        function maxLiftN = FindMaxLift(wingObj, airSpeed, rho)
            CL = wingObj.clFlap; if ~(isfinite(CL) && CL>0), CL = wingObj.clm; end
            S  = max(1e-6, wingObj.planformArea);
            maxLiftN = 0.5 * rho * (airSpeed.^2) * S * max(0, CL);
        end
    end
end
