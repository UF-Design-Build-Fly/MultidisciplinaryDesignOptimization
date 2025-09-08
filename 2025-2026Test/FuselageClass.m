classdef FuselageClass
    % All SI units unless noted
    properties
        frontalSurfaceArea = -1   % [m^2]
        height = -1               % [m]
        width  = -1               % [m]
        length = -1               % [m]
        totalSA = -1              % [m^2] wetted area

        weight = -1               % [kg] fuselage shell mass

        gearWeight = -1           % [kg]
        gearSA = -1               % [m^2] planform of gear strips
        gearFrontalSA = -1        % [m^2]

        wheelWidth = -1           % [m]
        wheelRadius = -1          % [m]
        wheelWeight = -1          % [kg]
        wheelSA = -1              % [m^2] wetted area
        wheelFrontalSA = -1       % [m^2]

        % Included in empty weight for all missions (same configuration rule)
        bannerMechanismMass = 0.20 % [kg] spool/tether/servo fairing etc.
    end

    methods (Static)
        function fuselage = CalcFuselageData(plane)
            % Geometry & mass estimate for the fuselage shell
            % Uses simple primitives (nose frustum + rectangular bay + boom)
            in2m = 0.0254;

            % Basic layout (adjust if you change bay dimensions)
            noseLength = (7 + 4) * in2m;     % battery + motor bay
            fusH = 4 * in2m;                 % fuselage height
            fusW = 4 * in2m;                 % fuselage width
            bayLen = 6 * in2m;               % payload compartment length

            % Fuselage length proportional to span, with floor
            totalLen = max(0.40, 0.90 * max(plane.wing.span, 0.1));
            tailLen  = max(0, totalLen - noseLength - bayLen);

            % Nose lateral area (frustum-like)
            l_slant = hypot(noseLength, (1.5*in2m - 2*in2m));
            SA_nose = pi * (1.5*in2m + 2*in2m) * l_slant;

            % Payload bay area (top + two sides + diagonal face)
            topA  = bayLen * fusW;
            sideA = 2 * 0.5 * bayLen * fusH;
            diagA = hypot(bayLen, fusH) * fusW;
            SA_bay = topA + sideA + diagA;

            % Tail boom (cylindrical shell)
            boomD = 1 * in2m;
            SA_tail = pi * boomD * tailLen;

            SA_total = SA_nose + SA_bay + SA_tail;

            % Skin thickness & material
            t_skin = (1/16) * in2m;   % 1/16 in carbon laminate (effective)
            rhoCF  = 1281;            % kg/m^3 (typical CF composite)

            % Start from existing object to preserve bannerMechanismMass if set
            fuselage = plane.fuselage;

            % If bannerMechanismMass is invalid/missing, apply default
            if ~(isprop(fuselage,'bannerMechanismMass') && isnumeric(fuselage.bannerMechanismMass) ...
                 && isfinite(fuselage.bannerMechanismMass) && fuselage.bannerMechanismMass >= 0)
                fuselage.bannerMechanismMass = 0.20; % kg
            end

            % Set computed fields
            fuselage.frontalSurfaceArea = fusW * fusH;
            fuselage.height = fusH;
            fuselage.width  = fusW;
            fuselage.length = totalLen;
            fuselage.totalSA = SA_total;

            % Shell mass only (mechanisms handled separately)
            fuselage.weight = SA_total * t_skin * rhoCF;
        end

        function fuselage = GenLandingGear(plane)
            % Simple triangular strut + two wheels estimate
            in2m = 0.0254;

            numWheels   = 2;
            wheelWidth  = 0.9063 * in2m;
            wheelRadius = 1.5 * in2m;

            gearStripWidth = 1.5 * in2m;   % strip width (spanwise)
            gearThickness  = (1/8) * in2m; % laminate thickness (frontal thickness)
            base           = 1.5 * in2m;   % vertical straight section
            rhoCF          = 1281;         % kg/m^3

            % Ground clearance target: prop radius - fuselage half-height + extra
            widthSpan = 0.25 * max(plane.wing.span, 0.1);
            heightReq = plane.powerSystem.propDiameter * in2m / 2 - 2*in2m + 3*in2m;

            fuselage = plane.fuselage;
            fusW = max(fuselage.width, 0.05);  % guard

            a = 0.5 * (widthSpan - fusW);
            b = max(0, heightReq - base);
            c = hypot(a, b);                   % diagonal leg

            % Total centerline length of strut path
            L = fusW + 2*(c + base);

            % Areas & masses
            fuselage.gearSA        = L * gearStripWidth;         % wetted (planar) area
            fuselage.gearFrontalSA = L * gearThickness;          % frontal area for Cd*A
            fuselage.gearWeight    = rhoCF * (fuselage.gearSA * gearThickness);

            fuselage.wheelWeight    = numWheels * 0.0241;        % ~0.85 oz each → 0.0241 kg
            fuselage.wheelSA        = numWheels * (2*pi*wheelRadius^2 + 2*pi*wheelRadius*wheelWidth);
            fuselage.wheelFrontalSA = numWheels * (2*wheelRadius*wheelWidth);

            fuselage.wheelWidth  = wheelWidth;
            fuselage.wheelRadius = wheelRadius;
        end
    end
end
