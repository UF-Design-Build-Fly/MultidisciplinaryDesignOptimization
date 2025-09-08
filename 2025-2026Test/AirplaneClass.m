classdef AirplaneClass
    properties
        wing = struct(WingClass)
        empennage = struct(EmpennageClass)
        powerSystem = struct(PowerClass)
        fuselage = struct(FuselageClass)
        performance = struct(PerformanceClass)
        failureReason string = "Not checked"
        sanityFlag = 1
        takeoffFail = 0; momentFail = 0; convergeFail = 0
    end
end
