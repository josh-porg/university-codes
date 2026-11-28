classdef PlyFailure
    %PlyFailure An object To contain the details of a failure in a ply
    %   keeps trak of what has failed and how
    %   FEILDS:
    %   METHODS:

    properties
        plyNumber
        failureMode % 1-> axial tensile; 2-> transverse tensile; 3-> axial compressive; 4-> transverse compressive; 5-> shear
        failureCriterion
        failureLoad
    end

    methods
        function this = PlyFailure(plyNumber, failureMode, failureCriterion, failureLoad)
            %PlyFailure Construct an instance of this class
            %   Initializes the properties with the given input arguments
            this.plyNumber = plyNumber;
            this.failureMode = failureMode;
            this.failureCriterion = failureCriterion;
            this.failureLoad = failureLoad;
        end

        function str = toString(this)
            %toString Converts the object to a string
            %   Used to print out a summary of the object

            % cheap hack to avoid enums
            failureCritStrings = ["strain", "stress"];
            failureModeStrings = ["tensile", "compressive"];
            failureDirStrings = ["axial", "transverse", "shear"];
            % convert ply code to string
            failCritStr = failureCritStrings((this.failureCriterion == "Max_Stress") +1);

            if this.failureMode == 5 % shear failure
                failureModeStr = failureDirStrings(3);
            else
                %failureModeStr = failureDirStrings(ceil(this.failureMode/2)) + " " + failureModeStrings(mod(this.failureMode-1,2)+1); % this doesnt work
                failureModeStr = failureDirStrings(mod(this.failureMode-1,2)+1) + " " + failureModeStrings(ceil(this.failureMode/2)); % this line should wordk correctly
            end

            loadstr = sprintf('%.3g, ',this.failureLoad); loadstr = loadstr(1:end-2);
            str = sprintf("Ply %g has failed in failed in max %s %s at load [%s]", this.plyNumber, failureModeStr, failCritStr, loadstr);
        end
    end
end