%TODO: implement this as a super class implemented either as major or minor
classdef PipeSection
    properties
        %TODO: implement length if want head loss major
        %length {mustBeNumeric}
        %isStrait section {must be logical}
        Kl {mustBeNumeric}
        areaInitial {mustBeNumeric}
        areaFinal {mustBeNumeric}
        
        Vinitial {mustBeNumericOrLogical} %false if unkown
        Vfinal {mustBeNumericOrLogical} %false if unkown
        
        %Hl
        Plminor
        
        deltaP
    end
    
    methods
        function obj = PipeSection(areaInitial,areaFinal,Kl,rho,Vinitial,Vfinal)
            obj.areaInitial = areaInitial;
            obj.areaFinal = areaFinal;
            obj.Kl = Kl;
            
            if Vfinal == false && Vinitial == false
                disp("either vinital or vfinal must be provided")
            elseif Vfinal == false
                obj.Vinitial = Vinitial;
                obj.Vfinal = obj.getUnknownV(true); %continuity eq
            elseif Vinitial == false
                obj.Vfinal = Vfinal;
                obj.Vinitial = obj.getUnknownV(false); %continuity eq
            end
            
            obj.Plminor = obj.getPLoss(obj.Vinitial,rho);
            obj.deltaP = obj.getDeltaP(rho);
        end
        
        function V = getUnknownV(self,initialVKnown)
            %Uses incompressible continuity equation
            if initialVKnown %we dont know VInitial
                V = (self.areaInitial*self.Vinitial) / self.areaFinal;
            else %we dont know VFinal
                V = (self.areaFinal*self.Vfinal) / self.areaInitial;
            end
        end
        
        function deltaP = getDeltaP(self,rho)
            deltaP = (1/2) * rho * (self.Vinitial^2 - self.Vfinal^2) - self.Plminor;
            
        end
        
        function Hl = getHLoss(self,V)
            Hl = getHLossMinor(self,V,g);
        end
        
        function Dp = getPLoss(self,V,rho)
            Dp = getPLossMinor(self,V,rho);
        end
    end
    
    methods (Access = private)
        function Hl = getHLossMinor(self,V,g)
            Hl = self.Kl * (V^2) / (2*g);
        end
        
        function Dp = getPLossMinor(self,V,rho)
            Dp = self.Kl * rho * (V^2) / 2;
        end
    end
end