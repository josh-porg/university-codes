%% create MPC controller object with sample time
mpc1 = mpc(plant_C, 0.05);
%% specify prediction horizon
mpc1.PredictionHorizon = 100;
%% specify control horizon
mpc1.ControlHorizon = 20;
%% specify nominal values for inputs and outputs
mpc1.Model.Nominal.U = [0;5];
mpc1.Model.Nominal.Y = [336.8;0.0872664625997164;0;0;0];
%% specify scale factors for inputs and outputs
mpc1.MV(2).ScaleFactor = 10;
mpc1.OV(1).ScaleFactor = 1000;
mpc1.OV(5).ScaleFactor = 10000;
%% specify constraints for MV and MV Rate
mpc1.MV(1).Min = -0.349065850398866;
mpc1.MV(1).Max = 0.349065850398866;
mpc1.MV(2).Min = 0;
mpc1.MV(2).Max = 10;
%% specify constraints for OV
mpc1.OV(1).Min = 0;
mpc1.OV(1).Max = 500;
mpc1.OV(2).Min = -0.279252680319093;
mpc1.OV(2).Max = 0.279252680319093;
%% specify constraint softening for OV
mpc1.OV(1).MinECR = 0;
mpc1.OV(1).MaxECR = 0.5;
mpc1.OV(2).MinECR = 0.2;
mpc1.OV(2).MaxECR = 0.5;
%% specify overall adjustment factor applied to weights
beta = 0.13534;
%% specify weights
mpc1.Weights.MV = [0 0]*beta;
mpc1.Weights.MVRate = [0.25 0.1]/beta;
mpc1.Weights.OV = [0.75 0.01 0.25 0.01 1]*beta;
mpc1.Weights.ECR = 100000;
%% specify simulation options
options = mpcsimopt();
options.MVSignal = mpc1_MVSignal;
options.RefLookAhead = 'off';
options.MDLookAhead = 'off';
options.Constraints = 'on';
options.OpenLoop = 'off';
%% run simulation
sim(mpc1, 201, mpc1_RefSignal, mpc1_MDSignal, options);
