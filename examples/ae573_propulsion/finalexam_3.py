"""AE 573 final exam problem 3 (``FinalExam_3``): axial compressor stage, mean-line analysis.

Mean radius 0.4 m, shaft speed 774.926 rad/s (U = 310 m/s), axial velocity
150 m/s held through the stage, no inlet swirl (c1 = c3 = 150 m/s axial),
degree of reaction 0.75, solidity 1.5, total-pressure loss coefficient
0.05; inlet Tt1 = 287 K, pt1 = 100 kPa; gamma = 1.4, R = 287.

Errors in the MATLAB relations: ``Rdeg = 1 - (ct1-ct2)/2*U`` (should be
``1 - (ct1 + ct2) / (2 U)``), ``a2 = (gamma*R*T2)`` without the square
root, ``T3 = Tt3 - c3/2*cp`` (should be ``c3^2 / (2 cp)``),
``Tt1/T1 = 1+(r-1)/2*M1^2`` used the radius ``r`` as gamma, and
``M1 = Ccz1/(gamma*R*T)^.5`` used undefined names. The loss coefficient is
applied to the rotor (relative frame) and to the stator.
"""

from unicodes.propulsion_cycles import compressor_stage

for label, loss in (("isentropic", 0.0), ("with losses (omega = 0.05)", 0.05)):
    s = compressor_stage(0.4, 774.926, 150.0, 0.75, 287.0, 100e3, 1.5, loss_rotor=loss, loss_stator=loss)
    print(label)
    print(f"  U = {s.U:.2f} m/s, ct1 = {s.ct1:.1f}, ct2 = {s.ct2:.2f} m/s, work = {s.work / 1e3:.3f} kJ/kg")
    print(f"  beta1 = {s.beta1:.2f} deg, beta2 = {s.beta2:.2f} deg, alpha2 = {s.alpha2:.2f} deg")
    print(f"  M1 = {s.M1:.4f}, M1 relative = {s.M1_relative:.4f}, M2 = {s.M2:.4f}, M3 = {s.M3:.4f}")
    print(f"  Tt3 = {s.Tt3:.2f} K, pt3 = {s.pt3:.0f} Pa, stage pressure ratio = {s.pressure_ratio:.4f}, "
          f"efficiency = {s.efficiency:.4f}")
    print(f"  diffusion factors: rotor {s.diffusion_factor_rotor:.4f}, stator {s.diffusion_factor_stator:.4f}")
