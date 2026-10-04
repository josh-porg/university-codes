# MATLAB to Python conversion map

Every `.m` and `.mlx` file under `matlab_source/` and where its Python port lives, or why it was not
ported. Library code is under `src/unicodes/`; scripts (homework, labs, exams, projects) are under
`examples/<course>/`. Generated from the docstrings of the Python files plus a table of rules, so a file
listed against several Python files is mentioned in each.

Totals: 937 ported or folded into a port, 31 not ported (empty, plotting-only, drafts, unfinished), 209 excluded missile work.

Corrections made to the originals are listed in the docstring of each Python function or script.

## AE 345

| MATLAB | Python / note |
|---|---|
| `AE 345/PipeSection.m` | `examples/ae345_fluids/homework.py` |
| `AE 345/Wind_Tunnel_Ex8_5.mlx` | `examples/ae345_fluids/homework.py` |
| `AE 345/homework_8.mlx` | `examples/ae345_fluids/homework.py` |
| `AE 345/wind_tunnel_driver.mlx` | `examples/ae345_fluids/homework.py` |

## AE 445

| MATLAB | Python / note |
|---|---|
| `AE 445/HW_11.mlx` | Not ported: only reads the problem 3.1 spreadsheets (used by examples/ae445_aerodynamics/hw12_critical_mach.py) |
| `AE 445/HW_12.mlx` | examples/ae445_aerodynamics/hw12_critical_mach.py |
| `AE 445/HW_14_and_15_problem_4_4.mlx` | `examples/ae445_aerodynamics/hw14_drag_polar.py` |
| `AE 445/Homework_4.mlx` | examples/ae445_aerodynamics/hw04_free_vortex.py |
| `AE 445/Homework_7.mlx` | examples/ae445_aerodynamics/hw07_flat_plate.py |
| `AE 445/Test_2.mlx` | examples/ae445_aerodynamics/exams.py |
| `AE 445/Test_3.mlx` | examples/ae445_aerodynamics/exams.py |
| `AE 445/deltafromThetaTropo.mlx` | src/unicodes/atmosphere.py |
| `AE 445/findClosestVectorIntersection.mlx` | Replaced by root finding (`brentq`) in examples/ae445_aerodynamics/hw12_critical_mach.py |
| `AE 445/getAlphaEffecive.mlx` | src/unicodes/aero/wing.py (`induced_angle`; alpha_eff = alpha - alpha_i) |
| `AE 445/getCLAlphaWIngFromAirfoil.mlx` | `src/unicodes/aero/wing.py` |
| `AE 445/getCLalpha.mlx` | `src/unicodes/aero/wing.py` |
| `AE 445/getCLalphafromClAlpha.mlx` | `src/unicodes/aero/wing.py` |
| `AE 445/getCLaplhaRectangularIncompresable.mlx` | `src/unicodes/aero/wing.py` |
| `AE 445/getCfLaminar.mlx` | `src/unicodes/aero/boundary_layer.py` |
| `AE 445/getCfTurbulentIncompressable.mlx` | `src/unicodes/aero/boundary_layer.py` |
| `AE 445/getCpCrit.mlx` | `src/unicodes/gasdynamics.py` |
| `AE 445/getCriticalMachNumber2D.mlx` | `src/unicodes/aero/wing.py` |
| `AE 445/getCriticalMachNumber3D.mlx` | `src/unicodes/aero/wing.py` |
| `AE 445/getEffectiveAspectRatio.mlx` | `src/unicodes/aero/wing.py` |
| `AE 445/getInducedAngleOfAttack.mlx` | `src/unicodes/aero/wing.py` |
| `AE 445/getInducedDragCoefficient.mlx` | `src/unicodes/aero/wing.py` |
| `AE 445/getLaminarFlatPlateDrag.mlx` | src/unicodes/aero/boundary_layer.py |
| `AE 445/getLaminarTurbulentBoundryLayerDistance.mlx` | `src/unicodes/aero/boundary_layer.py` |
| `AE 445/getMachNumberEffective.mlx` | `src/unicodes/aero/wing.py` |
| `AE 445/getMuAir.mlx` | `src/unicodes/atmosphere.py` |
| `AE 445/getParabolicDragPolar.mlx` | `src/unicodes/aero/wing.py` |
| `AE 445/getSatndardAtmosphericValues.mlx` | src/unicodes/atmosphere.py |
| `AE 445/getSatndardAtmosphericValuesStrato.mlx` | src/unicodes/atmosphere.py |
| `AE 445/getSatndardAtmosphericValuesTropo.mlx` | src/unicodes/atmosphere.py |
| `AE 445/getSpdSoundAir.mlx` | `src/unicodes/atmosphere.py` |
| `AE 445/getStallSpeed.mlx` | `src/unicodes/aero/wing.py` |
| `AE 445/getStandardAtmoRatiosStrato.mlx` | src/unicodes/atmosphere.py |
| `AE 445/getStandardAtmoRatiosTropo.mlx` | src/unicodes/atmosphere.py |
| `AE 445/getStandardAtmosphericValues.mlx` | `src/unicodes/atmosphere.py` |
| `AE 445/getStandardAtmosphericValuesBG.mlx` | src/unicodes/atmosphere.py |
| `AE 445/getTempRatioTropo.mlx` | src/unicodes/atmosphere.py |
| `AE 445/getTurbulentIncompressableFlatPlateDrag.mlx` | src/unicodes/aero/boundary_layer.py |
| `AE 445/ias2tas.mlx` | `src/unicodes/atmosphere.py` |
| `AE 445/karmanTsienTransformation.mlx` | src/unicodes/aero/wing.py |
| `AE 445/laitoneTransformation.mlx` | src/unicodes/aero/wing.py |
| `AE 445/pitostaicPressures2velocity.mlx` | src/unicodes/gasdynamics.py (`mach_from_pitot_*`) and examples/ae445_aerodynamics/exams.py |
| `AE 445/polhamus.mlx` | `src/unicodes/aero/wing.py` |
| `AE 445/prandtlGlauertTransformation.mlx` | src/unicodes/aero/wing.py |
| `AE 445/quiz_4.mlx` | examples/ae445_aerodynamics/exams.py |
| `AE 445/sigmafromThetaTropo.mlx` | src/unicodes/atmosphere.py |
| `AE 445/standardAtmosphericTableGenerator.mlx` | `examples/ae445_aerodynamics/standard_atmosphere_table.py`, `src/unicodes/atmosphere.py` |

## AE 507

| MATLAB | Python / note |
|---|---|
| `AE 507/Design_project_thickness_solver.mlx` | `examples/ae507_structures/homework.py` |
| `AE 507/HW_1.mlx` | `examples/ae507_structures/homework.py` |
| `AE 507/HW_2_2_4.mlx` | `examples/ae507_structures/homework.py` |
| `AE 507/HW_2_5.mlx` | `examples/ae507_structures/homework.py` |
| `AE 507/in_class_2_9_2022.mlx` | `examples/ae507_structures/homework.py` |

## AE 546

| MATLAB | Python / note |
|---|---|
| `AE 546/CpFromHead.mlx` | `src/unicodes/gasdynamics.py` |
| `AE 546/CpFromPressure.mlx` | `src/unicodes/gasdynamics.py` |
| `AE 546/HW_2.mlx` | examples/ae546_aero_lab/hw02_venturi.py |
| `AE 546/HW_3.mlx` | examples/ae546_aero_lab/hw03_lift.py |
| `AE 546/HW_4.mlx` | examples/ae546_aero_lab/hw04_thin_airfoil.py |
| `AE 546/Lab_1.mlx` | `examples/ae546_aero_lab/lab1_pressure_distribution.py` |
| `AE 546/Lab_1_altered_data.mlx` | examples/ae546_aero_lab/lab1_pressure_distribution.py |
| `AE 546/Lab_1_preliminary.m` | `examples/ae546_aero_lab/lab1_pressure_distribution.py` |
| `AE 546/Lab_1_version_two_this_time_with_fewer_functions_hopefully.mlx` | examples/ae546_aero_lab/lab1_pressure_distribution.py |
| `AE 546/Lab_2_driver.mlx` | `examples/ae546_aero_lab/lab2_supersonic_wedge.py` |
| `AE 546/Lab_2_function.mlx` | `examples/ae546_aero_lab/lab2_supersonic_wedge.py` |
| `AE 546/Lab_3.mlx` | examples/ae546_aero_lab/lab3_boundary_layer.py |
| `AE 546/Lab_3_equations.mlx` | examples/ae546_aero_lab/lab3_boundary_layer.py |
| `AE 546/MachFromPitoCompresible.mlx` | `src/unicodes/gasdynamics.py` |
| `AE 546/MachFromPitoSupersonic.mlx` | `src/unicodes/gasdynamics.py` |
| `AE 546/debuggingExpansionFan.mlx` | src/unicodes/gasdynamics.py (`expansion_fan`) |
| `AE 546/expansionFanProperties.mlx` | `src/unicodes/gasdynamics.py` |
| `AE 546/honours_project_xfoil_analysis.mlx` | Not ported: plots a saved XFOIL polar / camber line only (see src/unicodes/aero/xfoil.py, aero/thin_airfoil.py) |
| `AE 546/honours_thin_airfoil.mlx` | Not ported: plots a saved XFOIL polar / camber line only (see src/unicodes/aero/xfoil.py, aero/thin_airfoil.py) |
| `AE 546/isentropicFlowProperties.mlx` | `src/unicodes/gasdynamics.py` |
| `AE 546/lab1_ploter_maybe.m` | `examples/ae546_aero_lab/lab1_pressure_distribution.py` |
| `AE 546/laminarBoundryLayerThickness.mlx` | src/unicodes/aero/boundary_layer.py |
| `AE 546/latexEquations.mlx` | Not ported: LaTeX equation formatting only |
| `AE 546/machWaveAngle.mlx` | `src/unicodes/gasdynamics.py` |
| `AE 546/naca5gen.m` | `src/unicodes/aero/naca.py` |
| `AE 546/normalShockProperties.mlx` | `src/unicodes/gasdynamics.py` |
| `AE 546/obliqueShockProperties.mlx` | `src/unicodes/gasdynamics.py` |
| `AE 546/polhamus.mlx` | src/unicodes/aero/wing.py (`polhamus_lift_slope`) |
| `AE 546/prandtlMeyerFunction.mlx` | `src/unicodes/gasdynamics.py` |
| `AE 546/saveAllFig.mlx` | Not ported: saves open MATLAB figures (use `--save` options / `savefig`) |
| `AE 546/table2latex.m` | `src/unicodes/io/latex.py` |
| `AE 546/tst_naca5gen.m` | src/unicodes/aero/naca.py (tested in tests/test_aero_basics.py) |
| `AE 546/tubulentBoundryLayerThickness.mlx` | src/unicodes/aero/boundary_layer.py |

## AE 550

| MATLAB | Python / note |
|---|---|
| `AE 550/HIB.mlx` | `src/unicodes/flight_dynamics/kinematics.py` |
| `AE 550/HW_2.mlx` | `examples/ae550_flight_dynamics/homework.py` |
| `AE 550/HW_3.mlx` | `examples/ae550_flight_dynamics/homework.py` |
| `AE 550/HW_4.mlx` | Symbolic linearisation; see examples/ae550_flight_dynamics/homework.py docstring |
| `AE 550/HW_4_part_2.mlx` | Symbolic linearisation; see examples/ae550_flight_dynamics/homework.py docstring |
| `AE 550/aircraft_aerodynamic_center_derivation.mlx` | `src/unicodes/flight_dynamics/kinematics.py` |
| `AE 550/axisAngle2qaternion.mlx` | `src/unicodes/flight_dynamics/kinematics.py` |
| `AE 550/downwashGradient.mlx` | `src/unicodes/aero/wing.py` |
| `AE 550/eulerAngles2quaternion.mlx` | `src/unicodes/flight_dynamics/kinematics.py` |
| `AE 550/final_exam.mlx` | `examples/ae550_flight_dynamics/final_exam.py` |
| `AE 550/final_exam_part_2.mlx` | `examples/ae550_flight_dynamics/final_exam_part_2.py` |
| `AE 550/getAirFlowAngles.mlx` | `src/unicodes/flight_dynamics/kinematics.py` |
| `AE 550/getFlightPathAngles.mlx` | `src/unicodes/flight_dynamics/kinematics.py` |
| `AE 550/polhamus.mlx` | `src/unicodes/aero/wing.py` |
| `AE 550/quatConj.mlx` | src/unicodes/flight_dynamics/kinematics.py |
| `AE 550/quatInv.mlx` | src/unicodes/flight_dynamics/kinematics.py |
| `AE 550/quatMult.mlx` | `src/unicodes/flight_dynamics/kinematics.py` |
| `AE 550/quatRot.mlx` | `src/unicodes/flight_dynamics/kinematics.py` |
| `AE 550/quatRotPassive.mlx` | `src/unicodes/flight_dynamics/kinematics.py` |
| `AE 550/quaternion2eulerAngles.mlx` | `src/unicodes/flight_dynamics/kinematics.py` |
| `AE 550/rotateVectorWithQuaternion.mlx` | `src/unicodes/flight_dynamics/kinematics.py` |

## AE 571

| MATLAB | Python / note |
|---|---|
| `AE 571/final project/AE571_FinalProject_Team8.m` | `examples/ae571_combustion/final_project.py` |
| `AE 571/final project/other teams codes/AE571_FinalProject_Group9.m` | `examples/ae571_combustion/other_teams/group9.py` |
| `AE 571/final project/other teams codes/AE_571_Group7_Final_Code.m` | `examples/ae571_combustion/other_teams/group7.py` |
| `AE 571/final project/other teams codes/Group3_AE571_Final_Project.m` | `examples/ae571_combustion/other_teams/group3.py` |
| `AE 571/final project/other teams codes/OneDrive_1_12-16-2022/AE571_FinalProject_PartI.m` | `examples/ae571_combustion/other_teams/onedrive_team.py` |
| `AE 571/final project/other teams codes/OneDrive_1_12-16-2022/AE571_FinalProject_PartII.m` | examples/ae571_combustion/other_teams/onedrive_team.py |
| `AE 571/final project/other teams codes/OneDrive_1_12-16-2022/AE571_functions/atombalancefinal.m` | `examples/ae571_combustion/other_teams/onedrive_team.py` |
| `AE 571/final project/other teams codes/OneDrive_1_12-16-2022/AE571_functions/enthalpycalculator2.m` | `examples/ae571_combustion/other_teams/onedrive_team.py` |
| `AE 571/final project/other teams codes/OneDrive_1_12-16-2022/AE571_functions/enthalpycalculator2H2.m` | examples/ae571_combustion/other_teams/onedrive_team.py |
| `AE 571/final project/other teams codes/OneDrive_1_12-16-2022/AE571_functions/enthalpycalculatorfinal.m` | `examples/ae571_combustion/other_teams/onedrive_team.py` |
| `AE 571/final project/other teams codes/OneDrive_1_12-16-2022/AE571_functions/enthalpycalculatorfinalH2.m` | examples/ae571_combustion/other_teams/onedrive_team.py |
| `AE 571/final project/other teams codes/OneDrive_1_12-16-2022/AE571_functions/specificheats.m` | `examples/ae571_combustion/other_teams/onedrive_team.py` |
| `AE 571/final project/other teams codes/OneDrive_2022-12-16/Team 2/AE571_Group2_Final_Code.m` | `examples/ae571_combustion/other_teams/team2.py` |
| `AE 571/final project/other teams codes/team 10/Part_I.m` | `examples/ae571_combustion/other_teams/team10.py` |
| `AE 571/final project/other teams codes/team 10/Part_II.m` | `examples/ae571_combustion/other_teams/team10.py` |
| `AE 571/final project/other teams codes/team 4/Team 4/AE_571_Final_Project_Code.m` | `examples/ae571_combustion/other_teams/team4.py` |
| `AE 571/final project/other teams codes/team 5/Corrected Code that works/AE571_FinalProject_Part1.m` | `examples/ae571_combustion/other_teams/team5.py` |
| `AE 571/final project/other teams codes/team 5/Corrected Code that works/AE571_Part1_Calculators.m` | `examples/ae571_combustion/other_teams/team5.py` |
| `AE 571/final project/other teams codes/team 5/Corrected Code that works/Part 2/AE571_FinalProject_Part2_Complete.m` | `examples/ae571_combustion/other_teams/team5.py` |
| `AE 571/final project/other teams codes/team 5/Corrected Code that works/Part 2/AtomBalanceLab1.m` | `examples/ae571_combustion/other_teams/team5.py` |
| `AE 571/final project/other teams codes/team 5/Corrected Code that works/Part 2/enthalpycalculatorAir.m` | `examples/ae571_combustion/other_teams/team5.py` |
| `AE 571/final project/other teams codes/team 5/Corrected Code that works/Part 2/enthalpycalculatorFuel.m` | examples/ae571_combustion/other_teams/team5.py |
| `AE 571/final project/other teams codes/team 5/Corrected Code that works/Part 2/enthalpycalculatorO2.m` | examples/ae571_combustion/other_teams/team5.py |
| `AE 571/final project/other teams codes/team 5/Corrected Code that works/Part 2/enthalpycalculatorProducts.m` | examples/ae571_combustion/other_teams/team5.py |
| `AE 571/final project/other teams codes/team 6/Team 6/Gasoline.m` | `examples/ae571_combustion/other_teams/team6.py` |
| `AE 571/final project/other teams codes/team 6/Team 6/Hydrogen.m` | `examples/ae571_combustion/other_teams/team6.py` |
| `AE 571/final project/other teams codes/team 6/Team 6/atombalance.m` | `examples/ae571_combustion/other_teams/team6.py` |
| `AE 571/final project/other teams codes/team 6/Team 6/propertycalculator.m` | `examples/ae571_combustion/other_teams/team6.py` |
| `AE 571/final project/poznanski_main.mlx` | `examples/ae571_combustion/final_project.py` |
| `AE 571/final project/team_code.mlx` | `examples/ae571_combustion/final_project.py` |
| `AE 571/getAdiabFlameTemp.mlx` | `examples/ae571_combustion/lab1_flame_temperature.py` |
| `AE 571/getHFromTables.mlx` | `examples/ae571_combustion/lab1_flame_temperature.py` |
| `AE 571/getTable1.mlx` | `examples/ae571_combustion/lab1_flame_temperature.py` |
| `AE 571/getTable2.mlx` | `examples/ae571_combustion/lab1_flame_temperature.py` |
| `AE 571/lab1_hFromTable.mlx` | `examples/ae571_combustion/lab1_flame_temperature.py` |
| `AE 571/lab_1_Read_Data_Table_1.mlx` | Excel table readers; replaced by NASA polynomials in src/unicodes/thermo/combustion.py |
| `AE 571/lab_1_Read_Data_Table_2.mlx` | Excel table readers; replaced by NASA polynomials in src/unicodes/thermo/combustion.py |
| `AE 571/lab_1_poznanski.mlx` | `examples/ae571_combustion/lab1_flame_temperature.py` |
| `AE 571/lab_2/AE_571_lab_2_Extra_credit_Poznanski.mlx` | `examples/ae571_combustion/lab2_compression_work.py` |
| `AE 571/lab_2/AE_571_lab_2_Poznanski.mlx` | `examples/ae571_combustion/lab2_compression_work.py` |
| `AE 571/savePathTest.mlx` | Not ported: file-path experiment |
| `AE 571/str2Elements.mlx` | `examples/ae571_combustion/lab1_flame_temperature.py`, `src/unicodes/thermo/combustion.py` |
| `AE 571/string_chem_formula.mlx` | `examples/ae571_combustion/lab1_flame_temperature.py` |

## Autumn 2024/AE 709 Composites

| MATLAB | Python / note |
|---|---|
| `Autumn 2024/AE 709 Composites/Design_Project.mlx` | `examples/ae709_composites/design_project.py` |
| `Autumn 2024/AE 709 Composites/Exam_1_Problem_4.mlx` | `examples/ae709_composites/exam_1.py` |
| `Autumn 2024/AE 709 Composites/Exam_1_Problem_5.mlx` | `examples/ae709_composites/exam_1.py` |
| `Autumn 2024/AE 709 Composites/Exam_1_Problem_7.mlx` | `examples/ae709_composites/exam_1.py` |
| `Autumn 2024/AE 709 Composites/HIB.mlx` | `examples/ae709_composites/design_project.py` |
| `Autumn 2024/AE 709 Composites/HW_4.mlx` | `examples/ae709_composites/homework.py` |
| `Autumn 2024/AE 709 Composites/HW_5a.mlx` | `examples/ae709_composites/homework.py` |
| `Autumn 2024/AE 709 Composites/HW_5b.mlx` | `examples/ae709_composites/homework.py` |
| `Autumn 2024/AE 709 Composites/Laminate.m` | src/unicodes/composites/ |
| `Autumn 2024/AE 709 Composites/Ply.m` | `examples/ae709_composites/exam_1.py` |
| `Autumn 2024/AE 709 Composites/PlyFailure.m` | src/unicodes/composites/ |

## Autumn 2024/AE 725 Numerical Optimization

| MATLAB | Python / note |
|---|---|
| `Autumn 2024/AE 725 Numerical Optimization/BuckStresses.mlx` | `src/unicodes/structures/wingbox.py` |
| `Autumn 2024/AE 725 Numerical Optimization/ConstrainedSteepestDescentTest.mlx` | `examples/ae725_optimization/global_optimization.py` |
| `Autumn 2024/AE 725 Numerical Optimization/DescentFunction.mlx` | `src/unicodes/optimize/descent.py` |
| `Autumn 2024/AE 725 Numerical Optimization/FINAL_CODE_ALL_ASSEMBLED_PART_I_12_12_24_AT_11_PM.mlx` | examples/ae725_optimization/wingbox_optimization.py |
| `Autumn 2024/AE 725 Numerical Optimization/Final_Exam_Problem_2.mlx` | examples/ae725_optimization/final_exam.py |
| `Autumn 2024/AE 725 Numerical Optimization/Final_Exam_Problem_4.mlx` | `examples/ae725_optimization/final_exam.py` |
| `Autumn 2024/AE 725 Numerical Optimization/Final_Exam_Problem_4_FEA.mlx` | `examples/ae725_optimization/final_exam.py` |
| `Autumn 2024/AE 725 Numerical Optimization/Final_Exam_Problem_5.mlx` | examples/ae725_optimization/final_exam.py |
| `Autumn 2024/AE 725 Numerical Optimization/Final_Exam_Problem_5_to_print.mlx` | examples/ae725_optimization/final_exam.py |
| `Autumn 2024/AE 725 Numerical Optimization/GeneticAlgorithm.mlx` | `src/unicodes/optimize/genetic.py` |
| `Autumn 2024/AE 725 Numerical Optimization/GlobalOptimization_CSD_Test.mlx` | `examples/ae725_optimization/global_optimization.py` |
| `Autumn 2024/AE 725 Numerical Optimization/GlobalOptimization_ContinuousSimulatedAnealing_Test.mlx` | `examples/ae725_optimization/global_optimization.py` |
| `Autumn 2024/AE 725 Numerical Optimization/Homework_10_Problem_13_1.mlx` | examples/ae725_optimization/homework.py |
| `Autumn 2024/AE 725 Numerical Optimization/Homework_1_Problem_3_51.mlx` | `examples/ae725_optimization/homework.py` |
| `Autumn 2024/AE 725 Numerical Optimization/Homework_3_Problem_4_43.mlx` | examples/ae725_optimization/homework.py |
| `Autumn 2024/AE 725 Numerical Optimization/Homework_3_Problem_4_59.mlx` | examples/ae725_optimization/homework.py |
| `Autumn 2024/AE 725 Numerical Optimization/Homework_5_Problem_5_4.mlx` | examples/ae725_optimization/homework.py |
| `Autumn 2024/AE 725 Numerical Optimization/Homework_5_Problem_5_62.mlx` | examples/ae725_optimization/homework.py |
| `Autumn 2024/AE 725 Numerical Optimization/Homework_5_Problem_7_10.mlx` | examples/ae725_optimization/homework.py |
| `Autumn 2024/AE 725 Numerical Optimization/Homework_5_Problem_7_5.mlx` | examples/ae725_optimization/homework.py |
| `Autumn 2024/AE 725 Numerical Optimization/Homework_5_Problem_SP1.mlx` | `examples/ae725_optimization/homework.py` |
| `Autumn 2024/AE 725 Numerical Optimization/Homework_7_Problem_10_52.mlx` | examples/ae725_optimization/homework.py |
| `Autumn 2024/AE 725 Numerical Optimization/Homework_7_Problem_10_57.mlx` | examples/ae725_optimization/homework.py |
| `Autumn 2024/AE 725 Numerical Optimization/MasterCodeLiveScript (1).mlx` | examples/ae725_optimization/wingbox_optimization.py |
| `Autumn 2024/AE 725 Numerical Optimization/ParticleSwarmOptimization.mlx` | `examples/ae725_optimization/wingbox_optimization.py` |
| `Autumn 2024/AE 725 Numerical Optimization/ParticleSwarmOptimization_2.mlx` | `examples/ae725_optimization/wingbox_optimization.py` |
| `Autumn 2024/AE 725 Numerical Optimization/Results_Summary_Latex_Converter.mlx` | `examples/ae725_optimization/wingbox_optimization.py` |
| `Autumn 2024/AE 725 Numerical Optimization/Stresses.mlx` | `src/unicodes/structures/wingbox.py` |
| `Autumn 2024/AE 725 Numerical Optimization/Test_goldenSearch.mlx` | `examples/ae725_optimization/global_optimization.py` |
| `Autumn 2024/AE 725 Numerical Optimization/Wingbox_SA_V_0.mlx` | `examples/ae725_optimization/wingbox_optimization.py`, `src/unicodes/optimize/descent.py` |
| `Autumn 2024/AE 725 Numerical Optimization/bucklebox_constraints.mlx` | `src/unicodes/structures/wingbox.py` |
| `Autumn 2024/AE 725 Numerical Optimization/computeJacobian.mlx` | `src/unicodes/optimize/descent.py` |
| `Autumn 2024/AE 725 Numerical Optimization/constrained_steepest_descent.mlx` | `src/unicodes/optimize/descent.py` |
| `Autumn 2024/AE 725 Numerical Optimization/continuous_simulated_anealing.mlx` | `src/unicodes/optimize/descent.py` |
| `Autumn 2024/AE 725 Numerical Optimization/debugJacobian.mlx` | src/unicodes/numerics.py (`jacobian`) |
| `Autumn 2024/AE 725 Numerical Optimization/genetic_algorithm_optimization.mlx` | `src/unicodes/optimize/genetic.py` |
| `Autumn 2024/AE 725 Numerical Optimization/golden_search.mlx` | `src/unicodes/optimize/descent.py` |
| `Autumn 2024/AE 725 Numerical Optimization/inexactStepSize.mlx` | `src/unicodes/optimize/descent.py` |
| `Autumn 2024/AE 725 Numerical Optimization/particle_swarm_optimization.mlx` | `src/unicodes/optimize/swarm.py` |
| `Autumn 2024/AE 725 Numerical Optimization/results/bucklebox_constraints.mlx` | `src/unicodes/structures/wingbox.py` |
| `Autumn 2024/AE 725 Numerical Optimization/string_Stresses.mlx` | `src/unicodes/structures/wingbox.py` |
| `Autumn 2024/AE 725 Numerical Optimization/stringbox_constraints.mlx` | `src/unicodes/structures/wingbox.py` |
| `Autumn 2024/AE 725 Numerical Optimization/stringy_boi.mlx` | `src/unicodes/structures/wingbox.py` |
| `Autumn 2024/AE 725 Numerical Optimization/wingbox_constraints.mlx` | `src/unicodes/structures/wingbox.py` |
| `Autumn 2024/AE 725 Numerical Optimization/wingbox_objective_function.mlx` | `src/unicodes/structures/wingbox.py` |
| `Autumn 2024/AE 725 Numerical Optimization/wingbox_optimization_CSD_to_submit.mlx` | `examples/ae725_optimization/wingbox_optimization.py` |
| `Autumn 2024/AE 725 Numerical Optimization/wingbox_optimization_CSD_to_submit_extra_cred_working.mlx` | examples/ae725_optimization/wingbox_optimization.py |
| `Autumn 2024/AE 725 Numerical Optimization/wingbox_optimization_CSD_to_submit_working_minimimums_1_2.mlx` | examples/ae725_optimization/wingbox_optimization.py |
| `Autumn 2024/AE 725 Numerical Optimization/wingbox_optimization_GA.mlx` | `examples/ae725_optimization/wingbox_optimization.py` |
| `Autumn 2024/AE 725 Numerical Optimization/wingbox_optimization_continuousSA.mlx` | `examples/ae725_optimization/wingbox_optimization.py` |

## Autumn 2024/Missile Design Continued

| MATLAB | Python / note |
|---|---|
| `Autumn 2024/Missile Design Continued/BendingFrequency_Body_first.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/HIB.mlx` | `src/unicodes/flight_dynamics/kinematics.py` |
| `Autumn 2024/Missile Design Continued/Integration_of_boost_eom.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/L_over_D_for_qbars.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/L_over_D_for_qbars_V2-ENGR-HWC66X3.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/L_over_D_for_qbars_V2.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/L_over_D_for_qbars_V3.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/L_over_D_for_qbars_V4.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/L_over_D_for_qbars_V4_fleeman_reproducer.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/L_over_D_for_qbars_V5_altitudes.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/Planar_Surface_Drag_Verification.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/Plot_Wave_drag_benchmark_Round.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/Propulsion/2D_quasisentropic_inlet_geometry.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/Propulsion/4bar_duct_sheet_cuts_CAD.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/Propulsion/FaringGeometry_Triangularization.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/Propulsion/Hermite_spline.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/Propulsion/InitialSizing.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/Propulsion/IsentropicSpike_3D_from_table.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/Propulsion/Propulsion_Function_Fleeman_Verification.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/Propulsion/QuasiIsentropic_Inlet_2D.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/Propulsion/Scale_Based_Sizing.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/Propulsion/THermodynamic_sizing.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/Propulsion/circle_interception.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/Propulsion/develop_oblique_conical_frustum.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/Propulsion/flame_propagation_v0.m` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/Propulsion/importAxisymmetricSpike.m` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/Propulsion/loop_fairing_geometery_derivation.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/Propulsion/loop_fairing_geometery_triangularization.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/Propulsion/pipe_sizing_derivation.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/Propulsion/quasisentropic_inlet_geometry_3D.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/Propulsion/ramjet_area_ratio_inletThroat_to_combustor.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/Propulsion/ramjet_area_ratio_inletThroat_to_freeStreemFlow.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/Propulsion/ramjet_mach_combustorEntrance_thermalChoking.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/Propulsion/ramjet_temperature_combustorExit.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/Propulsion/sliding_duct_linkage_angle.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/Propulsion/sliding_duct_sizing.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/RocketBullet.m` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/Spline_For_transonic_body_drag.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/aerodynamic_Center_Body.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/aerodynamic_Center_Body_Alone.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/aerodynamic_center_planar_surface.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/ballistic_trajectory_simulator.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/ballistic_trajectory_simulator_V2.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/ballistic_trajectory_simulator_V2_for_joesRound.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/ballistic_trajectory_simulator_V3.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/ballistic_trajectory_simulator_V3_ASAT.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/ballistic_trajectory_simulator_V3_RAIDER_Bass.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/ballistic_trajectory_simulator_V3_joes_round.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/ballistic_trajectory_simulator_range_angle_computer.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/componentBuildup_aerodynamicCenter.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/dragCoeff_friction_surface.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/dragCoeff_wave_surface.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/dragCoeff_wave_zeroLift_bluntedNose.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/dragCoeff_zeroLift_body.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/dragCoeff_zeroLift_body_smoothed_transonic.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/dragCoeff_zeroLift_surface.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/drag_body.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/file to send to bass group/RocketBullet.m` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/file to send to bass group/ballistic_trajectory_simulator_V3_RAIDER_Bass.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/file to send to bass group/dragCoeff_wave_zeroLift_bluntedNose.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/file to send to bass group/dragCoeff_zeroLift_body.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/file to send to bass group/dragCoeff_zeroLift_body_smoothed_transonic.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/file to send to bass group/drag_body.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/file to send to bass group/save_bullet_info.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/flare_aerodynamic_center.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/getAirFlowAngles.mlx` | `src/unicodes/flight_dynamics/kinematics.py` |
| `Autumn 2024/Missile Design Continued/getFlightPathAngles.mlx` | `src/unicodes/flight_dynamics/kinematics.py` |
| `Autumn 2024/Missile Design Continued/hinge_moment.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/lift_to_drag_ratio_body.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/normalForceCoeff_liftingBody.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/normalForceCoeff_planar_surface.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/normalForce_Flare_withRespectTo_alpha.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/normalForce_withRespectTo_alpha.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/normalForce_withRespectTo_alpha_planar_surface.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/plot_C_D_0_surface_vs_mach.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/plot_C_N_vs_alpha_for_lifting_body.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/plot_aerodynamicCeneteroverNoseLength_v_s_alpha.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/plot_lift_vs_drag_body_example.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/plot_tail_area_sizing.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/plot_x_AC_planar_surface.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/polt_C_D_0_body_vs_mach.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/save_bullet_info.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/simulate_balistic_Trajectory.mlx` | Not ported: missile design (excluded as requested) |
| `Autumn 2024/Missile Design Continued/tail_area_sizing.mlx` | Not ported: missile design (excluded as requested) |

## Autumn 2025/Math 796

| MATLAB | Python / note |
|---|---|
| `Autumn 2025/Math 796/Basic_NN/activate.m` | `src/unicodes/neural.py` |
| `Autumn 2025/Math 796/Basic_NN/inexactStepSize.mlx` | `src/unicodes/optimize/descent.py` |
| `Autumn 2025/Math 796/Basic_NN/netbp.m` | `examples/math796_optimization/basic_nn.py`, `src/unicodes/neural.py` |
| `Autumn 2025/Math 796/Basic_NN/netbpfull.m` | `examples/math796_optimization/basic_nn.py`, `src/unicodes/neural.py` |
| `Autumn 2025/Math 796/Basic_NN/netbpfull_modified_other_optimizers.m` | `examples/math796_optimization/basic_nn.py` |
| `Autumn 2025/Math 796/Basic_NN/netbpfull_modified_step.m` | `examples/math796_optimization/basic_nn.py` |
| `Autumn 2025/Math 796/Basic_NN/nlsrun.m` | `examples/math796_optimization/basic_nn.py`, `src/unicodes/neural.py` |
| `Autumn 2025/Math 796/HW_1/HW1_script.m` | `examples/math796_optimization/hw1_optimizer_comparison.py` |
| `Autumn 2025/Math 796/HW_1/Homework_1.mlx` | `examples/math796_optimization/hw1_optimizer_comparison.py` |
| `Autumn 2025/Math 796/HW_1/optimization_test_functions/bowl_shaped/boha1.m` | src/unicodes/optimize/test_functions.py |
| `Autumn 2025/Math 796/HW_1/optimization_test_functions/bowl_shaped/perm0db.m` | src/unicodes/optimize/test_functions.py |
| `Autumn 2025/Math 796/HW_1/optimization_test_functions/bowl_shaped/rothyp.m` | src/unicodes/optimize/test_functions.py |
| `Autumn 2025/Math 796/HW_1/optimization_test_functions/bowl_shaped/sumsqu.m` | src/unicodes/optimize/test_functions.py |
| `Autumn 2025/Math 796/HW_1/optimization_test_functions/bowl_shaped/trid.m` | src/unicodes/optimize/test_functions.py |
| `Autumn 2025/Math 796/HW_1/optimization_test_functions/many_local_minima/ackley.m` | src/unicodes/optimize/test_functions.py |
| `Autumn 2025/Math 796/HW_1/optimization_test_functions/many_local_minima/drop.m` | src/unicodes/optimize/test_functions.py |
| `Autumn 2025/Math 796/HW_1/optimization_test_functions/many_local_minima/griewank.m` | src/unicodes/optimize/test_functions.py |
| `Autumn 2025/Math 796/HW_1/optimization_test_functions/many_local_minima/langer.m` | src/unicodes/optimize/test_functions.py |
| `Autumn 2025/Math 796/HW_1/optimization_test_functions/many_local_minima/levy.m` | src/unicodes/optimize/test_functions.py |
| `Autumn 2025/Math 796/HW_1/optimization_test_functions/many_local_minima/levy13.m` | src/unicodes/optimize/test_functions.py |
| `Autumn 2025/Math 796/HW_1/optimization_test_functions/many_local_minima/rastr.m` | src/unicodes/optimize/test_functions.py |
| `Autumn 2025/Math 796/HW_1/optimization_test_functions/many_local_minima/schaffer2.m` | src/unicodes/optimize/test_functions.py |
| `Autumn 2025/Math 796/HW_1/optimization_test_functions/many_local_minima/shubert.m` | src/unicodes/optimize/test_functions.py |
| `Autumn 2025/Math 796/HW_1/optimization_test_functions/other/beale.m` | src/unicodes/optimize/test_functions.py |
| `Autumn 2025/Math 796/HW_1/optimization_test_functions/other/branin.m` | src/unicodes/optimize/test_functions.py |
| `Autumn 2025/Math 796/HW_1/optimization_test_functions/other/braninmodif.m` | src/unicodes/optimize/test_functions.py |
| `Autumn 2025/Math 796/HW_1/optimization_test_functions/other/goldpr.m` | src/unicodes/optimize/test_functions.py |
| `Autumn 2025/Math 796/HW_1/optimization_test_functions/other/permdb.m` | src/unicodes/optimize/test_functions.py |
| `Autumn 2025/Math 796/HW_1/optimization_test_functions/other/stybtang.m` | src/unicodes/optimize/test_functions.py |
| `Autumn 2025/Math 796/HW_1/optimization_test_functions/plate_shaped/booth.m` | src/unicodes/optimize/test_functions.py |
| `Autumn 2025/Math 796/HW_1/optimization_test_functions/plate_shaped/matya.m` | src/unicodes/optimize/test_functions.py |
| `Autumn 2025/Math 796/HW_1/optimization_test_functions/plate_shaped/mccorm.m` | src/unicodes/optimize/test_functions.py |
| `Autumn 2025/Math 796/HW_1/optimization_test_functions/plate_shaped/powersum.m` | src/unicodes/optimize/test_functions.py |
| `Autumn 2025/Math 796/HW_1/optimization_test_functions/plate_shaped/zakharov.m` | src/unicodes/optimize/test_functions.py |
| `Autumn 2025/Math 796/HW_1/optimization_test_functions/steep_ridges_and_drops/dejong5.m` | src/unicodes/optimize/test_functions.py |
| `Autumn 2025/Math 796/HW_1/optimization_test_functions/steep_ridges_and_drops/easom.m` | src/unicodes/optimize/test_functions.py |
| `Autumn 2025/Math 796/HW_1/optimization_test_functions/steep_ridges_and_drops/michal.m` | src/unicodes/optimize/test_functions.py |
| `Autumn 2025/Math 796/HW_1/optimization_test_functions/valley_shaped/camel3.m` | src/unicodes/optimize/test_functions.py |
| `Autumn 2025/Math 796/HW_1/optimization_test_functions/valley_shaped/camel6.m` | src/unicodes/optimize/test_functions.py |
| `Autumn 2025/Math 796/HW_1/optimization_test_functions/valley_shaped/dixonpr.m` | src/unicodes/optimize/test_functions.py |
| `Autumn 2025/Math 796/HW_1/optimization_test_functions/valley_shaped/rosen.m` | src/unicodes/optimize/test_functions.py |
| `Autumn 2025/Math 796/HW_1/table2latex.m` | `src/unicodes/io/latex.py` |

## DMD

| MATLAB | Python / note |
|---|---|
| `DMD/CODE/CH01_INTRO/Algorithm_1_2.m` | `examples/dmd_book/ch01_intro.py` |
| `DMD/CODE/CH01_INTRO/Algorithm_1_3.m` | examples/dmd_book/ch01_intro.py |
| `DMD/CODE/CH01_INTRO/Algorithm_1_4.m` | examples/dmd_book/ch01_intro.py |
| `DMD/CODE/CH01_INTRO/Algorithm_1_5.m` | `examples/dmd_book/ch01_intro.py` |
| `DMD/CODE/CH01_INTRO/DMD.m` | `examples/data_driven/combustor_dmd.py`, `examples/data_driven/gluhareff_dmd.py`, `src/unicodes/decomposition/dmd.py` |
| `DMD/CODE/CH01_INTRO/DMDfull.m` | `examples/dmd_book/ch01_intro.py` |
| `DMD/CODE/CH02_FLUIDS/computeDMD.m` | `examples/dmd_book/ch02_fluids.py` |
| `DMD/CODE/CH02_FLUIDS/computePOD.m` | `examples/dmd_book/ch02_fluids.py` |
| `DMD/CODE/CH02_FLUIDS/loadDATA.m` | `examples/dmd_book/ch02_fluids.py` |
| `DMD/CODE/CH02_FLUIDS/loadIBPM.m` | `examples/dmd_book/ch02_fluids.py` |
| `DMD/CODE/CH02_FLUIDS/plotCylinder.m` | `examples/dmd_book/cylinder_plot.py` |
| `DMD/CODE/CH03_KOOPMAN/Algorithm_3_1.m` | `examples/dmd_book/ch03_koopman.py` |
| `DMD/CODE/CH04_VIDEO/Algorithm_4_1.m` | `examples/dmd_book/ch04_video.py` |
| `DMD/CODE/CH04_VIDEO/Algorithm_4_2.m` | examples/dmd_book/ch04_video.py |
| `DMD/CODE/CH04_VIDEO/Algorithm_4_3.m` | examples/dmd_book/ch04_video.py |
| `DMD/CODE/CH04_VIDEO/Algorithm_4_4.m` | `examples/dmd_book/ch04_video.py` |
| `DMD/CODE/CH05_MULTIRESOLUTION/mrDMD.m` | `examples/dmd_book/ch05_multiresolution.py`, `src/unicodes/decomposition/variants.py` |
| `DMD/CODE/CH05_MULTIRESOLUTION/mrDMD_demo.m` | `examples/dmd_book/ch05_multiresolution.py` |
| `DMD/CODE/CH05_MULTIRESOLUTION/mrDMD_map.m` | `examples/dmd_book/ch05_multiresolution.py`, `src/unicodes/decomposition/variants.py` |
| `DMD/CODE/CH06_DMDC/Algorithm_6_1.m` | `examples/dmd_book/ch06_dmdc.py` |
| `DMD/CODE/CH06_DMDC/Algorithm_Sec_6_2.m` | `examples/dmd_book/ch06_dmdc.py` |
| `DMD/CODE/CH07_TIMEDELAY/DMD_standingwave.m` | `examples/dmd_book/ch07_time_delay.py` |
| `DMD/CODE/CH07_TIMEDELAY/ERA.m` | `examples/dmd_book/ch07_time_delay.py`, `src/unicodes/decomposition/variants.py` |
| `DMD/CODE/CH07_TIMEDELAY/ERA_test01.m` | `examples/dmd_book/ch07_time_delay.py` |
| `DMD/CODE/CH07_TIMEDELAY/HMM_DMD.m` | `examples/dmd_book/ch07_time_delay.py` |
| `DMD/CODE/CH08_NOISEPOWER/DMD_eig.m` | `examples/dmd_book/ch08_noise_power.py`, `src/unicodes/decomposition/variants.py` |
| `DMD/CODE/CH08_NOISEPOWER/FFTDMD_spectrum.m` | `examples/dmd_book/ch08_noise_power.py` |
| `DMD/CODE/CH08_NOISEPOWER/LDS_DMD_eig.m` | `examples/dmd_book/ch08_noise_power.py` |
| `DMD/CODE/CH08_NOISEPOWER/SVHT_cylinder.m` | `examples/dmd_book/ch08_noise_power.py` |
| `DMD/CODE/CH08_NOISEPOWER/optimal_SVHT_coef.m` | `examples/dmd_book/ch08_noise_power.py` |
| `DMD/CODE/CH09_SPARSITY/Algorithm_9_1.m` | `examples/dmd_book/ch09_sparsity.py` |
| `DMD/CODE/CH09_SPARSITY/EX1_TORUS/compressedDMD.m` | `examples/dmd_book/ch09_sparsity.py` |
| `DMD/CODE/CH09_SPARSITY/EX1_TORUS/computeDMDModes.m` | `examples/dmd_book/ch09_sparsity.py` |
| `DMD/CODE/CH09_SPARSITY/EX1_TORUS/computeFFTModes.m` | `examples/dmd_book/ch09_sparsity.py` |
| `DMD/CODE/CH09_SPARSITY/EX1_TORUS/computePODModes.m` | `examples/dmd_book/ch09_sparsity.py` |
| `DMD/CODE/CH09_SPARSITY/EX1_TORUS/getParms.m` | `examples/dmd_book/ch09_sparsity.py` |
| `DMD/CODE/CH09_SPARSITY/EX1_TORUS/getSparseData.m` | `examples/dmd_book/ch09_sparsity.py` |
| `DMD/CODE/CH09_SPARSITY/EX1_TORUS/plotData.m` | `examples/dmd_book/ch09_sparsity.py` |
| `DMD/CODE/CH09_SPARSITY/EX1_TORUS/projectData.m` | `examples/dmd_book/ch09_sparsity.py` |
| `DMD/CODE/CH09_SPARSITY/EX1_TORUS/runExample.m` | `examples/dmd_book/ch09_sparsity.py` |
| `DMD/CODE/CH09_SPARSITY/EX2_CYLINDER/cDMD_p40_r21_VORT.m` | `examples/dmd_book/ch09_sparsity.py` |
| `DMD/CODE/CH09_SPARSITY/EX2_CYLINDER/csDMD_p1000_r21_VORT.m` | `examples/dmd_book/ch09_sparsity.py` |
| `DMD/CODE/CH09_SPARSITY/EX2_CYLINDER/plotCylinderNoSave.m` | `examples/dmd_book/ch09_sparsity.py`, `examples/dmd_book/cylinder_plot.py` |
| `DMD/CODE/CH09_SPARSITY/utils/cosamp.m` | `examples/dmd_book/ch09_sparsity.py`, `src/unicodes/decomposition/linalg.py` |
| `DMD/CODE/CH09_SPARSITY/utils/freezeColors.m` | `examples/dmd_book/ch09_sparsity.py` |
| `DMD/CODE/CH10_NONLINEAROBSV/Algorithm_10_1.m` | `examples/dmd_book/ch10_nonlinear_observables.py` |
| `DMD/CODE/CH10_NONLINEAROBSV/Algorithm_10_3.m` | examples/dmd_book/ch10_nonlinear_observables.py |
| `DMD/CODE/CH10_NONLINEAROBSV/Algorithm_10_4.m` | examples/dmd_book/ch10_nonlinear_observables.py |
| `DMD/CODE/CH10_NONLINEAROBSV/Algorithm_10_5.m` | examples/dmd_book/ch10_nonlinear_observables.py |
| `DMD/CODE/CH10_NONLINEAROBSV/Algorithm_10_6.m` | `examples/dmd_book/ch10_nonlinear_observables.py` |
| `DMD/CODE/CH10_NONLINEAROBSV/Algorithm_10_7through11.m` | `examples/dmd_book/ch10_nonlinear_observables.py`, `src/unicodes/decomposition/variants.py` |
| `DMD/CODE/CH10_NONLINEAROBSV/dmd_soliton_rhs.m` | `examples/dmd_book/ch10_nonlinear_observables.py` |
| `DMD/CODE/CH12_NEUROSCIENCE/DMD_ECoG.m` | `examples/dmd_book/ch12_neuroscience.py` |
| `DMD/DMD.mlx` | `examples/data_driven/combustor_dmd.py`, `examples/data_driven/gluhareff_dmd.py`, `src/unicodes/decomposition/dmd.py` |
| `DMD/DMD_combustor_driver.mlx` | `examples/data_driven/combustor_dmd.py` |
| `DMD/DMD_declarative_driver.mlx` | examples/data_driven/combustor_dmd.py and src/unicodes/decomposition/dmd.py |
| `DMD/DMD_driver.mlx` | `examples/data_driven/combustor_dmd.py`, `examples/data_driven/gluhareff_dmd.py` |
| `DMD/DMD_function_driver.mlx` | examples/data_driven/combustor_dmd.py and src/unicodes/decomposition/dmd.py |
| `DMD/DMD_function_driver_P.mlx` | examples/data_driven/combustor_dmd.py and src/unicodes/decomposition/dmd.py |
| `DMD/DMD_function_driver_T.mlx` | examples/data_driven/combustor_dmd.py and src/unicodes/decomposition/dmd.py |
| `DMD/DMD_function_driver_U.mlx` | examples/data_driven/combustor_dmd.py and src/unicodes/decomposition/dmd.py |
| `DMD/DMD_function_driver_V.mlx` | examples/data_driven/combustor_dmd.py and src/unicodes/decomposition/dmd.py |
| `DMD/DMD_function_driver_Y_CH4.mlx` | examples/data_driven/combustor_dmd.py and src/unicodes/decomposition/dmd.py |
| `DMD/DMD_function_driver_Y_CO2.mlx` | examples/data_driven/combustor_dmd.py and src/unicodes/decomposition/dmd.py |
| `DMD/DMD_function_driver_Y_H2O.mlx` | examples/data_driven/combustor_dmd.py and src/unicodes/decomposition/dmd.py |
| `DMD/DMD_function_driver_Y_O2.mlx` | examples/data_driven/combustor_dmd.py and src/unicodes/decomposition/dmd.py |
| `DMD/DMDc/DMDc.mlx` | `examples/data_driven/dmdc_flight.py`, `src/unicodes/decomposition/dmdc.py` |
| `DMD/Exporter.m` | examples/data_driven/combustor_dmd.py and src/unicodes/decomposition/dmd.py |
| `DMD/Gluhareff_Data/Data_26_2_2025/DMD.mlx` | `examples/data_driven/combustor_dmd.py`, `examples/data_driven/gluhareff_dmd.py`, `src/unicodes/decomposition/dmd.py` |
| `DMD/Gluhareff_Data/Data_26_2_2025/DMD_driver.mlx` | `examples/data_driven/combustor_dmd.py`, `examples/data_driven/gluhareff_dmd.py` |
| `DMD/Gluhareff_Data/Data_26_2_2025/Gluhareff_26_2_2025_rawData/data_extractor_V_0_Backup.mlx` | Not ported: empty backup of data_extracter (see examples/data_driven/gluhareff_dmd.py) |
| `DMD/Gluhareff_Data/Data_26_2_2025/data_extracter.mlx` | `examples/data_driven/gluhareff_dmd.py` |
| `DMD/Gluhareff_Data/Data_26_2_2025/data_extractor_V_0_Backup.mlx` | Not ported: empty backup of data_extracter (see examples/data_driven/gluhareff_dmd.py) |
| `DMD/Gluhareff_Data/Data_26_2_2025/optimal_SVHT_coef.m` | `examples/dmd_book/ch08_noise_power.py` |
| `DMD/Gluhareff_Data/Data_26_2_2025/perform_DMD_script.mlx` | `examples/data_driven/gluhareff_dmd.py` |
| `DMD/Gluhareff_Data/Data_26_2_2025/reconstruct.m` | src/unicodes/decomposition/dmd.py (`DMD.reconstruct`) |
| `DMD/InOutFormat.m` | src/unicodes/io/tecplot.py |
| `DMD/OutputTecASCIIdata.m` | `src/unicodes/io/tecplot.py` |
| `DMD/PSD_2025.mlx` | `examples/data_driven/dmd_spectrum.py`, `src/unicodes/decomposition/spectral.py` |
| `DMD/SINDy/Example Code/polynomial_generation_experiment.mlx` | Superseded by `src/unicodes/decomposition/sindy.py` (`library`) |
| `DMD/SINDy/Example Code/sparsedynamics/EX01a_Linear2D.m` | `examples/sindy_examples/ex01_linear_cubic.py` |
| `DMD/SINDy/Example Code/sparsedynamics/EX01b_Cubic2D.m` | `examples/sindy_examples/ex01_linear_cubic.py` |
| `DMD/SINDy/Example Code/sparsedynamics/EX01c_Linear3D.m` | `examples/sindy_examples/ex01_linear_cubic.py` |
| `DMD/SINDy/Example Code/sparsedynamics/EX02_Lorenz.m` | `examples/sindy_examples/ex02_lorenz.py` |
| `DMD/SINDy/Example Code/sparsedynamics/EX02_LorenzTVDiff.m` | `examples/sindy_examples/ex02_lorenz.py` |
| `DMD/SINDy/Example Code/sparsedynamics/EX03_Cylinder.m` | `examples/sindy_examples/ex03_cylinder.py` |
| `DMD/SINDy/Example Code/sparsedynamics/EX04a_LogisticMap.m` | `examples/sindy_examples/ex04_logistic_hopf.py` |
| `DMD/SINDy/Example Code/sparsedynamics/EX04b_Hopf_TVRegDiff.m` | `examples/sindy_examples/ex04_logistic_hopf.py` |
| `DMD/SINDy/Example Code/sparsedynamics/EX05_LorenzTimeDelay.m` | `examples/sindy_examples/ex05_lorenz_time_delay.py` |
| `DMD/SINDy/Example Code/sparsedynamics/EXappA_Sine.m` | `examples/sindy_examples/exappA_sine.py` |
| `DMD/SINDy/Example Code/sparsedynamics/utils/InterpretResults.mlx` | `src/unicodes/decomposition/sindy.py` |
| `DMD/SINDy/Example Code/sparsedynamics/utils/TVRegDiff.m` | `examples/data_driven/sindy_examples.py`, `examples/sindy_examples/ex02_lorenz.py`, `src/unicodes/numerics.py` |
| `DMD/SINDy/Example Code/sparsedynamics/utils/color_line3.m` | `examples/sindy_examples/common.py` |
| `DMD/SINDy/Example Code/sparsedynamics/utils/generateLibraryList.mlx` | `src/unicodes/decomposition/sindy.py` |
| `DMD/SINDy/Example Code/sparsedynamics/utils/hopf.m` | `examples/sindy_examples/ex04_logistic_hopf.py` |
| `DMD/SINDy/Example Code/sparsedynamics/utils/logistic.m` | `examples/sindy_examples/ex04_logistic_hopf.py` |
| `DMD/SINDy/Example Code/sparsedynamics/utils/lorenz.m` | `examples/sindy_examples/ex02_lorenz.py`, `src/unicodes/decomposition/sindy.py` |
| `DMD/SINDy/Example Code/sparsedynamics/utils/poolData.m` | `examples/sindy_examples/common.py` |
| `DMD/SINDy/Example Code/sparsedynamics/utils/poolDataLIST.m` | `examples/sindy_examples/common.py` |
| `DMD/SINDy/Example Code/sparsedynamics/utils/sparseGalerkin.m` | `examples/sindy_examples/common.py` |
| `DMD/SINDy/Example Code/sparsedynamics/utils/sparseGalerkinDiscrete.m` | examples/sindy_examples/ex04_logistic_hopf.py (the identified logistic map is iterated directly) |
| `DMD/SINDy/Example Code/sparsedynamics/utils/sparsifyDynamics.m` | `examples/sindy_examples/common.py` |
| `DMD/SINDy/InterpretResults.mlx` | `src/unicodes/decomposition/sindy.py` |
| `DMD/SINDy/SINDyGalerkin.mlx` | `src/unicodes/decomposition/sindy.py` |
| `DMD/SINDy/SINDy_Gluhareff_v_0.mlx` | `examples/data_driven/sindy_gluhareff.py` |
| `DMD/SINDy/SINDy_Version_0.mlx` | `examples/data_driven/sindy_examples.py` |
| `DMD/SINDy/SINDy_Version_0_backup.mlx` | examples/data_driven/sindy_examples.py |
| `DMD/SINDy/SINDy_Version_1.mlx` | examples/data_driven/sindy_examples.py |
| `DMD/SINDy/SINDy_Version_2.mlx` | examples/data_driven/sindy_examples.py |
| `DMD/SINDy/SINDy_Version_2_5.mlx` | examples/data_driven/sindy_examples.py |
| `DMD/SINDy/SINDy_Version_2_6.mlx` | `examples/data_driven/sindy_examples.py` |
| `DMD/SINDy/STLS.mlx` | `src/unicodes/decomposition/sindy.py` |
| `DMD/SINDy/STR.mlx` | Unfinished STRidge attempt; src/unicodes/decomposition/sindy.py (`stridge`) |
| `DMD/SINDy/STRidge.mlx` | `src/unicodes/decomposition/sindy.py` |
| `DMD/SINDy/TVRegDiff.m` | `examples/data_driven/sindy_examples.py`, `examples/sindy_examples/ex02_lorenz.py`, `src/unicodes/numerics.py` |
| `DMD/SINDy/generateLibrary.mlx` | `src/unicodes/decomposition/sindy.py` |
| `DMD/SINDy/generateLibraryList.mlx` | `src/unicodes/decomposition/sindy.py` |
| `DMD/SINDy/lorenz.mlx` | `examples/sindy_examples/ex02_lorenz.py`, `src/unicodes/decomposition/sindy.py` |
| `DMD/SINDy/simLorenz.mlx` | `examples/data_driven/sindy_examples.py` |
| `DMD/SINDy/simVanderpol.mlx` | `examples/data_driven/sindy_examples.py` |
| `DMD/SINDy/vanderpol.mlx` | `src/unicodes/decomposition/sindy.py` |
| `DMD/convert_mat_two_bin_for_combust_data.mlx` | examples/data_driven/combustor_dmd.py and src/unicodes/decomposition/dmd.py |
| `DMD/convert_mat_two_bin_for_pressure.mlx` | examples/data_driven/combustor_dmd.py and src/unicodes/decomposition/dmd.py |
| `DMD/create_bin_from_mat_data.mlx` | `examples/data_driven/combustor_dmd.py` |
| `DMD/dmd_psd_playground.mlx` | `examples/data_driven/dmd_spectrum.py`, `src/unicodes/decomposition/spectral.py` |
| `DMD/export_development.mlx` | examples/data_driven/combustor_dmd.py and src/unicodes/decomposition/dmd.py |
| `DMD/export_driver.mlx` | `examples/data_driven/combustor_dmd.py` |
| `DMD/export_snapshots.mlx` | `examples/data_driven/combustor_dmd.py` |
| `DMD/importGridFile.m` | `src/unicodes/io/tecplot.py` |
| `DMD/importTecASCIIdata.m` | `src/unicodes/io/tecplot.py` |
| `DMD/load_driver.mlx` | `examples/data_driven/combustor_dmd.py` |
| `DMD/load_snapshots.mlx` | `examples/data_driven/combustor_dmd.py` |
| `DMD/load_snapshots_P.mlx` | examples/data_driven/combustor_dmd.py and src/unicodes/decomposition/dmd.py |
| `DMD/load_snapshots_T.mlx` | examples/data_driven/combustor_dmd.py and src/unicodes/decomposition/dmd.py |
| `DMD/load_snapshots_U.mlx` | examples/data_driven/combustor_dmd.py and src/unicodes/decomposition/dmd.py |
| `DMD/load_snapshots_V.mlx` | examples/data_driven/combustor_dmd.py and src/unicodes/decomposition/dmd.py |
| `DMD/load_snapshots_Y_CH4.mlx` | examples/data_driven/combustor_dmd.py and src/unicodes/decomposition/dmd.py |
| `DMD/load_snapshots_Y_CO2.mlx` | examples/data_driven/combustor_dmd.py and src/unicodes/decomposition/dmd.py |
| `DMD/load_snapshots_Y_H2O.mlx` | examples/data_driven/combustor_dmd.py and src/unicodes/decomposition/dmd.py |
| `DMD/load_snapshots_Y_O2.mlx` | examples/data_driven/combustor_dmd.py and src/unicodes/decomposition/dmd.py |
| `DMD/psd_development.mlx` | `src/unicodes/decomposition/spectral.py` |
| `DMD/pyDMD Data Analysis/analyze_spectrum.mlx` | `examples/data_driven/pydmd_results.py` |
| `DMD/pyDMD Data Analysis/data_visualization.mlx` | `examples/data_driven/pydmd_results.py` |
| `DMD/pyDMD Data Analysis/data_visualization_10_30_2022.mlx` | examples/data_driven/pydmd_results.py |
| `DMD/pyDMD Data Analysis/getModeData.mlx` | `examples/data_driven/pydmd_results.py` |
| `DMD/pyDMD Data Analysis/local under sampled py dmd results 15-5-2023/analyze_spectrum.mlx` | `examples/data_driven/pydmd_results.py` |
| `DMD/pyDMD Data Analysis/local under sampled py dmd results 15-5-2023/data_visualization.mlx` | `examples/data_driven/pydmd_results.py` |
| `DMD/pyDMD Data Analysis/local under sampled py dmd results 15-5-2023/getModeData.mlx` | `examples/data_driven/pydmd_results.py` |
| `DMD/pyDMD Data Analysis/local under sampled py dmd results 15-5-2023/plot_spectrum.mlx` | `examples/data_driven/pydmd_results.py` |
| `DMD/pyDMD Data Analysis/local under sampled py dmd results 15-5-2023/read_data.mlx` | `examples/data_driven/pydmd_results.py` |
| `DMD/pyDMD Data Analysis/messy_graphs_AAAAAAA.mlx` | `examples/data_driven/pydmd_results.py` |
| `DMD/pyDMD Data Analysis/plot_spectrum.mlx` | `examples/data_driven/pydmd_results.py` |
| `DMD/pyDMD Data Analysis/read_data.mlx` | `examples/data_driven/pydmd_results.py` |
| `DMD/pyDMD Data Analysis/read_data_test.mlx` | examples/data_driven/pydmd_results.py |
| `DMD/pyDMD Data Analysis/read_modes.mlx` | `examples/data_driven/pydmd_results.py` |
| `DMD/pyDMD Data Analysis/read_modes_reverse.mlx` | examples/data_driven/pydmd_results.py |
| `DMD/pythonComparisonDMD.mlx` | examples/data_driven/combustor_dmd.py and src/unicodes/decomposition/dmd.py |
| `DMD/reconstruct.m` | examples/data_driven/combustor_dmd.py and src/unicodes/decomposition/dmd.py |
| `DMD/saveModesForExport.m` | `examples/data_driven/combustor_dmd.py` |
| `DMD/undersample_data_for_combustor.mlx` | `examples/data_driven/combustor_dmd.py` |

## GIT Test

| MATLAB | Python / note |
|---|---|
| `GIT Test/Hello_world.mlx` | Not ported: git test script |

## Math 290

| MATLAB | Python / note |
|---|---|
| `Math 290/Linear_Algebra_Experimentation.mlx` | `examples/math_courses/linear_algebra.py` |
| `Math 290/find_basis.mlx` | `examples/math_courses/linear_algebra.py` |
| `Math 290/linear_algebra_test_2.mlx` | `examples/math_courses/linear_algebra.py` |
| `Math 290/linear_algebra_test_final.mlx` | `examples/math_courses/linear_algebra.py` |

## Math 590

| MATLAB | Python / note |
|---|---|
| `Math 590/linear_algebra_hw_7.mlx` | `examples/math_courses/linear_algebra.py` |

## Math 650

| MATLAB | Python / note |
|---|---|
| `Math 650/HW_1.mlx` | examples/math_courses/nonlinear_dynamics.py |
| `Math 650/HW_2.mlx` | examples/math_courses/nonlinear_dynamics.py |
| `Math 650/Logistic_Map_Unstable_Fixed_Points.mlx` | `examples/math_courses/nonlinear_dynamics.py` |
| `Math 650/itterative_map_plotter.mlx` | `examples/math_courses/nonlinear_dynamics.py` |

## Spring 2025/AE 700 vehilce design requirements

| MATLAB | Python / note |
|---|---|
| `Spring 2025/AE 700 vehilce design requirements/DescentFunction.mlx` | Copy of the AE 725 function: src/unicodes/optimize/ and src/unicodes/numerics.py |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/Detect_IR_v3.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/Detect_Radar.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/Fuze.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/Improved_proportional_guidance_algorithm.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/Improved_proportional_guidance_algorithm_using_sensor.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/Improved_proportional_guidance_algorithm_using_sensor_v1.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/MissileGuidanceEnv.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/NN_RL_V0.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/NN_guidance.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/NN_guidance_v0.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/NN_guidance_v1.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/NN_guidance_v10.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/NN_guidance_v11.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/NN_guidance_v12.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/NN_guidance_v13.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/NN_guidance_v14.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/NN_guidance_v14_works_well.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/NN_guidance_v15.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/NN_guidance_v16.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/NN_guidance_v17.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/NN_guidance_v17_B.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/NN_guidance_v17_B2.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/NN_guidance_v17_intresting.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/NN_guidance_v18.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/NN_guidance_v18_B.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/NN_guidance_v18_B2.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/NN_guidance_v19.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/NN_guidance_v2.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/NN_guidance_v20.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/NN_guidance_v20_PID_poroportional.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/NN_guidance_v21.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/NN_guidance_v21_longtrain.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/NN_guidance_v22.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/NN_guidance_v23.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/NN_guidance_v24.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/NN_guidance_v24_faster_simulate.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/NN_guidance_v25.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/NN_guidance_v25_integrated_fusing.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/NN_guidance_v25_self_fusing.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/NN_guidance_v26.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/NN_guidance_v3.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/NN_guidance_v4.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/NN_guidance_v5.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/NN_guidance_v6.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/NN_guidance_v6_works.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/NN_guidance_v7.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/NN_guidance_v8.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/NN_guidance_v9.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/Pd_estimator.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/Pd_estimator_V2.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/Project_4_algegbra.mlx` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/Project_4_algegbra_alternative_0.mlx` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/Reference.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/Try1.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/Try2_Proportional.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/Try3_Pure.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/estimatePd.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/mySimpleNN.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/outputs.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/runner_final_V0.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/simAndPlot_requested.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project 4/simulateAndPlot.m` | Not ported: missile guidance project (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Project1Runthrough_5.m` | `examples/ae700_sensor_design/project1_optical.py`, `src/unicodes/remote_sensing.py` |
| `Spring 2025/AE 700 vehilce design requirements/Project_2_algegbra.mlx` | `examples/ae700_sensor_design/project2_thermal.py` |
| `Spring 2025/AE 700 vehilce design requirements/Project_3_algrebra.mlx` | `examples/ae700_sensor_design/project3_sar.py`, `src/unicodes/remote_sensing.py` |
| `Spring 2025/AE 700 vehilce design requirements/Radar_reflected_returns_2_object_V0.mlx` | Not ported: weapons effects / proportional navigation (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Radar_reflected_returns_2_object_V1.mlx` | Not ported: weapons effects / proportional navigation (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/Radar_reflected_returns_2_object_V2.mlx` | Not ported: weapons effects / proportional navigation (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/computeJacobian.mlx` | Copy of the AE 725 function: src/unicodes/optimize/ and src/unicodes/numerics.py |
| `Spring 2025/AE 700 vehilce design requirements/continuous_simulated_anealing.mlx` | Copy of the AE 725 function: src/unicodes/optimize/ and src/unicodes/numerics.py |
| `Spring 2025/AE 700 vehilce design requirements/final_exam_problem_1.mlx` | Not ported: symbolic derivation (IFOV rate) |
| `Spring 2025/AE 700 vehilce design requirements/final_exam_problem_2.mlx` | `src/unicodes/remote_sensing.py` |
| `Spring 2025/AE 700 vehilce design requirements/final_exam_problem_3.mlx` | examples/ae700_sensor_design/final_exam.py |
| `Spring 2025/AE 700 vehilce design requirements/final_exam_problem_4.mlx` | Not ported: weapons effects / proportional navigation (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/fragmentation_area.mlx` | Not ported: weapons effects / proportional navigation (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/genetic_algorithm_optimization.mlx` | Copy of the AE 725 function: src/unicodes/optimize/ and src/unicodes/numerics.py |
| `Spring 2025/AE 700 vehilce design requirements/inexactStepSize.mlx` | Copy of the AE 725 function: src/unicodes/optimize/ and src/unicodes/numerics.py |
| `Spring 2025/AE 700 vehilce design requirements/lambda_max_blackbody.mlx` | `src/unicodes/remote_sensing.py` |
| `Spring 2025/AE 700 vehilce design requirements/optimized_project_1_V_2.mlx` | `examples/ae700_sensor_design/project1_optical.py` |
| `Spring 2025/AE 700 vehilce design requirements/probability_of_damage.mlx` | Not ported: weapons effects / proportional navigation (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/project_1_v_1_0.m` | `examples/ae700_sensor_design/project1_optical.py`, `src/unicodes/remote_sensing.py` |
| `Spring 2025/AE 700 vehilce design requirements/pure_proportional_navigation.mlx` | Not ported: weapons effects / proportional navigation (excluded) |
| `Spring 2025/AE 700 vehilce design requirements/spectral_radiance_labertian_blackbody.mlx` | `src/unicodes/remote_sensing.py` |
| `Spring 2025/AE 700 vehilce design requirements/untitled.mlx` | Not ported: empty file |

## Spring 2025/AE 846 Advanced CFD

| MATLAB | Python / note |
|---|---|
| `Spring 2025/AE 846 Advanced CFD/FR_coefs.m` | `src/unicodes/cfd/flux_reconstruction.py` |
| `Spring 2025/AE 846 Advanced CFD/HW_1_question_1.mlx` | examples/ae746_cfd/ae846_homework.py |
| `Spring 2025/AE 846 Advanced CFD/HW_1_question_4.mlx` | examples/ae746_cfd/ae846_homework.py |
| `Spring 2025/AE 846 Advanced CFD/Reacting_flow_1D_equations_derivation.mlx` | Not ported: unfinished symbolic derivation |

## ae 360

| MATLAB | Python / note |
|---|---|
| `ae 360/HW_5.mlx` | examples/ae360_orbital/transfers.py |
| `ae 360/Joshua_Poznanski_HW3.mlx` | Not ported: empty file |
| `ae 360/Joshua_Poznanski_HW4.mlx` | `examples/ae360_orbital/euler_vs_kepler.py` |
| `ae 360/Joshua_Poznanski_HW6.mlx` | `examples/ae360_orbital/kepler_equation.py` |
| `ae 360/Joshua_Poznanski_HW9.m` | `examples/ae360_orbital/euler_vs_kepler.py` |
| `ae 360/Joshua_Poznanski_HW9.mlx` | `examples/ae360_orbital/euler_vs_kepler.py` |
| `ae 360/ae_360_analysis.mlx` | `examples/ae360_orbital/balloon_flight_data.py`, `examples/ae360_orbital/euler_vs_kepler.py` |
| `ae 360/coe2rv.mlx` | src/unicodes/orbital.py |
| `ae 360/debug_orbital_conversions.mlx` | src/unicodes/orbital.py |
| `ae 360/driver.mlx` | `examples/ae360_orbital/transfers.py`, `src/unicodes/orbital.py` |
| `ae 360/driver_final.mlx` | `examples/ae360_orbital/euler_vs_kepler.py` |
| `ae 360/getTranserOrbitHoffmamn.mlx` | `examples/ae360_orbital/transfers.py`, `src/unicodes/orbital.py` |
| `ae 360/hoffman.mlx` | `examples/ae360_orbital/transfers.py`, `src/unicodes/orbital.py` |
| `ae 360/integrate2Body.mlx` | `examples/ae360_orbital/euler_vs_kepler.py` |
| `ae 360/newtons_method_experimentation.mlx` | `examples/ae360_orbital/kepler_equation.py` |
| `ae 360/rendezvous.mlx` | `examples/ae360_orbital/transfers.py`, `src/unicodes/orbital.py` |
| `ae 360/rv2coe.mlx` | src/unicodes/orbital.py |

## autumn 2023/AE 510 Materials and Processes

| MATLAB | Python / note |
|---|---|
| `autumn 2023/AE 510 Materials and Processes/HW_4.mlx` | `examples/ae510_materials/material_selection.py` |
| `autumn 2023/AE 510 Materials and Processes/Material.m` | src/unicodes/materials.py |
| `autumn 2023/AE 510 Materials and Processes/save_materials.mlx` | `examples/ae510_materials/material_selection.py` |
| `autumn 2023/AE 510 Materials and Processes/supplimental_project_one-ENGR-85N8KQ3-2.mlx` | examples/ae510_materials/material_selection.py (drafts of the same script) |
| `autumn 2023/AE 510 Materials and Processes/supplimental_project_one-ENGR-85N8KQ3.mlx` | examples/ae510_materials/material_selection.py (drafts of the same script) |
| `autumn 2023/AE 510 Materials and Processes/supplimental_project_one.mlx` | `examples/ae510_materials/material_selection.py` |

## autumn 2023/AE 521 Aircraft Design

| MATLAB | Python / note |
|---|---|
| `autumn 2023/AE 521 Aircraft Design/BreguetEndurance.mlx` | `src/unicodes/aero/performance.py` |
| `autumn 2023/AE 521 Aircraft Design/BreguetRange.mlx` | `src/unicodes/aero/performance.py` |
| `autumn 2023/AE 521 Aircraft Design/Daedalus_preliminary_Sizing_V_0.mlx` | `examples/ae521_aircraft_design/daedalus_constraint_diagram.py` |
| `autumn 2023/AE 521 Aircraft Design/DensityRatio.mlx` | `src/unicodes/atmosphere.py` |
| `autumn 2023/AE 521 Aircraft Design/DragCoefficientZeroLift.mlx` | `src/unicodes/aero/performance.py` |
| `autumn 2023/AE 521 Aircraft Design/DragPolar_parabolic.mlx` | src/unicodes/aero/wing.py (`parabolic_drag_polar`) |
| `autumn 2023/AE 521 Aircraft Design/EquivalentParasiteArea.mlx` | `src/unicodes/aero/performance.py` |
| `autumn 2023/AE 521 Aircraft Design/FuelFraction_Endurance.mlx` | `src/unicodes/aero/performance.py` |
| `autumn 2023/AE 521 Aircraft Design/FuelFraction_Range.mlx` | `src/unicodes/aero/performance.py` |
| `autumn 2023/AE 521 Aircraft Design/FuelFraction_Suggested.mlx` | Reads Roskam's CSV; the table values are in examples/ae521_aircraft_design/daedalus_weight_sizing.py |
| `autumn 2023/AE 521 Aircraft Design/Internal_Combustion_Engine_Predictor_Equation.mlx` | `src/unicodes/aero/performance.py` |
| `autumn 2023/AE 521 Aircraft Design/Power_Lapse_With_Alt_IC_engine.mlx` | `src/unicodes/aero/performance.py` |
| `autumn 2023/AE 521 Aircraft Design/Report_3_voyager.mlx` | `examples/ae521_aircraft_design/voyager_power.py` |
| `autumn 2023/AE 521 Aircraft Design/Report_4_DC_3.mlx` | `examples/ae521_aircraft_design/dc3_schrenk.py` |
| `autumn 2023/AE 521 Aircraft Design/Required_Power_Subsonic.mlx` | `src/unicodes/aero/performance.py` |
| `autumn 2023/AE 521 Aircraft Design/Schrenk_Lift_aproximation.mlx` | src/unicodes/aero/wing.py (`schrenk_lift_distribution`) |
| `autumn 2023/AE 521 Aircraft Design/StallSpeed.mlx` | src/unicodes/aero/wing.py (`stall_speed`) |
| `autumn 2023/AE 521 Aircraft Design/T2W_ClimbGradient_FAR25.mlx` | `src/unicodes/aero/performance.py` |
| `autumn 2023/AE 521 Aircraft Design/T2W_Climb_FAR25.mlx` | `src/unicodes/aero/performance.py` |
| `autumn 2023/AE 521 Aircraft Design/T2W_CruiseSpeed.mlx` | `src/unicodes/aero/performance.py` |
| `autumn 2023/AE 521 Aircraft Design/T2W_Takeoff_FAR25.mlx` | `src/unicodes/aero/performance.py` |
| `autumn 2023/AE 521 Aircraft Design/T2W_TimeToClimb.mlx` | `src/unicodes/aero/performance.py` |
| `autumn 2023/AE 521 Aircraft Design/W2S_Landing_FAR25.mlx` | `src/unicodes/aero/performance.py` |
| `autumn 2023/AE 521 Aircraft Design/W2S_TimeToClimb.mlx` | `src/unicodes/aero/performance.py` |
| `autumn 2023/AE 521 Aircraft Design/Weight_Crew.mlx` | `src/unicodes/aero/performance.py` |
| `autumn 2023/AE 521 Aircraft Design/Weight_Empty.mlx` | `src/unicodes/aero/performance.py` |
| `autumn 2023/AE 521 Aircraft Design/Weight_PassengerAndBaggage.mlx` | `src/unicodes/aero/performance.py` |
| `autumn 2023/AE 521 Aircraft Design/Weight_Sizing_Daedalus_V_0.mlx` | `examples/ae521_aircraft_design/daedalus_weight_sizing.py` |
| `autumn 2023/AE 521 Aircraft Design/WettedArea.mlx` | `src/unicodes/aero/performance.py` |
| `autumn 2023/AE 521 Aircraft Design/prelim_code_to_to_run.m` | `examples/ae521_aircraft_design/daedalus_constraint_diagram.py` |
| `autumn 2023/AE 521 Aircraft Design/vehicle_weight_sizing_1.m` | `examples/ae521_aircraft_design/hermes_weight_sizing.py` |
| `autumn 2023/AE 521 Aircraft Design/weight_Sizing_Hermes_V_0.mlx` | `examples/ae521_aircraft_design/hermes_weight_sizing.py` |

## autumn 2023/AE 721 Missile Design

| MATLAB | Python / note |
|---|---|
| `autumn 2023/AE 721 Missile Design/BendingFrequency_Body_first.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/HIB.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/L_over_D_for_qbars.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/L_over_D_for_qbars_V2-ENGR-HWC66X3.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/L_over_D_for_qbars_V2.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/L_over_D_for_qbars_V3.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/L_over_D_for_qbars_V4.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/L_over_D_for_qbars_V4_fleeman_reproducer.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/L_over_D_for_qbars_V5_altitudes.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/Planar_Surface_Drag_Verification.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/Plot_Wave_drag_benchmark_Round.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/RocketBullet.m` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/Spline_For_transonic_body_drag.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/aerodynamic_Center_Body.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/aerodynamic_Center_Body_Alone.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/aerodynamic_center_planar_surface.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/ballistic_trajectory_simulator.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/ballistic_trajectory_simulator_V2.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/ballistic_trajectory_simulator_V2_for_joesRound.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/ballistic_trajectory_simulator_V3.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/ballistic_trajectory_simulator_V3_ASAT.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/ballistic_trajectory_simulator_V3_RAIDER_Bass.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/ballistic_trajectory_simulator_V3_joes_round.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/ballistic_trajectory_simulator_range_angle_computer.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/componentBuildup_aerodynamicCenter.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/dragCoeff_friction_surface.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/dragCoeff_wave_surface.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/dragCoeff_wave_zeroLift_bluntedNose.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/dragCoeff_zeroLift_body.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/dragCoeff_zeroLift_body_smoothed_transonic.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/dragCoeff_zeroLift_surface.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/drag_body.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/file to send to bass group/RocketBullet.m` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/file to send to bass group/ballistic_trajectory_simulator_V3_RAIDER_Bass.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/file to send to bass group/dragCoeff_wave_zeroLift_bluntedNose.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/file to send to bass group/dragCoeff_zeroLift_body.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/file to send to bass group/dragCoeff_zeroLift_body_smoothed_transonic.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/file to send to bass group/drag_body.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/file to send to bass group/save_bullet_info.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/flare_aerodynamic_center.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/getAirFlowAngles.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/getFlightPathAngles.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/hinge_moment.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/lift_to_drag_ratio_body.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/normalForceCoeff_liftingBody.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/normalForceCoeff_planar_surface.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/normalForce_Flare_withRespectTo_alpha.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/normalForce_withRespectTo_alpha.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/normalForce_withRespectTo_alpha_planar_surface.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/plot_C_D_0_surface_vs_mach.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/plot_C_N_vs_alpha_for_lifting_body.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/plot_aerodynamicCeneteroverNoseLength_v_s_alpha.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/plot_lift_vs_drag_body_example.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/plot_tail_area_sizing.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/plot_x_AC_planar_surface.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/polt_C_D_0_body_vs_mach.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/save_bullet_info.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/simulate_balistic_Trajectory.mlx` | Not ported: missile design (excluded as requested) |
| `autumn 2023/AE 721 Missile Design/tail_area_sizing.mlx` | Not ported: missile design (excluded as requested) |

## autumn 2023/AE 746 CFD

| MATLAB | Python / note |
|---|---|
| `autumn 2023/AE 746 CFD/HW_1_func.m` | `examples/ae746_cfd/homework.py` |
| `autumn 2023/AE 746 CFD/Homework_1.mlx` | `examples/ae746_cfd/homework.py` |
| `autumn 2023/AE 746 CFD/Homework_2_numericalSolver.mlx` | `examples/ae746_cfd/homework.py` |
| `autumn 2023/AE 746 CFD/Homework_2_p3_2.mlx` | examples/ae746_cfd/homework.py |
| `autumn 2023/AE 746 CFD/Homework_2_p3_7.mlx` | examples/ae746_cfd/homework.py |
| `autumn 2023/AE 746 CFD/RK4.mlx` | `src/unicodes/ode.py` |
| `autumn 2023/AE 746 CFD/SSP_RK3.mlx` | `src/unicodes/ode.py` |
| `autumn 2023/AE 746 CFD/explicitEuler.mlx` | `src/unicodes/ode.py` |
| `autumn 2023/AE 746 CFD/modifiedEuler.mlx` | `src/unicodes/ode.py` |
| `autumn 2023/AE 746 CFD/project 2/Euler1DResult.m` | `examples/ae746_cfd/project2_nozzle.py` |
| `autumn 2023/AE 746 CFD/project 2/ExactSolu.mlx` | `src/unicodes/cfd/nozzle.py` |
| `autumn 2023/AE 746 CFD/project 2/Exact_Nozzle.mlx` | `examples/ae746_cfd/project2_nozzle.py` |
| `autumn 2023/AE 746 CFD/project 2/Initialize.mlx` | `examples/ae746_cfd/project4_euler2d.py`, `src/unicodes/cfd/nozzle.py` |
| `autumn 2023/AE 746 CFD/project 2/MachArea.mlx` | `src/unicodes/cfd/nozzle.py` |
| `autumn 2023/AE 746 CFD/project 2/Project_2_MUSCL.mlx` | `examples/ae746_cfd/project2_nozzle.py`, `src/unicodes/cfd/nozzle.py` |
| `autumn 2023/AE 746 CFD/project 2/Project_2_Main_V_0.mlx` | examples/ae746_cfd/project2_nozzle.py and src/unicodes/cfd/nozzle.py (earlier versions of the same solver) |
| `autumn 2023/AE 746 CFD/project 2/Project_2_Main_V_1.mlx` | examples/ae746_cfd/project2_nozzle.py and src/unicodes/cfd/nozzle.py (earlier versions of the same solver) |
| `autumn 2023/AE 746 CFD/project 2/Project_2_Main_V_1_BACKUP_0.mlx` | examples/ae746_cfd/project2_nozzle.py and src/unicodes/cfd/nozzle.py (earlier versions of the same solver) |
| `autumn 2023/AE 746 CFD/project 2/Project_2_Main_V_1_working_Bakcup.mlx` | examples/ae746_cfd/project2_nozzle.py and src/unicodes/cfd/nozzle.py (earlier versions of the same solver) |
| `autumn 2023/AE 746 CFD/project 2/Project_2_Main_V_2.mlx` | examples/ae746_cfd/project2_nozzle.py and src/unicodes/cfd/nozzle.py (earlier versions of the same solver) |
| `autumn 2023/AE 746 CFD/project 2/Project_2_Main_V_3.mlx` | examples/ae746_cfd/project2_nozzle.py and src/unicodes/cfd/nozzle.py (earlier versions of the same solver) |
| `autumn 2023/AE 746 CFD/project 2/Project_2_Main_V_3_5.mlx` | examples/ae746_cfd/project2_nozzle.py and src/unicodes/cfd/nozzle.py (earlier versions of the same solver) |
| `autumn 2023/AE 746 CFD/project 2/Project_2_Main_V_4.mlx` | examples/ae746_cfd/project2_nozzle.py and src/unicodes/cfd/nozzle.py (earlier versions of the same solver) |
| `autumn 2023/AE 746 CFD/project 2/Project_2_Main_V_4_5.mlx` | examples/ae746_cfd/project2_nozzle.py and src/unicodes/cfd/nozzle.py (earlier versions of the same solver) |
| `autumn 2023/AE 746 CFD/project 2/Project_2_Main_V_4_6.mlx` | `examples/ae746_cfd/project2_nozzle.py` |
| `autumn 2023/AE 746 CFD/project 2/Property_testing.mlx` | examples/ae746_cfd/project2_nozzle.py and src/unicodes/cfd/nozzle.py (earlier versions of the same solver) |
| `autumn 2023/AE 746 CFD/project 2/area.mlx` | `src/unicodes/cfd/nozzle.py` |
| `autumn 2023/AE 746 CFD/project 2/characteristicBC_algebraicValidation.mlx` | examples/ae746_cfd/project2_nozzle.py and src/unicodes/cfd/nozzle.py (earlier versions of the same solver) |
| `autumn 2023/AE 746 CFD/project 2/getFlowProperties.m` | examples/ae746_cfd/project2_nozzle.py and src/unicodes/cfd/nozzle.py (earlier versions of the same solver) |
| `autumn 2023/AE 746 CFD/project 3/arcLenParamterizationPlaybox.mlx` | `examples/ae746_cfd/project3_mesh.py` |
| `autumn 2023/AE 746 CFD/project 3/project_3_main_v_0.mlx` | examples/ae746_cfd/project3_mesh.py |
| `autumn 2023/AE 746 CFD/project 3/project_3_main_v_0_backup.mlx` | examples/ae746_cfd/project3_mesh.py |
| `autumn 2023/AE 746 CFD/project 3/project_3_main_v_1.mlx` | `examples/ae746_cfd/project3_mesh.py` |
| `autumn 2023/AE 746 CFD/project 3/thirdArc.mlx` | `examples/ae746_cfd/project3_mesh.py` |
| `autumn 2023/AE 746 CFD/project 4/Initialize.mlx` | `examples/ae746_cfd/project4_euler2d.py`, `src/unicodes/cfd/nozzle.py` |
| `autumn 2023/AE 746 CFD/project 4/Project_4_V_0.mlx` | examples/ae746_cfd/project4_euler2d.py and src/unicodes/cfd/ |
| `autumn 2023/AE 746 CFD/project 4/Project_4_V_1.mlx` | `examples/ae746_cfd/project4_euler2d.py` |
| `autumn 2023/AE 746 CFD/project 4/Project_4_V_1_Backup.m` | examples/ae746_cfd/project4_euler2d.py and src/unicodes/cfd/ |
| `autumn 2023/AE 746 CFD/project 4/Project_4_V_1_GPU_0.mlx` | examples/ae746_cfd/project4_euler2d.py and src/unicodes/cfd/ |
| `autumn 2023/AE 746 CFD/project 4/configure_parameters_and_BC.mlx` | `examples/ae746_cfd/project4_euler2d.py` |
| `autumn 2023/AE 746 CFD/project 4/getFlowProperties.mlx` | examples/ae746_cfd/project4_euler2d.py and src/unicodes/cfd/ |
| `autumn 2023/AE 746 CFD/project 4/getQ.mlx` | `examples/ae746_cfd/project4_euler2d.py` |
| `autumn 2023/AE 746 CFD/project 4/untitled4.mlx` | examples/ae746_cfd/project4_euler2d.py and src/unicodes/cfd/ |
| `autumn 2023/AE 746 CFD/project_1/Project_1.mlx` | `examples/ae746_cfd/project1_advection.py`, `src/unicodes/cfd/advection.py` |
| `autumn 2023/AE 746 CFD/project_1/Project_1_part_2.mlx` | `examples/ae746_cfd/project1_advection.py`, `src/unicodes/cfd/advection.py` |
| `autumn 2023/AE 746 CFD/test_solving_hw_3.mlx` | `examples/ae746_cfd/homework.py` |

## autumn 2023/orthogonal_vectors.mlx

| MATLAB | Python / note |
|---|---|
| `autumn 2023/orthogonal_vectors.mlx` | `examples/misc/small_scripts.py` |

## data

| MATLAB | Python / note |
|---|---|
| `data/sampleO_parallelized.m` | src/unicodes/io/tecplot.py |

## data and processed data

| MATLAB | Python / note |
|---|---|
| `data and processed data/InOutFormat.m` | src/unicodes/io/tecplot.py |
| `data and processed data/OutputTecASCIIdata.m` | `src/unicodes/io/tecplot.py` |
| `data and processed data/importGridFile.m` | `src/unicodes/io/tecplot.py` |
| `data and processed data/importTecASCIIdata.m` | `src/unicodes/io/tecplot.py` |
| `data and processed data/sampleIO.m` | src/unicodes/io/tecplot.py |
| `data and processed data/sampleO_comented.m` | src/unicodes/io/tecplot.py |
| `data and processed data/sampleO_parallelized.m` | src/unicodes/io/tecplot.py |

## first semester/ae 245

| MATLAB | Python / note |
|---|---|
| `first semester/ae 245/Analysis.m` | `examples/data_driven/pydmd_results.py` |
| `first semester/ae 245/ae245quickanalysis.mlx` | `examples/ae245_flight_logs/flight_log_analysis.py` |
| `first semester/ae 245/ae_245_analysis.mlx` | `examples/ae245_flight_logs/flight_log_analysis.py` |
| `first semester/ae 245/ae_245_hw_3.mlx` | `examples/ae245_flight_logs/flight_log_analysis.py` |
| `first semester/ae 245/deg_min_sec_to_deci_deg.m` | examples/ae245_flight_logs/flight_log_analysis.py (`lla_to_flat` replaces geodetic2ned) |
| `first semester/ae 245/lab2_1.mlx` | `examples/ae245_flight_logs/flight_log_analysis.py` |
| `first semester/ae 245/lab2_1_script.m` | examples/ae245_flight_logs/flight_log_analysis.py (`lla_to_flat` replaces geodetic2ned) |

## first semester/diff eq

| MATLAB | Python / note |
|---|---|
| `first semester/diff eq/diff_eq_homework_9.m` | `examples/math_courses/differential_equations.py` |

## first semester/statics and dynamics

| MATLAB | Python / note |
|---|---|
| `first semester/statics and dynamics/Statics_and_dynamics_group_porject_rollercoaster.mlx` | `examples/misc/rollercoaster.py` |

## gluhareff pressure jet

| MATLAB | Python / note |
|---|---|
| `gluhareff pressure jet/Bramlette_ku_0099D_14625_DATA_2IR.m` | `examples/gluhareff_pressure_jet/ir_temperature_contours.py` |
| `gluhareff pressure jet/Bramlette_ku_0099D_14625_DATA_3gn.m` | `examples/gluhareff_pressure_jet/pressure_jet_sizing.py` |

## jacksons 211 programs

| MATLAB | Python / note |
|---|---|
| `jacksons 211 programs/AE360_Torok_HW2.m` | `examples/jackson_torok/ae360_hw02.py` |
| `jacksons 211 programs/AE360_Torok_HW4.m` | `examples/jackson_torok/ae360_hw04.py` |
| `jacksons 211 programs/AE360_Torok_HW6.m` | `examples/jackson_torok/ae360_hw06.py` |
| `jacksons 211 programs/AE360_Torok_HW9.m` | `examples/jackson_torok/ae360_hw09.py` |
| `jacksons 211 programs/CompiledCode771.m` | `examples/jackson_torok/rocket_propulsion_771.py` |
| `jacksons 211 programs/HW5 Functions/Torok_Jackson_HW5_MATLAB.m` | `examples/jackson_torok/ae211_hw05.py` |
| `jacksons 211 programs/HW5 Functions/distance.m` | `examples/jackson_torok/ae211_hw05.py` |
| `jacksons 211 programs/HW5 Functions/height.m` | `examples/jackson_torok/ae211_hw05.py` |
| `jacksons 211 programs/HW5 Functions/num_grain.m` | `examples/jackson_torok/ae211_hw05.py` |
| `jacksons 211 programs/HW5 Functions/plotndfhs.m` | `examples/jackson_torok/ae211_hw05.py` |
| `jacksons 211 programs/Matlab_Quiz5_Jackson_Torok.m` | `examples/jackson_torok/ae211_quiz5.py` |
| `jacksons 211 programs/Torok_Jackson_Exam2.m` | `examples/jackson_torok/ae211_exams.py` |
| `jacksons 211 programs/Torok_Jackson_HW10_MATLAB.m` | `examples/jackson_torok/ae211_hw10.py` |
| `jacksons 211 programs/Torok_Jackson_HW1_MATLAB.m` | `examples/jackson_torok/ae211_hw01.py` |
| `jacksons 211 programs/Torok_Jackson_HW2_MATLAB.m` | `examples/jackson_torok/ae211_hw02.py` |
| `jacksons 211 programs/Torok_Jackson_HW3_MATLAB.m` | `examples/jackson_torok/ae211_hw03.py` |
| `jacksons 211 programs/Torok_Jackson_HW4_MATLAB.m` | `examples/jackson_torok/ae211_hw04.py` |
| `jacksons 211 programs/Torok_Jackson_HW6_MATLAB.m` | `examples/jackson_torok/ae211_hw06.py` |
| `jacksons 211 programs/Torok_Jackson_HW7_MATLAB.m` | `examples/jackson_torok/ae211_hw07.py` |
| `jacksons 211 programs/Torok_Jackson_HW8_MATLAB.m` | `examples/jackson_torok/ae211_hw08.py` |
| `jacksons 211 programs/Torok_Jackson_HW9_MATLAB.m` | `examples/jackson_torok/ae211_hw09.py` |
| `jacksons 211 programs/Torok__Jackson_Exam1.m` | `examples/jackson_torok/ae211_exams.py` |
| `jacksons 211 programs/Untitled.m` | `examples/jackson_torok/intercept_game.py` |
| `jacksons 211 programs/Untitled4.m` | `examples/jackson_torok/ae211_quiz5.py` |
| `jacksons 211 programs/exam.m` | `examples/jackson_torok/ae211_exams.py` |

## sample_code

| MATLAB | Python / note |
|---|---|
| `sample_code/InOutFormat.m` | src/unicodes/io/tecplot.py |
| `sample_code/OutputTecASCIIdata.m` | `src/unicodes/io/tecplot.py` |
| `sample_code/importGridFile.m` | `src/unicodes/io/tecplot.py` |
| `sample_code/importTecASCIIdata.m` | `src/unicodes/io/tecplot.py` |
| `sample_code/sampleIO.m` | src/unicodes/io/tecplot.py |
| `sample_code/sampleO_parallelized.m` | src/unicodes/io/tecplot.py |

## spring 2023/AE 508

| MATLAB | Python / note |
|---|---|
| `spring 2023/AE 508/HW_11.mlx` | `examples/ae508_structures/hw11_stress_convergence.py` |
| `spring 2023/AE 508/HW_11_Rect_GradientPlot.m` | examples/ae508_structures/hw11_stress_convergence.py / final_exam.py (plot helpers) |
| `spring 2023/AE 508/HW_11_Rect_GradientPlot2.m` | examples/ae508_structures/hw11_stress_convergence.py / final_exam.py (plot helpers) |
| `spring 2023/AE 508/HW_11_Tri_GradientPlot2.mlx` | examples/ae508_structures/hw11_stress_convergence.py / final_exam.py (plot helpers) |
| `spring 2023/AE 508/HW_11_gradient.mlx` | `examples/ae508_structures/hw11_stress_convergence.py` |
| `spring 2023/AE 508/HW_11_tri_import.mlx` | `examples/ae508_structures/hw11_stress_convergence.py` |
| `spring 2023/AE 508/HW_6.mlx` | examples/ae508_structures/homework.py |
| `spring 2023/AE 508/HW_7.mlx` | examples/ae508_structures/homework.py |
| `spring 2023/AE 508/HW_7_part_2.mlx` | examples/ae508_structures/homework.py |
| `spring 2023/AE 508/HW_7_part_3.mlx` | examples/ae508_structures/homework.py |
| `spring 2023/AE 508/HW_9.mlx` | examples/ae508_structures/homework.py |
| `spring 2023/AE 508/final exam/problem 3/ConvergencePlots.mlx` | examples/ae508_structures/final_exam.py |
| `spring 2023/AE 508/final exam/problem 3/HW_11_Rect_GradientPlot.mlx` | examples/ae508_structures/hw11_stress_convergence.py / final_exam.py (plot helpers) |
| `spring 2023/AE 508/final exam/problem 3/MinPrincipalGradient.m` | examples/ae508_structures/final_exam.py |
| `spring 2023/AE 508/final exam/problem 4/Gradient_Analysis.mlx` | examples/ae508_structures/final_exam.py |
| `spring 2023/AE 508/final exam/problem 5/Gradient_Investigation.mlx` | examples/ae508_structures/final_exam.py |

## spring 2023/AE 551

| MATLAB | Python / note |
|---|---|
| `spring 2023/AE 551/HW_1_analysis.mlx` | Not ported: log-file plotting only (see examples/ae551_flight_controls/homework.py) |
| `spring 2023/AE 551/HW_2.mlx` | Not ported: symbolic only (see examples/ae551_flight_controls/homework.py docstring) |
| `spring 2023/AE 551/HW_3.mlx` | examples/ae551_flight_controls/homework.py |
| `spring 2023/AE 551/HW_4.mlx` | examples/ae551_flight_controls/homework.py |
| `spring 2023/AE 551/HW_6.mlx` | examples/ae551_flight_controls/homework.py |
| `spring 2023/AE 551/HW_7.mlx` | Not ported: symbolic only (see examples/ae551_flight_controls/homework.py docstring) |
| `spring 2023/AE 551/final project/AE551_Final_Project_TF_SS_Doublet.mlx` | `examples/ae551_flight_controls/final_project.py` |
| `spring 2023/AE 551/final project/MyPiWrap.m` | `examples/ae551_flight_controls/final_project.py`, `src/unicodes/controls.py` |
| `spring 2023/AE 551/final project/mpc1_script.m` | `examples/ae551_flight_controls/final_project.py`, `src/unicodes/controls.py` |
| `spring 2023/AE 551/final project/mpcscript_final.m` | `examples/ae551_flight_controls/final_project.py`, `src/unicodes/controls.py` |

## spring 2023/AE 573

| MATLAB | Python / note |
|---|---|
| `spring 2023/AE 573/EXAM_1_simple.mlx` | `examples/ae573_propulsion/exam_1_simple.py`, `src/unicodes/propulsion_cycles.py` |
| `spring 2023/AE 573/Exam_1.mlx` | Not ported: unfilled exam template (all known values blank); the same separate-exhaust turbofan with values is examples/ae573_propulsion/pratice_quiz_turbofan.py |
| `spring 2023/AE 573/Exam_2_part_1.mlx` | `examples/ae573_propulsion/exam_2_part_1.py`, `src/unicodes/propulsion_cycles.py` |
| `spring 2023/AE 573/Exam_2_part_2.mlx` | `examples/ae573_propulsion/gasdynamics_homework.py` |
| `spring 2023/AE 573/Exam_2_part_3.mlx` | `examples/ae573_propulsion/exam_2_part_3.py`, `src/unicodes/propulsion_cycles.py` |
| `spring 2023/AE 573/FinalExam_1-DESKTOP-BJID2FH.mlx` | Drafts of FinalExam_1: examples/ae573_propulsion/finalexam_1.py |
| `spring 2023/AE 573/FinalExam_1.mlx` | `examples/ae573_propulsion/finalexam_1.py`, `src/unicodes/propulsion_cycles.py` |
| `spring 2023/AE 573/FinalExam_1_v0.mlx` | Drafts of FinalExam_1: examples/ae573_propulsion/finalexam_1.py |
| `spring 2023/AE 573/FinalExam_1_v1.mlx` | Drafts of FinalExam_1: examples/ae573_propulsion/finalexam_1.py |
| `spring 2023/AE 573/FinalExam_1_v2.mlx` | Drafts of FinalExam_1: examples/ae573_propulsion/finalexam_1.py |
| `spring 2023/AE 573/FinalExam_1_v3.mlx` | Drafts of FinalExam_1: examples/ae573_propulsion/finalexam_1.py |
| `spring 2023/AE 573/FinalExam_1_v4.mlx` | Drafts of FinalExam_1: examples/ae573_propulsion/finalexam_1.py |
| `spring 2023/AE 573/FinalExam_1_v5.mlx` | Drafts of FinalExam_1: examples/ae573_propulsion/finalexam_1.py |
| `spring 2023/AE 573/FinalExam_2.mlx` | `examples/ae573_propulsion/gasdynamics_homework.py` |
| `spring 2023/AE 573/FinalExam_3.mlx` | `examples/ae573_propulsion/finalexam_3.py`, `src/unicodes/propulsion_cycles.py` |
| `spring 2023/AE 573/FinalExam_3_v0.mlx` | Draft of FinalExam_3: examples/ae573_propulsion/finalexam_3.py |
| `spring 2023/AE 573/GasLaw.m` | `src/unicodes/thermo/gas_law.py` |
| `spring 2023/AE 573/GasLawHW18.m` | `examples/ae573_propulsion/thermo_hw18.py` |
| `spring 2023/AE 573/GasLaw_Test.mlx` | `examples/ae573_propulsion/thermo_hw18.py` |
| `spring 2023/AE 573/HW_10.mlx` | `examples/ae573_propulsion/hw_10.py`, `src/unicodes/propulsion_cycles.py` |
| `spring 2023/AE 573/HW_15_diffusers.mlx` | `examples/ae573_propulsion/gasdynamics_homework.py` |
| `spring 2023/AE 573/HW_2.mlx` | `examples/ae573_propulsion/gasdynamics_homework.py` |
| `spring 2023/AE 573/HW_3.mlx` | `examples/ae573_propulsion/gasdynamics_homework.py` |
| `spring 2023/AE 573/HW_4.mlx` | `examples/ae573_propulsion/gasdynamics_homework.py` |
| `spring 2023/AE 573/HW_7.mlx` | `examples/ae573_propulsion/hw_7.py`, `examples/ae573_propulsion/hw_8_problem_2.py`, `src/unicodes/propulsion_cycles.py` |
| `spring 2023/AE 573/HW_8_problem_2.mlx` | `examples/ae573_propulsion/hw_8_problem_2.py`, `src/unicodes/propulsion_cycles.py` |
| `spring 2023/AE 573/H_8_part_1.mlx` | `examples/ae573_propulsion/gasdynamics_homework.py` |
| `spring 2023/AE 573/H_8_part_2.mlx` | `examples/ae573_propulsion/gasdynamics_homework.py` |
| `spring 2023/AE 573/H_8_part_3.mlx` | `examples/ae573_propulsion/h_8_part_3.py`, `src/unicodes/propulsion_cycles.py` |
| `spring 2023/AE 573/HomeWork.mlx` | `examples/ae573_propulsion/homework.py`, `src/unicodes/propulsion_cycles.py` |
| `spring 2023/AE 573/InletEntropyRiseNondimensional.mlx` | src/unicodes/gasdynamics.py (`inlet_*`) |
| `spring 2023/AE 573/Inlet_Total_Pressure_Recovery_Relation.m` | src/unicodes/thermo/gas_law.py (`Relation`, `solve_relations`) |
| `spring 2023/AE 573/Inlet_Total_Pressure_Recovery_relation.mlx` | src/unicodes/thermo/gas_law.py (`Relation`, `solve_relations`) |
| `spring 2023/AE 573/Pratice_quiz_turbofan.mlx` | `examples/ae573_propulsion/pratice_quiz_turbofan.py`, `src/unicodes/propulsion_cycles.py` |
| `spring 2023/AE 573/Quiz_2.mlx` | `examples/ae573_propulsion/quiz_2.py`, `src/unicodes/propulsion_cycles.py` |
| `spring 2023/AE 573/Relation.m` | `src/unicodes/thermo/gas_law.py` |
| `spring 2023/AE 573/RelationV0.m` | src/unicodes/thermo/gas_law.py (`Relation`, `solve_relations`) |
| `spring 2023/AE 573/Relation_Solver.mlx` | src/unicodes/thermo/gas_law.py (`Relation`, `solve_relations`) |
| `spring 2023/AE 573/Thermo_HW_15.mlx` | `examples/ae573_propulsion/gasdynamics_homework.py` |
| `spring 2023/AE 573/Thermo_HW_18.mlx` | `examples/ae573_propulsion/thermo_hw18.py` |
| `spring 2023/AE 573/ThermodynamicDevice.m` | Not ported: unfinished relational engine framework (superseded by `solve_relations`) |
| `spring 2023/AE 573/ThermodynamicState.m` | Not ported: unfinished relational engine framework (superseded by `solve_relations`) |
| `spring 2023/AE 573/ThermodynamicStateTester.m` | Not ported: unfinished relational engine framework (superseded by `solve_relations`) |
| `spring 2023/AE 573/debug_Relation_class.mlx` | src/unicodes/thermo/gas_law.py (`Relation`, `solve_relations`) |
| `spring 2023/AE 573/equationRepeater.mlx` | `src/unicodes/thermo/gas_law.py` |
| `spring 2023/AE 573/expansionFanProperties.mlx` | `src/unicodes/gasdynamics.py` |
| `spring 2023/AE 573/getMuAir.mlx` | src/unicodes/atmosphere.py (`sutherland_viscosity`) |
| `spring 2023/AE 573/getSatndardAtmosphericValues.mlx` | src/unicodes/atmosphere.py |
| `spring 2023/AE 573/getSatndardAtmosphericValuesStrato.mlx` | src/unicodes/atmosphere.py |
| `spring 2023/AE 573/getSatndardAtmosphericValuesTropo.mlx` | src/unicodes/atmosphere.py |
| `spring 2023/AE 573/getSpdSoundAir.mlx` | src/unicodes/atmosphere.py |
| `spring 2023/AE 573/getStandardAtmoRatiosStrato.mlx` | src/unicodes/atmosphere.py |
| `spring 2023/AE 573/getStandardAtmoRatiosTropo.mlx` | src/unicodes/atmosphere.py |
| `spring 2023/AE 573/getStandardAtmosphericValues.mlx` | src/unicodes/atmosphere.py |
| `spring 2023/AE 573/getStandardAtmosphericValuesBG.mlx` | src/unicodes/atmosphere.py |
| `spring 2023/AE 573/getTempRatioTropo.mlx` | src/unicodes/atmosphere.py |
| `spring 2023/AE 573/inletAdiabaticEfficiency.mlx` | src/unicodes/gasdynamics.py (`inlet_*`) |
| `spring 2023/AE 573/inletTotalPressureRecovery.mlx` | src/unicodes/gasdynamics.py (`inlet_*`) |
| `spring 2023/AE 573/inletTotalPressureRecoveryFromMach.mlx` | src/unicodes/gasdynamics.py (`inlet_*`) |
| `spring 2023/AE 573/isentropicFlowProperties.mlx` | `src/unicodes/gasdynamics.py` |
| `spring 2023/AE 573/mexwell_relation_experimentor.mlx` | src/unicodes/thermo/gas_law.py (`GasLaw` implicit derivatives / Maxwell relations) |
| `spring 2023/AE 573/normalShockProperties.mlx` | `src/unicodes/gasdynamics.py` |
| `spring 2023/AE 573/obliqueShockProperties.mlx` | `src/unicodes/gasdynamics.py` |
| `spring 2023/AE 573/relation_test.mlx` | src/unicodes/thermo/gas_law.py (`Relation`, `solve_relations`) |
| `spring 2023/AE 573/sigmafromThetaTropo.mlx` | src/unicodes/atmosphere.py |
| `spring 2023/AE 573/solveRelations.mlx` | `src/unicodes/thermo/gas_law.py` |
| `spring 2023/AE 573/totalPressureFromStaticAndMach.mlx` | `src/unicodes/gasdynamics.py` |

## spring 2023/EESC 316

| MATLAB | Python / note |
|---|---|
| `spring 2023/EESC 316/EECS_HW_4.mlx` | `examples/misc/other_courses.py` |

## spring 2023/MATH 526

| MATLAB | Python / note |
|---|---|
| `spring 2023/MATH 526/MedicalDataAnalysis/DMD.mlx` | `src/unicodes/decomposition/dmd.py` |
| `spring 2023/MATH 526/MedicalDataAnalysis/DataLoader.mlx` | `examples/math526_statistics/diagnosis_dmd.py` |
| `spring 2023/MATH 526/MedicalDataAnalysis/DiagnosisDMD.mlx` | `examples/math526_statistics/diagnosis_dmd.py` |
| `spring 2023/MATH 526/correlation_demonstration.mlx` | `examples/math526_statistics/correlation_demo.py` |
| `spring 2023/MATH 526/didgital_distribution_chi_squared.mlx` | `examples/math526_statistics/pi_digit_chi_squared.py` |
| `spring 2023/MATH 526/didgital_distribution_chi_squared_graphs.mlx` | `examples/math526_statistics/pi_digit_chi_squared.py` |
| `spring 2023/MATH 526/plotSingalVTime.m` | `examples/math526_statistics/correlation_demo.py` |
| `spring 2023/MATH 526/saveAllFig.mlx` | `examples/math526_statistics/correlation_demo.py` |

## spring 2023/ME 712

| MATLAB | Python / note |
|---|---|
| `spring 2023/ME 712/HW_9.mlx` | examples/misc/other_courses.py |

## spring 2024/AE 430 Instramentation

| MATLAB | Python / note |
|---|---|
| `spring 2024/AE 430 Instramentation/HW 1/Problem_1.mlx` | `examples/ae430_instrumentation/signal_analysis.py` |
| `spring 2024/AE 430 Instramentation/HW 2/HW_2.mlx` | `examples/ae430_instrumentation/signal_analysis.py` |
| `spring 2024/AE 430 Instramentation/Lab 7/Lab_7_Data_analysis_version_1.mlx` | `examples/ae430_instrumentation/signal_analysis.py` |

## spring 2024/AE 722 AIAA graduate Electric Sailplane Design

| MATLAB | Python / note |
|---|---|
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Archytas_Superimposed_Design_Constraints.mlx` | `examples/ae722_sailplane/constraint_diagram.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Blade_Element_Momentum_Theory.mlx` | `examples/ae722_sailplane/propeller_bemt.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/BreguetEndurance.mlx` | `src/unicodes/aero/performance.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/BreguetRange.mlx` | `src/unicodes/aero/performance.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Class_I_Drag_Polar.mlx` | `examples/ae722_sailplane/drag_polar_and_speed_polar.py`, `examples/ae722_sailplane/geometry.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/ClimbIndex.mlx` | `src/unicodes/aero/performance.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Climb_Coefficients.mlx` | `examples/ae722_sailplane/aspect_ratio_constraints.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Climb_Power.mlx` | `examples/ae722_sailplane/climb_power.py`, `src/unicodes/aero/performance.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Climb_Sizing.mlx` | `examples/ae722_sailplane/constraint_diagram.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Climb_Sizing_and_Drag_Polars.mlx` | examples/ae722_sailplane/constraint_diagram.py and src/unicodes/aero/performance.py (FAR 23 constraints) |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Climb_Sizing_and_Drag_Polars_EDITED.mlx` | examples/ae722_sailplane/constraint_diagram.py and src/unicodes/aero/performance.py (FAR 23 constraints) |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Climb_Sizing_and_Drag_Polars_EDITED_2.mlx` | examples/ae722_sailplane/constraint_diagram.py and src/unicodes/aero/performance.py (FAR 23 constraints) |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Combined_CrossCountry_Sizing.mlx` | `examples/ae722_sailplane/aspect_ratio_constraints.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/CrossCountryFigureFormatter.mlx` | Not ported: figure formatting only |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Daedalus_preliminary_Sizing_V_0.mlx` | AE 521 copy: examples/ae521_aircraft_design/daedalus_constraint_diagram.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/DensityRatio.mlx` | `src/unicodes/atmosphere.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Design_Line_Plot.mlx` | Not ported: figure formatting only |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/DragCoefficent_RateOfClimb_Max.mlx` | `src/unicodes/aero/performance.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/DragCoefficientZeroLift.mlx` | `src/unicodes/aero/performance.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/DragCoefficient_Airbrakes.mlx` | `src/unicodes/aero/performance.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/DragPolar_parabolic.mlx` | `src/unicodes/aero/wing.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Electrical Systems/PowerDiagram.mlx` | `examples/ae722_sailplane/power_budget.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Empenage_Geometery.mlx` | `examples/ae722_sailplane/geometry.py`, `src/unicodes/aero/wing.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/EquivalentParasiteArea.mlx` | `src/unicodes/aero/performance.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/FuelFraction_Endurance.mlx` | `src/unicodes/aero/performance.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/FuelFraction_Range.mlx` | `src/unicodes/aero/performance.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/FuelFraction_Suggested.mlx` | Reads Roskam's CSV; the table values are in examples/ae521_aircraft_design/daedalus_weight_sizing.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/HilghLift_Sizing.mlx` | `examples/ae722_sailplane/climb_power.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/HorstmannThermal.mlx` | `src/unicodes/aero/soaring.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Internal_Combustion_Engine_Predictor_Equation.mlx` | `src/unicodes/aero/performance.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Landing_Requirements_Reida.mlx` | `examples/ae722_sailplane/constraint_diagram.py`, `src/unicodes/aero/performance.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Landing_Requirements_Reida_edited.mlx` | `examples/ae722_sailplane/constraint_diagram.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/LiftCoefficent_RateOfClimb_Max.mlx` | `src/unicodes/aero/performance.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Lift_to_Drag_Requirements.mlx` | `examples/ae722_sailplane/aspect_ratio_constraints.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Lift_to_Drag_Requirements_V2.mlx` | examples/ae722_sailplane/aspect_ratio_constraints.py (`w2s_for_lift_to_drag`) |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Lift_to_drag_algebra.mlx` | examples/ae722_sailplane/aspect_ratio_constraints.py (`w2s_for_lift_to_drag`) |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Lift_to_drag_constraint_animated_plot.m` | examples/ae722_sailplane/aspect_ratio_constraints.py (`w2s_for_lift_to_drag`) |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/LoadFactor_Gust_CS22.mlx` | `src/unicodes/aero/loads.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/LoadFactor_ManuverLimit_FAR23.mlx` | `src/unicodes/aero/loads.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/LoadFactor_Manuver_CS22.mlx` | `src/unicodes/aero/loads.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/LoadFactor_Manuver_FAR23.mlx` | `src/unicodes/aero/loads.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/MeanGeometricCord.mlx` | `src/unicodes/aero/wing.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Nonstandardrho.mlx` | `src/unicodes/atmosphere.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/AETHER_Ben_6DOF_with_optimizer.mlx` | `examples/ae722_sailplane/pah/pah_optimizer.py`, `src/unicodes/flight_dynamics/sixdof.py`, `src/unicodes/optimize/genetic.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/AETHER_Ben_6DOF_with_optimizer_AAA_V0.mlx` | examples/ae722_sailplane/pah/pah_optimizer.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/AETHER_Ben_6DOF_with_optimizer_AAA_V1.mlx` | examples/ae722_sailplane/pah/pah_optimizer.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/AETHER_Ben_6DOF_with_optimizer_AAA_V1_0.mlx` | examples/ae722_sailplane/pah/pah_optimizer.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/AETHER_Ben_6DOF_with_optimizer_AAA_V1_1.mlx` | examples/ae722_sailplane/pah/pah_optimizer.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/AETHER_Ben_6DOF_with_optimizer_AAA_V1_2.mlx` | examples/ae722_sailplane/pah/pah_optimizer.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/AETHER_Ben_6DOF_with_optimizer_AAA_V1_3.mlx` | examples/ae722_sailplane/pah/pah_optimizer.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/AETHER_Ben_6DOF_with_optimizer_AAA_V1_4.mlx` | examples/ae722_sailplane/pah/pah_optimizer.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/AETHER_Ben_6DOF_with_optimizer_AAA_V1_5.mlx` | examples/ae722_sailplane/pah/pah_optimizer.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/AETHER_Ben_6DOF_with_optimizer_AAA_V1_6.mlx` | examples/ae722_sailplane/pah/pah_optimizer.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/AETHER_Ben_6DOF_with_optimizer_AAA_V1_7.mlx` | examples/ae722_sailplane/pah/pah_optimizer.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/AETHER_Ben_6DOF_with_optimizer_AAA_V2.mlx` | examples/ae722_sailplane/pah/pah_optimizer.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/Adaptive_Aerocompliant_Airfoil_Dynamics.mlx` | `src/unicodes/flight_dynamics/sixdof.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/Aether_6DOF.mlx` | `examples/ae722_sailplane/pah/pah_optimizer.py`, `src/unicodes/flight_dynamics/sixdof.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/Airfoil.m` | src/unicodes/aero/airfoil.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/BenMays6DOF.mlx` | `examples/ae722_sailplane/pah/sixdof_doublet.py`, `src/unicodes/flight_dynamics/sixdof.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/BenMays6DOF_PAH.mlx` | `examples/ae722_sailplane/pah/sixdof_doublet.py`, `src/unicodes/flight_dynamics/sixdof.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/BenMays6DOF_PAH_cessna182_2Flap.mlx` | `examples/ae722_sailplane/pah/sixdof_doublet.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/BenMays6DOF_PAH_sailplane.mlx` | `examples/ae722_sailplane/pah/sixdof_doublet.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/BenMays6DOF_PAH_sailplane_1Flap.mlx` | `examples/ae722_sailplane/pah/sixdof_doublet.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/BenMays6DOF_PAH_sailplane_2Flap.mlx` | `examples/ae722_sailplane/pah/sixdof_doublet.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/BenMays_6DOF_jacobian.mlx` | `examples/ae722_sailplane/pah/pah_modes.py`, `src/unicodes/flight_dynamics/sixdof.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/Ben_6DOF_with_optimizer.mlx` | `examples/ae722_sailplane/pah/pah_optimizer.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/Ben_6DOF_with_optimizer_fast_script.m` | `examples/ae722_sailplane/pah/pah_optimizer.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/Ben_6DOF_with_optimizer_fast_script_SA.m` | examples/ae722_sailplane/pah/pah_optimizer.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/Cessna_182_PAH.mlx` | `examples/ae722_sailplane/pah/sixdof_doublet.py`, `src/unicodes/flight_dynamics/sixdof.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/Getting_Interpolated_Mesh_Sheet_Intersections.mlx` | src/unicodes/aero/pah.py and examples/ae722_sailplane/pah/pah_surfaces.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/HIB.mlx` | src/unicodes/flight_dynamics/kinematics.py (`dcm_inertial_to_body`) |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/Lift_mesh_Xfoil.mlx` | `src/unicodes/aero/xfoil.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/Miscelanious_scitetech_2025_figures.mlx` | Not ported: figure formatting only |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/Miscelanious_scitetech_2025_figures_2.mlx` | Not ported: figure formatting only |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/PAH_Lift_Surface_Gradients.mlx` | `examples/ae722_sailplane/pah/pah_surfaces.py`, `src/unicodes/aero/pah.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/PAH_Moment_Surface_Gradients.mlx` | `examples/ae722_sailplane/pah/pah_surfaces.py`, `src/unicodes/aero/pah.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/PAH_generation_ploter.mlx` | `examples/ae722_sailplane/pah/pah_surfaces.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/PAH_generator.mlx` | `examples/ae722_sailplane/pah/pah_surfaces.py`, `src/unicodes/aero/airfoil.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/PitchMomentCoefficient_withRespectTo_Alpha_horizontalTail.mlx` | src/unicodes/aero/wing.py (`cm_alpha_*`, `cm0_tail`, `cm_tail_incidence`) |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/PitchMomentCoefficient_withRespectTo_Alpha_wingFueslage.mlx` | src/unicodes/aero/wing.py (`cm_alpha_*`, `cm0_tail`, `cm_tail_incidence`) |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/PitchMomentCoefficient_withRespectTo_horizontalTailIncidence.mlx` | src/unicodes/aero/wing.py (`cm_alpha_*`, `cm0_tail`, `cm_tail_incidence`) |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/PitchMomentCoefficient_zeroAoA_horizontalTail.mlx` | src/unicodes/aero/wing.py (`cm_alpha_*`, `cm0_tail`, `cm_tail_incidence`) |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/RK4.mlx` | `src/unicodes/ode.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/ReadAirfoil.mlx` | `examples/ae722_sailplane/pah/pah_surfaces.py`, `src/unicodes/aero/airfoil.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/Roskam_6DOF_C182.mlx` | `examples/ae722_sailplane/pah/sixdof_doublet.py`, `src/unicodes/flight_dynamics/sixdof.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/Rotation_Derivation_for_Flaps.mlx` | Symbolic derivation; the resulting equations are in src/unicodes/flight_dynamics/sixdof.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/Vehicle_dynamics_derivation.mlx` | Symbolic derivation; the resulting equations are in src/unicodes/flight_dynamics/sixdof.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/adaptive_airfoil_dynamic_simulator.mlx` | `examples/ae722_sailplane/pah/flap_phase_portrait.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/adaptive_airfoil_dynamics_linearizedEvaluation.mlx` | `examples/ae722_sailplane/pah/flap_phase_portrait.py`, `src/unicodes/flight_dynamics/sixdof.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/airfoil_file_inspector.mlx` | `examples/ae722_sailplane/pah/pah_surfaces.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/assess_comfort.mlx` | `src/unicodes/flight_dynamics/ride_quality.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/calculate_crest_factors.mlx` | src/unicodes/flight_dynamics/ride_quality.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/calculate_rms.mlx` | src/unicodes/flight_dynamics/ride_quality.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/calculate_vdv.mlx` | `src/unicodes/flight_dynamics/ride_quality.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/computeJacobian.mlx` | `src/unicodes/numerics.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/design_iso2631_filter.mlx` | `src/unicodes/flight_dynamics/ride_quality.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/display_ISO2631_results.mlx` | `examples/ae722_sailplane/pah/ride_quality_comparison.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/downwashGradient.mlx` | `src/unicodes/aero/wing.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/eta_h_Function_Example.mlx` | `examples/ae722_sailplane/tail_blanking.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/generate_dryden_gust.mlx` | `src/unicodes/flight_dynamics/turbulence.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/generate_dryden_gust_with_rotational.mlx` | `src/unicodes/flight_dynamics/turbulence.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/getAirFlowAngles.mlx` | src/unicodes/flight_dynamics/kinematics.py (`air_flow_angles`) |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/getFlightPathAngles.mlx` | src/unicodes/flight_dynamics/kinematics.py (`flight_path_angles`) |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/hardening_spring_tester.mlx` | `examples/ae722_sailplane/pah/flap_phase_portrait.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/load_example_data.mlx` | `examples/ae722_sailplane/pah/ride_quality_comparison.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/plot_ISO2631_results.mlx` | `examples/ae722_sailplane/pah/ride_quality_comparison.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/polhamus.mlx` | `src/unicodes/aero/wing.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/random_smoothed_number_playground.m` | `src/unicodes/numerics.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/ride_quality_analsyis.mlx` | `examples/ae722_sailplane/pah/ride_quality_comparison.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/ride_quality_analsyis_V1.m` | `examples/ae722_sailplane/pah/ride_quality_comparison.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/ride_quality_analysis.mlx` | `examples/ae722_sailplane/pah/ride_quality_comparison.py`, `src/unicodes/flight_dynamics/ride_quality.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/xfoil.m` | `src/unicodes/aero/xfoil.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/PAH airfoils/xfoilCl.m` | `src/unicodes/aero/xfoil.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Power_Lapse_With_Alt_IC_engine.mlx` | `src/unicodes/aero/performance.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Propeller Optimization/CE_PrOptimization.m` | `examples/ae722_sailplane/propeller_optimization.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Propeller Optimization/ClimbIndex.mlx` | `src/unicodes/aero/performance.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Propeller Optimization/ClimbPowerCalculator.mlx` | `examples/ae722_sailplane/climb_power.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Propeller Optimization/LiftCoefficent_RateOfClimb_Max.mlx` | `src/unicodes/aero/performance.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Propeller Optimization/Nonstandardrho.mlx` | `src/unicodes/atmosphere.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Propeller Optimization/P_Climb_Prop.mlx` | examples/ae722_sailplane/climb_power.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Propeller Optimization/Plot_PrOptiamal_Solution.mlx` | `examples/ae722_sailplane/propeller_optimization.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Propeller Optimization/PropDesign_SetPower.m` | `examples/ae722_sailplane/propeller_optimization.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Propeller Optimization/V_sink.mlx` | `examples/ae722_sailplane/drag_polar_and_speed_polar.py`, `src/unicodes/aero/soaring.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Propeller Optimization/W2P_Climb_Propeller.mlx` | `src/unicodes/aero/performance.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Propeller Optimization/speedpolarplot.mlx` | `examples/ae722_sailplane/drag_polar_and_speed_polar.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Report_3_voyager.mlx` | AE 521 copy: examples/ae521_aircraft_design/voyager_power.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Report_4_DC_3.mlx` | AE 521 copy: examples/ae521_aircraft_design/dc3_schrenk.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Required_Power_Subsonic.mlx` | `src/unicodes/aero/performance.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Ride_quality.mlx` | `examples/ae722_sailplane/ride_quality_index.py`, `src/unicodes/aero/soaring.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/STAMPED_empenage.mlx` | src/unicodes/aero/sizing.py (`stamped_trend`) and examples/ae722_sailplane/geometry.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Schrenk_Lift_aproximation.mlx` | `src/unicodes/aero/wing.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Service_Ceiling_Requirements_Reida.mlx` | `examples/ae722_sailplane/constraint_diagram.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Service_Ceiling_Requirements_Reida_edited.mlx` | `examples/ae722_sailplane/constraint_diagram.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/SinkRateCircle.mlx` | examples/ae722_sailplane/cross_country_speed.py and src/unicodes/aero/soaring.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/SinkRate_Requirment.mlx` | `examples/ae722_sailplane/aspect_ratio_constraints.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/SinkRate_Requirment_algebra.mlx` | examples/ae722_sailplane/cross_country_speed.py and src/unicodes/aero/soaring.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Sizing.m` | `examples/ae722_sailplane/weight_sizing_electric.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Sizing_FORMATTED.m` | examples/ae722_sailplane/constraint_diagram.py and src/unicodes/aero/performance.py (FAR 23 constraints) |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Sizing_edited.m` | examples/ae722_sailplane/constraint_diagram.py and src/unicodes/aero/performance.py (FAR 23 constraints) |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/SpeedDive.mlx` | `src/unicodes/aero/performance.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Speed_Polar.mlx` | `examples/ae722_sailplane/drag_polar_and_speed_polar.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Speed_Polars.mlx` | `examples/ae722_sailplane/drag_polar_and_speed_polar.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/StallSpeed.mlx` | `src/unicodes/aero/wing.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/T2W_ClimbGradient_FAR25.mlx` | `src/unicodes/aero/performance.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/T2W_Climb_FAR25.mlx` | `src/unicodes/aero/performance.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/T2W_CruiseSpeed.mlx` | `src/unicodes/aero/performance.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/T2W_Takeoff_FAR25.mlx` | `src/unicodes/aero/performance.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/T2W_TimeToClimb.mlx` | `src/unicodes/aero/performance.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Takeoff_Requirements_Reworked.mlx` | `examples/ae722_sailplane/constraint_diagram.py`, `src/unicodes/aero/performance.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Thermal Modeling/Algebraic_solution_to_sink_rate_constraint.mlx` | examples/ae722_sailplane/cross_country_speed.py and src/unicodes/aero/soaring.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Thermal Modeling/AverageCrossCountrySpeed.mlx` | examples/ae722_sailplane/cross_country_speed.py and src/unicodes/aero/soaring.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Thermal Modeling/AverageCrossCountrySpeed_Horstmnan.mlx` | `src/unicodes/aero/soaring.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Thermal Modeling/AverageCrossCountrySpeed_Quast.mlx` | `src/unicodes/aero/soaring.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Thermal Modeling/CLoptimum.m` | examples/ae722_sailplane/cross_country_speed.py and src/unicodes/aero/soaring.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Thermal Modeling/HorstmannThermal.mlx` | `src/unicodes/aero/soaring.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Thermal Modeling/LiftCoefficient_Optimal_Interthermal_Glide.mlx` | `src/unicodes/aero/soaring.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Thermal Modeling/MeanGeometricCord.mlx` | `src/unicodes/aero/wing.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Thermal Modeling/Nonstandardrho.mlx` | `src/unicodes/atmosphere.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Thermal Modeling/OptimumCL.mlx` | `examples/ae722_sailplane/cross_country_speed.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Thermal Modeling/SinkRateCircle.mlx` | examples/ae722_sailplane/cross_country_speed.py and src/unicodes/aero/soaring.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Thermal Modeling/Themral_Modeling_Part_2.mlx` | `examples/ae722_sailplane/cross_country_speed.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Thermal Modeling/Themral_Modeling_and_CL_opt_Poz.mlx` | examples/ae722_sailplane/cross_country_speed.py and src/unicodes/aero/soaring.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Thermal Modeling/Thermal_Modeling.mlx` | `examples/ae722_sailplane/cross_country_speed.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Thermal Modeling/Thermal_Modeling_V1.mlx` | examples/ae722_sailplane/cross_country_speed.py and src/unicodes/aero/soaring.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Thermal Modeling/V_sc.mlx` | `src/unicodes/aero/soaring.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Thermal Modeling/Vavgfsd.m` | examples/ae722_sailplane/cross_country_speed.py and src/unicodes/aero/soaring.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Thermal Modeling/Vclimbspeed.m` | `src/unicodes/aero/soaring.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Thermal Modeling/Vsinkrate.m` | examples/ae722_sailplane/cross_country_speed.py and src/unicodes/aero/soaring.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Thermal Modeling/VsinkrateW2S.m` | examples/ae722_sailplane/cross_country_speed.py and src/unicodes/aero/soaring.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Thermal Modeling/Vsinkrate_phi.m` | examples/ae722_sailplane/cross_country_speed.py and src/unicodes/aero/soaring.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Thermal Modeling/Vsinkratephi.m` | examples/ae722_sailplane/cross_country_speed.py and src/unicodes/aero/soaring.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Thermal Modeling/Vsinkratevk.m` | examples/ae722_sailplane/cross_country_speed.py and src/unicodes/aero/soaring.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Thermal Modeling/Vthermal.m` | `src/unicodes/aero/soaring.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Thermal Modeling/untitled.m` | `examples/ae722_sailplane/cross_country_speed.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Thermal_Modeling.mlx` | `examples/ae722_sailplane/cross_country_speed.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/V_sc.mlx` | `src/unicodes/aero/soaring.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/V_sink.mlx` | `examples/ae722_sailplane/drag_polar_and_speed_polar.py`, `src/unicodes/aero/soaring.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Vclimbspeed.m` | `src/unicodes/aero/soaring.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Vn_Diagram_V_0.mlx` | `examples/ae722_sailplane/vn_diagram.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Vn_Diagram_V_1.mlx` | `examples/ae722_sailplane/vn_diagram.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Vn_Diagram_V_1_modified_forFAR23.mlx` | examples/ae722_sailplane/vn_diagram.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Vsinkrate.m` | examples/ae722_sailplane/cross_country_speed.py and src/unicodes/aero/soaring.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Vsinkrate_phi.m` | examples/ae722_sailplane/cross_country_speed.py and src/unicodes/aero/soaring.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Vthermal.m` | `src/unicodes/aero/soaring.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/W2P_Climb_Propeller.mlx` | `src/unicodes/aero/performance.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/W2S_Landing_FAR25.mlx` | `src/unicodes/aero/performance.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/W2S_LiftToDragAtSpeed.mlx` | `src/unicodes/aero/performance.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/W2S_SinkRateCircle.mlx` | `src/unicodes/aero/soaring.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/W2S_TimeToClimb.mlx` | `src/unicodes/aero/performance.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Weight Sizing/Archytas_Weight_Sizing.mlx` | `examples/ae722_sailplane/weight_sizing_electric.py`, `src/unicodes/aero/sizing.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Weight Sizing/LD.mlx` | `src/unicodes/aero/wing.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Weight Sizing/Nonstandardrho.mlx` | `src/unicodes/atmosphere.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Weight Sizing/battPED.mlx` | `examples/ae722_sailplane/weight_sizing_electric.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Weight Sizing/powerest.mlx` | `src/unicodes/aero/performance.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Weight Sizing/tottime.mlx` | examples/ae722_sailplane/weight_sizing_electric.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Weight_Crew.mlx` | `src/unicodes/aero/performance.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Weight_Empty.mlx` | `src/unicodes/aero/performance.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Weight_PassengerAndBaggage.mlx` | `src/unicodes/aero/performance.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Weight_Sizing_Daedalus_V_0.mlx` | AE 521 copy: examples/ae521_aircraft_design/daedalus_weight_sizing.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/WettedArea.mlx` | `src/unicodes/aero/performance.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/WettedArea_Fuselage_Class_I.mlx` | `src/unicodes/aero/performance.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/WettedArea_Planform_Class_I.mlx` | `src/unicodes/aero/performance.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/Windtunnel_Blanking_Analysis_V0.mlx` | `examples/ae722_sailplane/tail_blanking.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/cross-country_modeling_PAH.mlx` | examples/ae722_sailplane/cross_country_speed.py and src/unicodes/aero/soaring.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/polhamus.mlx` | `src/unicodes/aero/wing.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/prelim_code_to_to_run.m` | AE 521 copy: examples/ae521_aircraft_design/daedalus_constraint_diagram.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/speedpolarplot.mlx` | `examples/ae722_sailplane/drag_polar_and_speed_polar.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/vehicle_weight_sizing_1.m` | AE 521 copy: examples/ae521_aircraft_design/hermes_weight_sizing.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/weight_Sizing_Hermes_V_0.mlx` | AE 521 copy: examples/ae521_aircraft_design/hermes_weight_sizing.py |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/wing_geometery.mlx` | `examples/ae722_sailplane/geometry.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/wing_geometery2.mlx` | `examples/ae722_sailplane/geometry.py`, `src/unicodes/aero/wing.py` |
| `spring 2024/AE 722 AIAA graduate Electric Sailplane Design/wing_geometry_algebra.mlx` | examples/ae722_sailplane/geometry.py (`tapered_wing`) |

## unsorted

| MATLAB | Python / note |
|---|---|
| `unsorted/DMDTurbulentSeperationBubble.mlx` | `examples/data_driven/separation_bubble.py` |
| `unsorted/DMD_benchmarking_againgst_simple_functions.mlx` | `examples/data_driven/dmd_simple_functions.py` |
| `unsorted/DMD_combustor.mlx` | examples/data_driven/combustor_dmd.py |
| `unsorted/DMD_combustor_fast.m` | examples/data_driven/combustor_dmd.py |
| `unsorted/DMD_simple_functions.mlx` | `examples/data_driven/dmd_simple_functions.py` |
| `unsorted/DynamicModeDecomposer.m` | `examples/data_driven/dmd_simple_functions.py` |
| `unsorted/InOutFormat.m` | src/unicodes/io/tecplot.py |
| `unsorted/JacksonRocketPropullsionCompiledCode771.m` | `examples/jackson_torok/rocket_propulsion_771.py` |
| `unsorted/Lab_1.m.mlx` | Copy of AE 546 Lab 1: examples/ae546_aero_lab/lab1_pressure_distribution.py |
| `unsorted/Logistic.mlx` | `examples/misc/small_scripts.py` |
| `unsorted/OutputTecASCIIdata.m` | `src/unicodes/io/tecplot.py` |
| `unsorted/PODwavy.m` | `examples/data_driven/pod_video.py` |
| `unsorted/PODwavy.mlx` | `examples/data_driven/pod_video.py` |
| `unsorted/ProperOrthogonalDecomposer.m` | src/unicodes/decomposition/pod.py |
| `unsorted/ProperOrthogonalDecomposerVersion2.m` | src/unicodes/decomposition/pod.py |
| `unsorted/Proper_orthogonal_decomposition_Grey_Ignition_video.mlx` | `examples/data_driven/pod_video.py` |
| `unsorted/Proper_orthogonal_decomposition_RGB_Ignition.mlx` | `examples/data_driven/pod_video.py` |
| `unsorted/Proper_orthogonal_decomposition_turbulent_seperation_bubble.mlx` | `examples/data_driven/separation_bubble.py` |
| `unsorted/RealTimeFFTtst.m` | Not ported: live FFT animation demo (see src/unicodes/decomposition/spectral.py `fft_amplitude`) |
| `unsorted/SaveVideoToArray.mlx` | `examples/data_driven/pod_video.py` |
| `unsorted/Spectral_POD.mlx` | Not ported: empty file |
| `unsorted/blackBox2.m` | Lorenz test data: src/unicodes/decomposition/sindy.py (`lorenz`) |
| `unsorted/blackbox.m` | Pi-digit-sum test data: examples/misc/small_scripts.py |
| `unsorted/function_comparison_timing.m` | Not ported: MATLAB timing / parpool experiments |
| `unsorted/generate_from_blackbox.m` | Pi-digit-sum test data: examples/misc/small_scripts.py |
| `unsorted/getBlackboxData2D.m` | Lorenz test data: src/unicodes/decomposition/sindy.py (`lorenz`) |
| `unsorted/getBlackboxData3D.m` | Lorenz test data: src/unicodes/decomposition/sindy.py (`lorenz`) |
| `unsorted/getsnapshotmatrix.m` | src/unicodes/decomposition/snapshots.py |
| `unsorted/imageVisualizationTest.mlx` | `examples/data_driven/pod_video.py` |
| `unsorted/importGridFile.m` | `src/unicodes/io/tecplot.py` |
| `unsorted/importTecASCIIdata.m` | `src/unicodes/io/tecplot.py` |
| `unsorted/lorenz.mlx` | `examples/sindy_examples/ex02_lorenz.py`, `src/unicodes/decomposition/sindy.py` |
| `unsorted/messing_around_with_algebra.m` | Not ported: symbolic scratch work |
| `unsorted/parallel_test.m` | Not ported: MATLAB timing / parpool experiments |
| `unsorted/removemean.m` | src/unicodes/decomposition/snapshots.py |
| `unsorted/restoremean.mlx` | src/unicodes/decomposition/snapshots.py |
| `unsorted/rsvd.m` | src/unicodes/decomposition/linalg.py (`rsvd`) |
| `unsorted/sampleIO.m` | src/unicodes/io/tecplot.py |
| `unsorted/sampleO_comented.m` | src/unicodes/io/tecplot.py |
| `unsorted/sampleO_parallelized.m` | src/unicodes/io/tecplot.py |
| `unsorted/simLorenz.mlx` | `examples/data_driven/sindy_examples.py` |
| `unsorted/sumDidgetsOfPi.mlx` | `examples/misc/small_scripts.py` |
| `unsorted/svdms.mlx` | Not ported: unfinished (does not run) |
| `unsorted/test_parallelizing_for.m` | Not ported: MATLAB timing / parpool experiments |
| `unsorted/timeandSpacialWave.mlx` | `examples/data_driven/dmd_simple_functions.py` |
| `unsorted/timedynamics_experimentation.mlx` | `examples/data_driven/dmd_simple_functions.py` |
| `unsorted/trackGearCalculator.mlx` | `examples/misc/small_scripts.py` |
| `unsorted/video2array.mlx` | `examples/data_driven/pod_video.py` |
| `unsorted/weather_probability_advancing.mlx` | `examples/misc/small_scripts.py` |
