###########################################################################
#
# # Numerical Astrodynamics 2025/2026
#
# # Assignment 1 - Propagation Settings
#
###########################################################################


""" 
Copyright (c) 2010-2020, Delft University of Technology
All rights reserved

This file is part of Tudat. Redistribution and use in source and 
binary forms, with or without modification, are permitted exclusively
under the terms of the Modified BSD license. You should have received
a copy of the license with this file. If not, please or visit:
http://tudat.tudelft.nl/LICENSE.
"""

import numpy as np

from tudatpy import constants
from tudatpy.interface import spice
from tudatpy.dynamics import environment_setup, propagation_setup, simulator
from tudatpy.data import save2txt

import os


def run_juice_simulation(
    simulation_start_epoch, 
    simulation_end_epoch, 
    force_coeff, 
    radiation_pressure_coefficient,
    current_directory
):
    # Load spice kernels
    spice.load_standard_kernels()
    spice.load_kernel(current_directory + "/juice_mat_crema_5_1_150lb_v01.bsp")

    # Create settings for celestial bodies
    bodies_to_create = ['Ganymede', 'Jupiter', 'Sun', 'Saturn', 'Europa', 'Io', 'Callisto']
    global_frame_origin = "Ganymede"
    global_frame_orientation = "ECLIPJ2000"
    body_settings = environment_setup.get_default_body_settings(
        bodies_to_create, global_frame_origin, global_frame_orientation
    )

    # Atmosphere model
    body_settings.get('Ganymede').atmosphere_settings = environment_setup.atmosphere.exponential(
        scale_height=40000.0, surface_density=2e-9)

    # Sun radiation
    irradiance_at_1AU = 1367.0
    luminosity_model_settings = environment_setup.radiation_pressure.irradiance_based_constant_luminosity(
        irradiance_at_1AU, constants.ASTRONOMICAL_UNIT)
    radiation_source_settings_sun = environment_setup.radiation_pressure.isotropic_radiation_source(
        luminosity_model_settings)
    body_settings.get("Sun").radiation_source_settings = radiation_source_settings_sun

    eclipse_settings = {'Sun': ['Ganymede']}

    # Create Vehicle
    body_settings.add_empty_settings("JUICE")
    body_settings.get('JUICE').aerodynamic_coefficient_settings = environment_setup.aerodynamic_coefficients.constant(
        reference_area=100.0, constant_force_coefficient=force_coeff)
    body_settings.get('JUICE').radiation_pressure_target_settings = environment_setup.radiation_pressure.cannonball_radiation_target(
        reference_area=100, radiation_pressure_coefficient=radiation_pressure_coefficient, 
        per_source_occulting_bodies=eclipse_settings)

    bodies = environment_setup.create_system_of_bodies(body_settings)
    bodies.get('JUICE').mass = 2e3

    # Accelerations
    acceleration_settings_on_vehicle = dict(
        Ganymede=[propagation_setup.acceleration.spherical_harmonic_gravity(2,2),
                  propagation_setup.acceleration.aerodynamic()],
        Jupiter=[propagation_setup.acceleration.spherical_harmonic_gravity(4,0)],
        Sun=[propagation_setup.acceleration.point_mass_gravity(),
             propagation_setup.acceleration.radiation_pressure()],
        Saturn=[propagation_setup.acceleration.point_mass_gravity()],
        Europa=[propagation_setup.acceleration.point_mass_gravity()],
        Io=[propagation_setup.acceleration.point_mass_gravity()],
        Callisto=[propagation_setup.acceleration.point_mass_gravity()]
    )

    acceleration_models = propagation_setup.create_acceleration_models(
        bodies, {"JUICE": acceleration_settings_on_vehicle}, ["JUICE"], ["Ganymede"]
    )

    # Initial state
    system_initial_state = spice.get_body_cartesian_state_at_epoch(
        "JUICE", "Ganymede", "ECLIPJ2000", "NONE", simulation_start_epoch)

    # Output variables
    dependent_variables_to_save = [
        propagation_setup.dependent_variable.keplerian_state("JUICE", "Ganymede"),  #  Keplerian elements, 0-5
        propagation_setup.dependent_variable.single_acceleration(                   
            propagation_setup.acceleration.aerodynamic_type, 'JUICE', 'Ganymede'), # Aerodynamic acceleration, 6-8
        propagation_setup.dependent_variable.single_acceleration(
            propagation_setup.acceleration.radiation_pressure_type, 'JUICE', 'Sun'), # SRP acceleration, 9-11
        propagation_setup.dependent_variable.relative_position("JUICE", "Sun"), # Relative position to the Sun, 12-14
        propagation_setup.dependent_variable.intermediate_aerodynamic_rotation_matrix_variable( 
            'JUICE', environment_setup.aerodynamic_coefficients.AerodynamicsReferenceFrames.aerodynamic_frame, 
            environment_setup.aerodynamic_coefficients.AerodynamicsReferenceFrames.inertial_frame, 'Ganymede'), # Rotation matrix from aerodynamic to inertial frame, 15-23
        propagation_setup.dependent_variable.airspeed("JUICE", "Ganymede"), # Airspeed, 24
        propagation_setup.dependent_variable.density("JUICE", "Ganymede"), # Air density, 25
        propagation_setup.dependent_variable.rsw_to_inertial_rotation_matrix("JUICE", "Ganymede") # Rotation matrix from RSW to nertial frame, 26-34
    ]

    integrator_settings = propagation_setup.integrator.runge_kutta_fixed_step(
        10.0, coefficient_set=propagation_setup.integrator.CoefficientSets.rk_4)

    propagator_settings = propagation_setup.propagator.translational(
        ["Ganymede"], acceleration_models, ["JUICE"], system_initial_state,
        simulation_start_epoch, integrator_settings, 
        propagation_setup.propagator.time_termination(simulation_end_epoch),
        output_variables=dependent_variables_to_save)

    dynamics_simulator = simulator.create_dynamics_simulator(bodies, propagator_settings)
    
    return dynamics_simulator.propagation_results


current_task = 'question3'



# Retrieve current directory
current_directory = os.getcwd()

# # student number: 1244779 --> 1244ABC
A = 1
B = 5
C = 1

simulation_start_epoch = (
    35.4 * constants.JULIAN_YEAR
    + A * 7.0 * constants.JULIAN_DAY
    + B * constants.JULIAN_DAY
    + C * constants.JULIAN_DAY / 24.0
)
simulation_end_epoch = simulation_start_epoch + 344.0 * constants.JULIAN_DAY / 24.0

simulation_cases = ['complete', 'no_drag', 'no_srp']
all_case_results = {}

for case in simulation_cases:
    print(f"Running simulation case: {case}")
    
    # Toggle coefficients based on case
    force_coeff = [1.2, 0.0, 0.0] if case != 'no_drag' else [0.0, 0.0, 0.0]
    srp_coeff = 1.2 if case != 'no_srp' else 0.0
    
    # Run simulation
    results = run_juice_simulation(
        simulation_start_epoch, 
        simulation_end_epoch, 
        force_coeff, 
        srp_coeff, 
        current_directory
    )
    
    all_case_results[case] = results

    # Save results specific to each case
    save2txt(results.state_history,
              f"JUICEPropagationHistory_Q3_{case}.txt",
                f"./Data/{current_task}")
    
    save2txt(results.dependent_variable_history,
              f"JUICEPropagationHistory_DependentVariables_Q3_{case}.txt",
                f"./Data/{current_task}")


time = results.state_history.keys()
time_days = [
    t / constants.JULIAN_DAY - simulation_start_epoch / constants.JULIAN_DAY
    for t in time
]

np.savetxt(f'./Data/{current_task}/epochs_days.txt', time_days)



####################################################################################################
# Post-Processing 
####################################################################################################
complete_results = all_case_results['complete']
no_drag_results = all_case_results['no_drag']   
no_srp_results = all_case_results['no_srp']

complete_dependent_variables = np.vstack(list(complete_results.dependent_variable_history.values()))
no_drag_dependent_variables = np.vstack(list(no_drag_results.dependent_variable_history.values()))
no_srp_dependent_variables = np.vstack(list(no_srp_results.dependent_variable_history.values()))



############################################################################################################################################
#                                                               TASK (a)                                                                   # 
############################################################################################################################################
# Aerodynamic acceleration analytical computation
R_aer_inertial = complete_dependent_variables[:, 15:24]
R_aer_inertial = R_aer_inertial.reshape(R_aer_inertial.shape[0], 3, 3)
air_speed = complete_dependent_variables[:, 24]
air_density = complete_dependent_variables[:, 25]

drag_body = 1/2 * air_density[:, None] * air_speed[:, None]**2 * 100 * np.array([1.2, 0, 0])

acc_aer_analyt = - 1/(2e3) * np.einsum('ijk,ik->ij', R_aer_inertial, drag_body)

np.savetxt(f'./Data/{current_task}/acc_aer_analyt.txt', acc_aer_analyt)


# SRP acceleration analytical computation
r_sun_juice = complete_dependent_variables[:, 12:15]
r_sun_juice_norm = np.linalg.norm(r_sun_juice, axis=1)
power = 1367.0 * 4 * np.pi * constants.ASTRONOMICAL_UNIT**2

acc_rad_sun_analyt = power * 100 * 1.2 / (4 * np.pi * constants.SPEED_OF_LIGHT * 2e3) * r_sun_juice / r_sun_juice_norm[:, None]**3               

np.savetxt(f'./Data/{current_task}/acc_rad_sun_analyt.txt', acc_rad_sun_analyt)




############################################################################################################################################
#                                                               TASK (b)                                                                   # 
############################################################################################################################################
from tudatpy.astro import frame_conversion

R_rsw_inertial = complete_dependent_variables[:, 26:35]
R_rsw_inertial = R_rsw_inertial.reshape(R_rsw_inertial.shape[0], 3, 3)
R_inertial_rsw = R_rsw_inertial.transpose(0, 2, 1)

aer_acc_num = complete_dependent_variables[:, 6:9]
acc_aer_rsw = np.einsum('ijk,ik->ij', R_inertial_rsw, aer_acc_num)  

srp_acc_num = complete_dependent_variables[:, 9:12]
acc_srp_rsw = np.einsum('ijk,ik->ij', R_inertial_rsw, srp_acc_num)  


np.savetxt(f'./Data/{current_task}/acc_aer_rsw.txt', acc_aer_rsw)
np.savetxt(f'./Data/{current_task}/acc_srp_rsw.txt', acc_srp_rsw)



############################################################################################################################################
#                                                               TASK (c)                                                                   # 
############################################################################################################################################
np.savetxt(f'./Data/{current_task}/complete_cartesian_state.txt', np.vstack(list(complete_results.state_history.values())))
np.savetxt(f'./Data/{current_task}/no_drag_cartesian_state.txt', np.vstack(list(no_drag_results.state_history.values())))
np.savetxt(f'./Data/{current_task}/no_srp_cartesian_state.txt', np.vstack(list(no_srp_results.state_history.values())))







############################################################################################################################################
#                                                            DATA SAVE                                                                     # 
############################################################################################################################################
# Saving results for final state comparison 
final_car_state_no_drag = np.vstack(list(no_drag_results.state_history.values()))[-1, :]
final_time = np.array(list(no_drag_results.state_history.keys()))[-1]
cartesian_states_results = []; cartesian_states_results.append(np.hstack([final_time, final_car_state_no_drag]))

final_car_state_no_srp = np.vstack(list(no_srp_results.state_history.values()))[-1, :]
final_time = np.array(list(no_srp_results.state_history.keys()))[-1]
cartesian_states_results.append(np.hstack([final_time, final_car_state_no_srp]))


with open("cartesian results AE4868 2026 A1 6541151.txt", "a") as f:
    for row in cartesian_states_results:
        long_format_row = [f"{val:.16e}" for val in row]
        line = " ".join(map(str, long_format_row))
        f.write(f"{line}\n")






