###########################################################################
#
# # Numerical Astrodynamics 2025/2026
#
# # Assignment 1 - Propagation Settings
#
###########################################################################


""" 
"Jupiter"yright (c) 2010-2020, Delft University of Technology
All rights reserved

This file is part of Tudat. Redistribution and use in source and 
binary forms, with or without modification, are permitted exclusively
under the terms of the Modified BSD license. You should have received
a "Jupiter"y of the license with this file. If not, please or visit:
http://tudat.tudelft.nl/LICENSE.
"""

import numpy as np

from tudatpy import constants
from tudatpy.interface import spice
from tudatpy.dynamics import environment_setup, propagation_setup, simulator
from tudatpy.data import save2txt

import os




def run_juice_simulation_perturbed(
    simulation_start_epoch, 
    simulation_end_epoch, 
    fixed,
    current_directory
):
    # Load spice kernels
    spice.load_standard_kernels()
    spice.load_kernel(current_directory + "/juice_mat_crema_5_1_150lb_v01.bsp")

    # Create settings for celestial bodies
    bodies_to_create = ['Ganymede', 'Jupiter', 'Sun', 'Saturn', 'Europa', 'Io', 'Callisto']
    global_frame_origin = "Jupiter"
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
        reference_area=100.0, constant_force_coefficient=[1.2, 0.0, 0.0])
    body_settings.get('JUICE').radiation_pressure_target_settings = environment_setup.radiation_pressure.cannonball_radiation_target(
        reference_area=100, radiation_pressure_coefficient=1.2, 
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

    
   
    acceleration_settings_on_ganymede = dict(
        Jupiter=[propagation_setup.acceleration.spherical_harmonic_gravity(4,0)],
        Sun=[propagation_setup.acceleration.point_mass_gravity()],
        Saturn=[propagation_setup.acceleration.point_mass_gravity()],
        Europa=[propagation_setup.acceleration.point_mass_gravity()],
        Io=[propagation_setup.acceleration.point_mass_gravity()],
        Callisto=[propagation_setup.acceleration.point_mass_gravity()]
    )
    
     

    # Initial state
    juice_initial_state = spice.get_body_cartesian_state_at_epoch(
        "JUICE", "Jupiter", "ECLIPJ2000", "NONE", simulation_start_epoch)
     
    
    if not fixed :
        ganymede_initial_state = spice.get_body_cartesian_state_at_epoch(
        "Ganymede", "Jupiter", "ECLIPJ2000", "NONE", simulation_start_epoch) 
        system_initial_state = np.concatenate((juice_initial_state, ganymede_initial_state))
    else:
        system_initial_state = juice_initial_state


    integrator_settings = propagation_setup.integrator.runge_kutta_fixed_step(
        10.0, coefficient_set=propagation_setup.integrator.CoefficientSets.rk_4)

    
    

    # Output variables
    dependent_variables_to_save = [
        propagation_setup.dependent_variable.keplerian_state("JUICE", "Jupiter")
    ]

    dependent_variables_to_save.append(propagation_setup.dependent_variable.relative_position("Ganymede", "Jupiter"))
    dependent_variables_to_save.append(propagation_setup.dependent_variable.relative_velocity("Ganymede", "Jupiter"))
    
    if not fixed:

        acceleration_models = propagation_setup.create_acceleration_models(
            bodies, 
            {"JUICE": acceleration_settings_on_vehicle, "Ganymede": acceleration_settings_on_ganymede}, 
            ["JUICE", "Ganymede"], ["Jupiter", "Jupiter"])
    
        propagator_settings = propagation_setup.propagator.translational(
            ["Jupiter", "Jupiter"], acceleration_models, 
            ["JUICE", "Ganymede"], system_initial_state,
            simulation_start_epoch, integrator_settings, 
            propagation_setup.propagator.time_termination(simulation_end_epoch),
            output_variables=dependent_variables_to_save)
    else:

        acceleration_models = propagation_setup.create_acceleration_models(
        bodies, {"JUICE": acceleration_settings_on_vehicle}, ["JUICE"], ["Jupiter"]
    )
        propagator_settings = propagation_setup.propagator.translational(
            ["Jupiter"], acceleration_models, ["JUICE"], system_initial_state,
            simulation_start_epoch, integrator_settings, 
            propagation_setup.propagator.time_termination(simulation_end_epoch),
            output_variables=dependent_variables_to_save)

        

    dynamics_simulator = simulator.create_dynamics_simulator(bodies, propagator_settings)
    
    return dynamics_simulator.propagation_results


def run_juice_simulation_unperturbed(
    simulation_start_epoch, 
    simulation_end_epoch, 
    fixed,
    current_directory
):
    # Load spice kernels
    spice.load_standard_kernels()
    spice.load_kernel(current_directory + "/juice_mat_crema_5_1_150lb_v01.bsp")

    # Create settings for celestial bodies
    if not fixed:
        bodies_to_create = ["Jupiter", "Ganymede", "Sun", "Saturn", "Io", "Europa", "Callisto"]
    else:
        bodies_to_create = ['Ganymede', 'Jupiter']
    

    global_frame_origin = "Jupiter"
    global_frame_orientation = "ECLIPJ2000"
    body_settings = environment_setup.get_default_body_settings(
        bodies_to_create, global_frame_origin, global_frame_orientation
    )

    


    # Create Vehicle
    body_settings.add_empty_settings("JUICE")
    body_settings.get("JUICE").constant_mass = 2.0e3
    
    bodies = environment_setup.create_system_of_bodies(body_settings)



    # Accelerations
    acceleration_settings_on_vehicle = dict(
        Ganymede=[propagation_setup.acceleration.point_mass_gravity()]
    )


    # Initial state
    juice_initial_state = spice.get_body_cartesian_state_at_epoch(
        "JUICE", "Jupiter", "ECLIPJ2000", "NONE", simulation_start_epoch)
     
    
    if not fixed:
        ganymede_initial_state = spice.get_body_cartesian_state_at_epoch(
        "Ganymede", "Jupiter", "ECLIPJ2000", "NONE", simulation_start_epoch) 
        system_initial_state = np.concatenate((juice_initial_state, ganymede_initial_state))
    else:
        system_initial_state = juice_initial_state



    # Output variables
    dependent_variables_to_save = [
        propagation_setup.dependent_variable.keplerian_state("JUICE", "Jupiter")
    ]
    dependent_variables_to_save.append(propagation_setup.dependent_variable.relative_position("Ganymede", "Jupiter"))
    dependent_variables_to_save.append(propagation_setup.dependent_variable.relative_velocity("Ganymede", "Jupiter"))

    integrator_settings = propagation_setup.integrator.runge_kutta_fixed_step(
        10.0, coefficient_set=propagation_setup.integrator.CoefficientSets.rk_4)

    



    if not fixed:

        acceleration_settings_on_ganymede = {}
        

        # acceleration_settings_on_ganymede = dict(
        #     Jupiter=[propagation_setup.acceleration.spherical_harmonic_gravity(4,0)],
        #     Sun=[propagation_setup.acceleration.point_mass_gravity()],
        #     Saturn=[propagation_setup.acceleration.point_mass_gravity()],
        #     Europa=[propagation_setup.acceleration.point_mass_gravity()],
        #     Io=[propagation_setup.acceleration.point_mass_gravity()],
        #     Callisto=[propagation_setup.acceleration.point_mass_gravity()]
        #     )

        acceleration_models = propagation_setup.create_acceleration_models(
            bodies, 
            {"JUICE": acceleration_settings_on_vehicle, "Ganymede": acceleration_settings_on_ganymede}, 
            ["JUICE", "Ganymede"], ["Jupiter", "Jupiter"])

        # acceleration_models = propagation_setup.create_acceleration_models(
        #     bodies, 
        #     {"JUICE": acceleration_settings_on_vehicle}, 
        #     ["JUICE", "Ganymede"], ["Jupiter", "Jupiter"])

        propagator_settings = propagation_setup.propagator.translational(
            ["Jupiter", "Jupiter"], acceleration_models, ["JUICE", "Ganymede"], system_initial_state,
            simulation_start_epoch, integrator_settings, 
            propagation_setup.propagator.time_termination(simulation_end_epoch),
            output_variables=dependent_variables_to_save)
    else:

        acceleration_models = propagation_setup.create_acceleration_models(
        bodies, {"JUICE": acceleration_settings_on_vehicle}, ["JUICE"], ["Jupiter"])

        propagator_settings = propagation_setup.propagator.translational(
            ["Jupiter"], acceleration_models, ["JUICE"], system_initial_state,
            simulation_start_epoch, integrator_settings, 
            propagation_setup.propagator.time_termination(simulation_end_epoch),
            output_variables=dependent_variables_to_save)

    dynamics_simulator = simulator.create_dynamics_simulator(bodies, propagator_settings)
    
    return dynamics_simulator.propagation_results



current_task = 'question4'



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

simulation_cases = ['Unperturbed_Jupiter_centered', 
                    'Perturbed_Jupiter_centered', 
                    'Unperturbed_Jupiter_centered_Ganymede', 
                    'Perturbed_Jupiter_centered_Ganymede']



all_case_results = {}

for case in simulation_cases:

    print(f"Running simulation case: {case}")
    
    
    
    # Run simulation
    if "Ganymede" not in case and "Perturbed" not in case:
        results = run_juice_simulation_unperturbed(
            simulation_start_epoch, 
            simulation_end_epoch, 
            True, 
            current_directory
        )
    elif "Ganymede" not in case and "Perturbed" in case:
        results = run_juice_simulation_perturbed(
            simulation_start_epoch, 
            simulation_end_epoch, 
            True, 
            current_directory
        )
    elif "Ganymede" in case and "Perturbed" not in case:    
        results = run_juice_simulation_unperturbed(
            simulation_start_epoch, 
            simulation_end_epoch, 
            False, 
            current_directory
        )   
    elif "Ganymede" in case and "Perturbed" in case:
        results = run_juice_simulation_perturbed(
            simulation_start_epoch, 
            simulation_end_epoch, 
            False, 
            current_directory
    )
    
    all_case_results[case] = results

    # Save results specific to each case
    save2txt(results.state_history,
              f"JUICEPropagationHistory_Q4_{case}.txt",
                f"./Data/{current_task}")
    
    save2txt(results.dependent_variable_history,
              f"JUICEPropagationHistory_DependentVariables_Q4_{case}.txt",
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
UJ_results =   all_case_results['Unperturbed_Jupiter_centered']
PJ_results =   all_case_results['Perturbed_Jupiter_centered']   
UJ_G_results = all_case_results['Unperturbed_Jupiter_centered_Ganymede']
PJ_G_results = all_case_results['Perturbed_Jupiter_centered_Ganymede']

UJ_cartesian_state_Gcentered = np.array(list(UJ_results.state_history.values()))[:, :6] - np.array(list(UJ_results.dependent_variable_history.values()))[:, 6:12]

PJ_cartesian_state_Gcentered = np.array(list(PJ_results.state_history.values()))[:, :6] - np.array(list(PJ_results.dependent_variable_history.values()))[:, 6:12]

# UJ_G_cartesian_state_Gcentered = np.array(list(UJ_G_results.state_history.values()))[:, :6] \
#                                 - np.array(list(UJ_G_results.dependent_variable_history.values()))[:, 6:12]

# PJ_G_cartesian_state_Gcentered = np.array(list(PJ_G_results.state_history.values()))[:, :6] \
#                                 - np.array(list(PJ_G_results.dependent_variable_history.values()))[:, 6:12] 

UJ_G_cartesian_state_Gcentered = np.array(list(UJ_G_results.state_history.values()))[:, :6]  - np.array(list(UJ_G_results.state_history.values()))[:, 6:12]

PJ_G_cartesian_state_Gcentered = np.array(list(PJ_G_results.state_history.values()))[:, :6] - np.array(list(PJ_G_results.state_history.values()))[:, 6:12] 


np.savetxt(f'./Data/{current_task}/UJ_cartesian_state_Gcentered.txt', UJ_cartesian_state_Gcentered)
np.savetxt(f'./Data/{current_task}/PJ_cartesian_state_Gcentered.txt', PJ_cartesian_state_Gcentered)

np.savetxt(f'./Data/{current_task}/UJ_G_cartesian_state_Gcentered.txt', UJ_G_cartesian_state_Gcentered)
np.savetxt(f'./Data/{current_task}/PJ_G_cartesian_state_Gcentered.txt', PJ_G_cartesian_state_Gcentered)





############################################################################################################################################
#                                                            DATA SAVE                                                                     # 
############################################################################################################################################
# Saving results for final state comparison 
cartesian_states_results = []

final_car_state_UJ = UJ_cartesian_state_Gcentered[-1, :]
final_time = np.array(list(UJ_results.state_history.keys()))[-1]

cartesian_states_results.append(np.hstack([final_time, final_car_state_UJ]))

final_car_state_PJ = PJ_cartesian_state_Gcentered[-1, :]
final_time = np.array(list(PJ_results.state_history.keys()))[-1]

cartesian_states_results.append(np.hstack([final_time, final_car_state_PJ]))


with open("cartesian results AE4868 2026 A1 6541151.txt", "a") as f:
    for row in cartesian_states_results:
        long_format_row = [f"{val:.16e}" for val in row]
        line = " ".join(map(str, long_format_row))
        f.write(f"{line}\n")