###########################################################################
#
# # Numerical Astrodynamics 2024/2025
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

import os

import numpy as np
from matplotlib import pyplot as plt
from tudatpy import constants
from tudatpy.data import save2txt
from tudatpy.interface import spice
from tudatpy.dynamics import environment_setup, propagation_setup, simulator



current_task = 'question1'



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

###########################################################################
# CREATE ENVIRONMENT ######################################################
###########################################################################

# Load spice kernels.
spice.load_standard_kernels()
spice.load_kernel(current_directory + "/juice_mat_crema_5_1_150lb_v01.bsp")

# Create settings for celestial bodies
bodies_to_create = ['Ganymede']
global_frame_origin = "Ganymede"
global_frame_orientation = "ECLIPJ2000"
body_settings = environment_setup.get_default_body_settings(
    bodies_to_create, global_frame_origin, global_frame_orientation
)


###########################################################################
# CREATE VEHICLE ##########################################################
###########################################################################

body_settings.add_empty_settings("JUICE")

# Create environment
bodies = environment_setup.create_system_of_bodies(body_settings)

###########################################################################
# CREATE ACCELERATIONS ####################################################
###########################################################################

# Define bodies that are propagated, and their central bodies of propagation.
bodies_to_propagate = ["JUICE"]
central_bodies = ["Ganymede"]

# Define accelerations acting on vehicle.
acceleration_settings_on_vehicle = dict(
    Ganymede=[propagation_setup.acceleration.point_mass_gravity()]
)

# Create global accelerations dictionary.
acceleration_settings = {"JUICE": acceleration_settings_on_vehicle}

# Create acceleration models.
acceleration_models = propagation_setup.create_acceleration_models(
    bodies, acceleration_settings, bodies_to_propagate, central_bodies
)

###########################################################################
# CREATE PROPAGATION SETTINGS #############################################
###########################################################################

# Define initial state.
system_initial_state = spice.get_body_cartesian_state_at_epoch(
    target_body_name="JUICE",
    observer_body_name="Ganymede",
    reference_frame_name="ECLIPJ2000",
    aberration_corrections="NONE",
    ephemeris_time=simulation_start_epoch,
)

# Define required outputs
dependent_variables_to_save = [
    propagation_setup.dependent_variable.keplerian_state("JUICE", "Ganymede")
]

# Create numerical integrator settings.
fixed_step_size = 10.0
integrator_settings = propagation_setup.integrator.runge_kutta_fixed_step(
    fixed_step_size, coefficient_set=propagation_setup.integrator.CoefficientSets.rk_4
)

# Create propagation settings.
termination_settings = propagation_setup.propagator.time_termination(
    simulation_end_epoch
)
propagator_settings = propagation_setup.propagator.translational(
    central_bodies,
    acceleration_models,
    bodies_to_propagate,
    system_initial_state,
    simulation_start_epoch,
    integrator_settings,
    termination_settings,
    output_variables=dependent_variables_to_save,
)

propagator_settings.print_settings.print_initial_and_final_conditions = True


###########################################################################
# PROPAGATE ORBIT #########################################################
###########################################################################

# Create simulation object and propagate dynamics.
dynamics_simulator = simulator.create_dynamics_simulator(
    bodies, propagator_settings
)

# Retrieve all data produced by simulation
propagation_results = dynamics_simulator.propagation_results

# Extract numerical solution for states and dependent variables
state_history = propagation_results.state_history
dependent_variables = propagation_results.dependent_variable_history

###########################################################################
# SAVE RESULTS ############################################################
###########################################################################

save2txt(
    solution=state_history, filename="JUICEPropagationHistory_Q1.txt", directory=f"./Data/{current_task}"
)

save2txt(
    solution=dependent_variables,
    filename="JUICEPropagationHistory_DependentVariables_Q1.txt",
    directory=f"./Data/{current_task}",
)

###########################################################################
# PLOT RESULTS ############################################################
###########################################################################

# Extract time and Kepler elements from dependent variables
kepler_elements = np.vstack(list(dependent_variables.values()))
time = dependent_variables.keys()
time_days = [
    t / constants.JULIAN_DAY - simulation_start_epoch / constants.JULIAN_DAY
    for t in time
]




############################################################################################################################################
#                                                               TASK 1                                                                     # 
############################################################################################################################################
# My imports
import plotting
from tudatpy.astro import element_conversion

np.savetxt(f'./Data/{current_task}/epochs_days.txt', time_days)

final_car_state_PM = np.vstack(list(state_history.values()))[-1, :]
final_time_PM = np.array(list(state_history.keys()))[-1]

cartesian_states_results = []; cartesian_states_results.append(np.hstack([final_time_PM, final_car_state_PM]))






## Initial cartesian and Keplerian states:
initial_cartesian_state = system_initial_state
mu_ganymede=bodies.get("Ganymede").gravitational_parameter


initial_keplerian_state = element_conversion.cartesian_to_keplerian(initial_cartesian_state, mu_ganymede)

# Analytical solution 
kepler_ephemeris_settings = environment_setup.ephemeris.keplerian(
    initial_keplerian_state,
    simulation_start_epoch,
    mu_ganymede,
    'Ganymede',
    'ECLIPJ2000', 
)


kepler_ephemeris = environment_setup.create_body_ephemeris(kepler_ephemeris_settings, "JUICE")

epochs = np.array(list(state_history.keys()))

analytical_cartesian_state= {t: kepler_ephemeris.cartesian_state(t) for t in epochs}
analytical_cartesian_state = np.vstack(list(analytical_cartesian_state.values()))

keplerian_elements_analytical= []
length = analytical_cartesian_state.shape[0]

for idx in range(0, length):
    keplerian_elements_analytical.append(element_conversion.cartesian_to_keplerian(analytical_cartesian_state[idx, :], mu_ganymede))

keplerian_elements_analytical = np.array(keplerian_elements_analytical)
cartesian_integrated_state = np.vstack(list(state_history.values()))

# Computing the integration error
cartesian_error = {t: row for t, row in zip(epochs, cartesian_integrated_state - analytical_cartesian_state)}

keplerian_error = {t: row for t, row in zip(epochs, kepler_elements - keplerian_elements_analytical)}

save2txt(cartesian_error, 'cartesian_residual.txt', f'./Data/{current_task}')
save2txt(keplerian_error, 'keplerian_residual.txt', f'./Data/{current_task}')


r_a = np.max(np.linalg.norm(cartesian_integrated_state[:, 0:3], axis=1))
r_p = np.min(np.linalg.norm(cartesian_integrated_state[:, 0:3], axis=1))

print(f"Apocenter: {r_a:.16e} m")
print(f"Pericenter: {r_p:.16e} m")  








############################################################################################################################################
#                                                            DATA SAVE                                                                     # 
############################################################################################################################################
with open("cartesian results AE4868 2026 A1 6541151.txt", "a") as f:
    for row in cartesian_states_results:
        long_format_row = [f"{val:.16e}" for val in row]
        line = " ".join(map(str, long_format_row))
        f.write(f"{line}\n")

