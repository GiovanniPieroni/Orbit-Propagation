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

import os

import numpy as np
from matplotlib import pyplot as plt
from tudatpy import constants
from tudatpy.data import save2txt
from tudatpy.interface import spice
from tudatpy.dynamics import environment_setup, propagation_setup, simulator


current_task = 'question2'



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
bodies_to_create = ['Ganymede', 'Jupiter', 'Sun', 'Saturn', 'Europa', 'Io', 'Callisto']
global_frame_origin = "Ganymede"
global_frame_orientation = "ECLIPJ2000"
body_settings = environment_setup.get_default_body_settings(
    bodies_to_create, global_frame_origin, global_frame_orientation
)

# Modifying atmosphere model -> H = 40 km, rho_0 = 2*10^-9 kg/m^3
body_settings.get( 'Ganymede' ).atmosphere_settings = environment_setup.atmosphere.exponential(
    scale_height = 40000.0,
    surface_density = 2e-9 )



# Sun irradiance definition
irradiance_at_1AU = 1367.0  # W/m^2, Vallado 2013

luminosity_model_settings = (
    environment_setup.radiation_pressure.irradiance_based_constant_luminosity(
        irradiance_at_1AU, constants.ASTRONOMICAL_UNIT
    )
)
radiation_source_settings_sun = (
    environment_setup.radiation_pressure.isotropic_radiation_source(
        luminosity_model_settings
    )
)

body_settings.get("Sun").radiation_source_settings = radiation_source_settings_sun

eclipse_settings = {
    'Sun': ['Ganymede']
}

###########################################################################
# CREATE VEHICLE ##########################################################
###########################################################################

body_settings.add_empty_settings("JUICE")

# Add aerodynamic model settings -> Reference area = 100 m^2, C_D = 1.2
body_settings.get( 'JUICE' ).aerodynamic_coefficient_settings = environment_setup.aerodynamic_coefficients.constant(
    reference_area = 100.0,
    constant_force_coefficient = [1.2, 0.0, 0.0])


body_settings.get('JUICE').radiation_pressure_target_settings = environment_setup.radiation_pressure.cannonball_radiation_target(
    reference_area = 100,
    radiation_pressure_coefficient = 1.2, 
    per_source_occulting_bodies = eclipse_settings
)







# Create environment
bodies = environment_setup.create_system_of_bodies(body_settings)
bodies.get('JUICE').mass = 2e3

###########################################################################
# CREATE ACCELERATIONS ####################################################
###########################################################################

# Define bodies that are propagated, and their central bodies of propagation.
bodies_to_propagate = ["JUICE"]
central_bodies = ["Ganymede"]

# Define accelerations acting on vehicle.
acceleration_settings_on_vehicle = dict(
    Ganymede=[propagation_setup.acceleration.spherical_harmonic_gravity(2,2),
              propagation_setup.acceleration.aerodynamic()
              ],

    Jupiter=[propagation_setup.acceleration.spherical_harmonic_gravity(4,0)],

    Sun=[propagation_setup.acceleration.point_mass_gravity(),
         propagation_setup.acceleration.radiation_pressure()
         ],

    Saturn  =[propagation_setup.acceleration.point_mass_gravity()],
    Europa  =[propagation_setup.acceleration.point_mass_gravity()],
    Io      =[propagation_setup.acceleration.point_mass_gravity()],
    Callisto=[propagation_setup.acceleration.point_mass_gravity()]
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
    propagation_setup.dependent_variable.keplerian_state("JUICE", "Ganymede"),  # Keplerian state of JUICE
    propagation_setup.dependent_variable.single_acceleration(
        propagation_setup.acceleration.spherical_harmonic_gravity_type, 'JUICE', 'Ganymede' ),  # Gravitational acceleration of Ganymede
    propagation_setup.dependent_variable.single_acceleration(
        propagation_setup.acceleration.point_mass_gravity_type, 'JUICE', 'Io' ),    # Gravitational acceleration of Io
        propagation_setup.dependent_variable.single_acceleration(
        propagation_setup.acceleration.point_mass_gravity_type, 'JUICE', 'Sun' ),   # Gravitational acceleration of the Sun
    propagation_setup.dependent_variable.single_acceleration(
        propagation_setup.acceleration.aerodynamic_type, 'JUICE', 'Ganymede' ),       # Aerodynamic acceleration due to Ganymede exosphere
    propagation_setup.dependent_variable.single_acceleration(
        propagation_setup.acceleration.radiation_pressure_type, 'JUICE', 'Sun')     # Radiation pressure acceleration of the Sun
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
    solution=state_history, filename="JUICEPropagationHistory_Q2.txt", directory=f"./Data/{current_task}"
)

# Dependent variables order: a, e, i, Om, om, theta, a_grav_G (3), a_grav_Io (3), a_grav_S (3), a_aer_G (3), a_rad_S (3)
save2txt(
    solution=dependent_variables,
    filename="JUICEPropagationHistory_DependentVariables_Q2.txt",       
    directory=f"./Data/{current_task}",
)


# Extract time and Kepler elements from dependent variables
kepler_elements = np.vstack(list(dependent_variables.values()))[:, 0:6]
time = dependent_variables.keys()
time_days = [
    t / constants.JULIAN_DAY - simulation_start_epoch / constants.JULIAN_DAY
    for t in time
]

np.savetxt(f'./Data/{current_task}/epochs_days.txt', time_days)






############################################################################################################################################
#                                                               TASK C (i)                                                                 # 
############################################################################################################################################
final_car_state_PM = np.vstack(list(state_history.values()))[-1, :]
final_time_PM = np.array(list(state_history.keys()))[-1]

cartesian_states_results = []; cartesian_states_results.append(np.hstack([final_time_PM, final_car_state_PM]))





############################################################################################################################################
#                                                               TASK C (ii)                                                                # 
############################################################################################################################################
import scipy as sp

Io_grav_acc = np.vstack(list(dependent_variables.values()))[:, 10:13]
Sun_grav_acc = np.vstack(list(dependent_variables.values()))[:, 13:16]

Io_grav_acc_norm = np.linalg.norm(Io_grav_acc, axis=1)
Sun_grav_acc_norm = np.linalg.norm(Sun_grav_acc, axis=1)

Io_accel_no_mean = Io_grav_acc_norm -  np.mean(Io_grav_acc_norm)  # Removing mean acceleration to focus on variations
Sun_accel_no_mean = Sun_grav_acc_norm -  np.mean(Sun_grav_acc_norm)  # Removing mean acceleration to focus on variations


Io_fft = sp.fftpack.fft(Io_accel_no_mean)
Io_power_density = np.abs(Io_fft)**2
Io_freqs = sp.fftpack.fftfreq(len(Io_accel_no_mean), d=10)

Io_peak_idx = np.argmax(Io_power_density)
Io_peak_freq = Io_freqs[Io_peak_idx]
Io_period_seconds = 1.0 / Io_peak_freq
Io_period_hours = Io_period_seconds / 3600.0    
print(f"Io's period of the largest variation : {Io_period_hours:.2f} hours")


Sun_fft = sp.fftpack.fft(Sun_accel_no_mean)
Sun_power_density = np.abs(Sun_fft)**2
Sun_freqs = sp.fftpack.fftfreq(len(Sun_accel_no_mean), d=10)

Sun_peak_idx = np.argmax(Sun_power_density)
Sun_peak_freq = Sun_freqs[Sun_peak_idx]
Sun_period_seconds = 1.0 / Sun_peak_freq
Sun_period_hours = Sun_period_seconds / 3600.0    
print(f"Sun's period of the largest variation : {Sun_period_hours:.2f} hours")


np.savetxt(f'./Data/{current_task}/Io_frequency_spectrum.txt', np.column_stack((Io_freqs, Io_power_density)), header='Frequency(Hz) PowerDensity')
np.savetxt(f'./Data/{current_task}/Sun_frequency_spectrum.txt', np.column_stack((Sun_freqs, Sun_power_density)), header='Frequency(Hz) PowerDensity')



############################################################################################################################################
#                                                               TASK C (iii)                                                               # 
############################################################################################################################################
Io_period_hr = 42.5
G_period_hr = 171.7
sma_mean = np.mean(kepler_elements[:, 0])
mu_ganymede = bodies.get('Ganymede').gravitational_parameter
JUICE_period_hr = 2*np.pi*np.sqrt(sma_mean**3/mu_ganymede)/3600

print(f'JUICE orbital period: {JUICE_period_hr:.8f} hours')
print(f'Io orbital period around Jupiter: {Io_period_hr:.2f} hours')
print(f'Ganymede orbital period around Jupiter: {G_period_hr:.2f} hours')

T_syn_hr_Io_Ganymede = 1 / abs(1/Io_period_hr - 1/G_period_hr)
print('Synodic period between Io and Ganymede: {:.2f} hours'.format(T_syn_hr_Io_Ganymede))


















############################################################################################################################################
#                                                            DATA SAVE                                                                     # 
############################################################################################################################################
with open("cartesian results AE4868 2026 A1 6541151.txt", "a") as f:
    for row in cartesian_states_results:
        long_format_row = [f"{val:.16e}" for val in row]
        line = " ".join(map(str, long_format_row))
        f.write(f"{line}\n")
