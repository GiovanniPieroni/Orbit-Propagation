# import os
# import numpy as np
# from matplotlib import pyplot as plt

# from tudatpy import constants
# from tudatpy.data import save2txt
# from tudatpy.interface import spice
# from tudatpy.dynamics import environment_setup, propagation_setup, simulator
# from tudatpy.dynamics.environment_setup import rotation_model

# current_directory = os.getcwd()

# A = 7
# B = 5
# C = 3

# simulation_start_epoch = (
#     35.4 * constants.JULIAN_YEAR
#     + A * 7.0 * constants.JULIAN_DAY
#     + B * constants.JULIAN_DAY
#     + C * constants.JULIAN_DAY / 24.0
# )
# simulation_end_epoch = simulation_start_epoch + 344.0 * constants.JULIAN_DAY / 24.0

# fixed_step_size = 10.0

# ###########################################################################
# # SIMULATION FUNCTION SETUP ###############################################
# ###########################################################################

# def run_q4_simulation(model="Q1", center="Ganymede", concurrent=False):
#     """
#     Funzione per eseguire i diversi casi richiesti dalla Q4.
#     """
#     spice.load_standard_kernels()
#     spice.load_kernel(os.path.join(current_directory, "juice_mat_crema_5_1_150lb_v01.bsp"))

#     bodies_to_create = ["Jupiter", "Ganymede", "Sun", "Saturn", "Io", "Europa", "Callisto"]
#     global_frame_origin = center
#     global_frame_orientation = "ECLIPJ2000"

#     body_settings = environment_setup.get_default_body_settings(
#         bodies_to_create, global_frame_origin, global_frame_orientation
#     )

#     if model == "Q2":
#         rho0 = 2.0e-9
#         H = 40.0e3
#         body_settings.get("Ganymede").atmosphere_settings = environment_setup.atmosphere.exponential(H, rho0)
        
#         L_sun = 3.828e26
#         luminosity_settings = environment_setup.radiation_pressure.constant_luminosity(L_sun)
#         body_settings.get("Sun").radiation_source_settings = environment_setup.radiation_pressure.isotropic_radiation_source(luminosity_settings)

#     body_settings.add_empty_settings("JUICE")
#     body_settings.get("JUICE").constant_mass = 2.0e3

#     bodies = environment_setup.create_system_of_bodies(body_settings)

#     if model == "Q2":
#         reference_area = 100.0
#         CD = 1.2
#         Cr = 1.2
#         aero_coefficient_settings = environment_setup.aerodynamic_coefficients.constant(reference_area, [CD, 0.0, 0.0])
#         environment_setup.add_aerodynamic_coefficient_interface(bodies, "JUICE", aero_coefficient_settings)
        
#         occulting_bodies_dict = {"Sun": ["Ganymede"]}
#         juice_target_settings = environment_setup.radiation_pressure.cannonball_radiation_target(
#             reference_area=reference_area, radiation_pressure_coefficient=Cr, per_source_occulting_bodies=occulting_bodies_dict
#         )
#         environment_setup.add_radiation_pressure_target_model(bodies, "JUICE", juice_target_settings)

#     bodies_to_propagate = ["JUICE"]
#     central_bodies = [center]

#     if concurrent:
#         bodies_to_propagate.append("Ganymede")
#         central_bodies.append("Jupiter")

#     acc_dict = {"JUICE": {}}
#     if concurrent:
#         acc_dict["Ganymede"] = {}

#     if model == "Q1":
#         acc_dict["JUICE"]["Ganymede"] = [propagation_setup.acceleration.point_mass_gravity()]
        
#     elif model == "Q2":
#         acc_dict["JUICE"]["Ganymede"] = [propagation_setup.acceleration.spherical_harmonic_gravity(2, 2), propagation_setup.acceleration.aerodynamic()]
#         acc_dict["JUICE"]["Jupiter"] = [propagation_setup.acceleration.spherical_harmonic_gravity(4, 0)]
#         acc_dict["JUICE"]["Sun"] = [propagation_setup.acceleration.point_mass_gravity(), propagation_setup.acceleration.radiation_pressure()]
#         acc_dict["JUICE"]["Saturn"] = [propagation_setup.acceleration.point_mass_gravity()]
#         acc_dict["JUICE"]["Io"] = [propagation_setup.acceleration.point_mass_gravity()]
#         acc_dict["JUICE"]["Europa"] = [propagation_setup.acceleration.point_mass_gravity()]
#         acc_dict["JUICE"]["Callisto"] = [propagation_setup.acceleration.point_mass_gravity()]
        
#         if concurrent:
#             acc_dict["Ganymede"]["Jupiter"] = [propagation_setup.acceleration.spherical_harmonic_gravity(4, 0)]
#             acc_dict["Ganymede"]["Sun"] = [propagation_setup.acceleration.point_mass_gravity()]
#             acc_dict["Ganymede"]["Saturn"] = [propagation_setup.acceleration.point_mass_gravity()]
#             acc_dict["Ganymede"]["Io"] = [propagation_setup.acceleration.point_mass_gravity()]
#             acc_dict["Ganymede"]["Europa"] = [propagation_setup.acceleration.point_mass_gravity()]
#             acc_dict["Ganymede"]["Callisto"] = [propagation_setup.acceleration.point_mass_gravity()]

#     acceleration_models = propagation_setup.create_acceleration_models(
#         bodies, acc_dict, bodies_to_propagate, central_bodies
#     )

#     initial_state = spice.get_body_cartesian_state_at_epoch(
#         "JUICE", center, global_frame_orientation, "NONE", simulation_start_epoch
#     )
    
#     if concurrent:
#         initial_state_ganymede = spice.get_body_cartesian_state_at_epoch(
#             "Ganymede", "Jupiter", global_frame_orientation, "NONE", simulation_start_epoch
#         )
#         initial_state = np.concatenate((initial_state, initial_state_ganymede))

#     integrator_settings = propagation_setup.integrator.runge_kutta_fixed_step(
#         fixed_step_size, propagation_setup.integrator.CoefficientSets.rk_4
#     )
#     termination_settings = propagation_setup.propagator.time_termination(simulation_end_epoch)

#     propagator_settings = propagation_setup.propagator.translational(
#         central_bodies, acceleration_models, bodies_to_propagate, initial_state,
#         simulation_start_epoch, integrator_settings, termination_settings
#     )

#     dyn = simulator.create_dynamics_simulator(bodies, propagator_settings)
#     prop_results = dyn.propagation_results

#     # Data Extraction
#     times = np.array(list(prop_results.state_history.keys()))
#     state_history_array = np.vstack(list(prop_results.state_history.values()))
    
#     if center == "Ganymede":
#         states_wrt_ganymede = state_history_array[:, 0:6]
        
#     elif center == "Jupiter":
#         states_wrt_ganymede = np.zeros_like(state_history_array[:, 0:6])
#         if concurrent:
#             x_J_s = state_history_array[:, 0:6]
#             x_J_G_prop = state_history_array[:, 6:12]
#             states_wrt_ganymede = x_J_s - x_J_G_prop
#         else:
#             for i in range(len(times)):
#                 x_J_G_spice = spice.get_body_cartesian_state_at_epoch(
#                     "Ganymede", "Jupiter", global_frame_orientation, "NONE", times[i])
#                 states_wrt_ganymede[i, :] = state_history_array[i, 0:6] - x_J_G_spice

#     return times, states_wrt_ganymede, prop_results

# ###########################################################################
# # PROPAGATION RUNNING (in all cases)#######################################
# ###########################################################################

# t_i, s_i, _ = run_q4_simulation("Q1", "Ganymede", concurrent=False)
# t_iii, s_iii, _ = run_q4_simulation("Q2", "Ganymede", concurrent=False)

# t_ii_bad, s_ii_bad, _ = run_q4_simulation("Q1", "Jupiter", concurrent=False)
# t_iv_bad, s_iv_bad, _ = run_q4_simulation("Q2", "Jupiter", concurrent=False)

# t_ii_good, s_ii_good, _ = run_q4_simulation("Q1", "Jupiter", concurrent=True)
# t_iv_good, s_iv_good, results_q4_nominal = run_q4_simulation("Q2", "Jupiter", concurrent=True)

# ###########################################################################
# # DATA EXTRACTION AND ERRORS CALCULATION ##################################
# ###########################################################################

# time_hours = (t_i - t_i[0]) / 3600.0

# # Errors in i-ii case (first formulation)
# err_q1_bad = np.linalg.norm(s_i[:, 0:3] - s_ii_bad[:, 0:3], axis=1)
# err_q2_bad = np.linalg.norm(s_iii[:, 0:3] - s_iv_bad[:, 0:3], axis=1)

# # Errors in iii-iv case (second formulation)
# err_q1_good = np.linalg.norm(s_i[:, 0:3] - s_ii_good[:, 0:3], axis=1)
# err_q2_good = np.linalg.norm(s_iii[:, 0:3] - s_iv_good[:, 0:3], axis=1)









#######################################################################################################
import os

import numpy as np
from matplotlib import pyplot as plt
from tudatpy import constants
from tudatpy.data import save2txt
from tudatpy.interface import spice
from tudatpy.dynamics import environment_setup, propagation_setup, simulator
from tudatpy.util import result2array

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

def run_scenario_q1(cb):
    # Create settings for celestial bodies
    bodies_to_create = ["Ganymede", "Jupiter"]
    global_frame_origin = cb
    global_frame_orientation = "ECLIPJ2000"
    body_settings = environment_setup.get_default_body_settings(
        bodies_to_create, global_frame_origin, global_frame_orientation
    )

    body_settings.add_empty_settings("JUICE")

    # Create environment
    bodies = environment_setup.create_system_of_bodies(body_settings)

    # Define bodies that are propagated, and their central bodies of propagation.
    bodies_to_propagate = ["JUICE"]
    central_bodies = [cb]

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

    # Define initial state.
    system_initial_state = spice.get_body_cartesian_state_at_epoch(
        target_body_name="JUICE",
        observer_body_name= cb,
        reference_frame_name="ECLIPJ2000",
        aberration_corrections="NONE",
        ephemeris_time=simulation_start_epoch,
    )

    # Define required outputs
    dependent_variables_to_save = [
        propagation_setup.dependent_variable.relative_position("Ganymede", "Jupiter"),
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

    # Create simulation object and propagate dynamics.
    dynamics_simulator = simulator.create_dynamics_simulator(
        bodies, propagator_settings
    )

    # Retrieve all data produced by simulation
    propagation_results = dynamics_simulator.propagation_results

    # Extract numerical solution for states and dependent variables
    state_history = propagation_results.state_history
    dependent_variables = propagation_results.dependent_variable_history

    dep_vars_history = dynamics_simulator.propagation_results.dependent_variable_history
    dep_vars_array = result2array(dep_vars_history)
    states_history = dynamics_simulator.propagation_results.state_history
    states_array = result2array(states_history)
    time = dependent_variables.keys()
    time_days = [
        t / constants.JULIAN_DAY - simulation_start_epoch / constants.JULIAN_DAY
        for t in time
    ]


    if cb == "Ganymede" :
        position = states_array[:, 1:4]
    else :
        r_Js = states_array[:, 1:4]
        r_JG = dep_vars_array[:, 1:4]
        position = r_Js - r_JG

    position = np.linalg.norm(position, axis=1)

    return position, time_days

total_position_11, time_days = run_scenario_q1("Ganymede")
total_position_12, time_days = run_scenario_q1("Jupiter")

delta_total_position = total_position_11 - total_position_12

plt.figure(figsize=(10, 8))
plt.plot(time_days, np.linalg.norm(delta_total_position, axis=1)   )
plt.xlabel("Time (days)")
plt.ylabel("Distance (m)")
plt.title("Distance between JUICE and Ganymede vs Jupiter")
plt.grid(True)
plt.show()

def run_scenario_q2(cb):
    # Create settings for celestial bodies
    bodies_to_create = ["Ganymede", "Jupiter", "Saturn", "Sun", "Europa", "Io", "Callisto"]
    global_frame_origin = cb
    global_frame_orientation = "ECLIPJ2000"
    body_settings = environment_setup.get_default_body_settings(
        bodies_to_create, global_frame_origin, global_frame_orientation
    )

    # Ganymede atmosphere settings
    density_scale_height = 40.0E3
    density_at_zero_altitude = 2.0E-9
    body_settings.get("Ganymede").atmosphere_settings = environment_setup.atmosphere.exponential(
        density_scale_height, density_at_zero_altitude)

    # ----- Create vehicle -----

    body_settings.add_empty_settings("JUICE")

    body_settings.get("JUICE").constant_mass = 2.0E3

    # Aerodynamic settings of JUICE
    drag_coefficient = 1.2
    drag_reference_area = 100.0
    aerodynamic_settings = environment_setup.aerodynamic_coefficients.constant(
        drag_reference_area, [drag_coefficient, 0.0, 0.0])

    body_settings.get("JUICE").aerodynamic_coefficient_settings = aerodynamic_settings

    # Radiation pressure settings of JUICE
    radiation_pressure_coefficient = 1.2
    reference_area_radiation = 100.0
    occulting_bodies_dict = dict()
    occulting_bodies_dict["Sun"] = ["Ganymede"]
    vehicle_target_settings = environment_setup.radiation_pressure.cannonball_radiation_target(
        reference_area_radiation, radiation_pressure_coefficient, occulting_bodies_dict)

    body_settings.get("JUICE").radiation_pressure_target_settings = vehicle_target_settings

    # Create system of bodies
    bodies = environment_setup.create_system_of_bodies(body_settings)

    # ----- Create Accelerations -----

    bodies_to_propagate = ["JUICE"]
    central_bodies = [cb]

    acceleration_settings_on_vehicle = dict(
        Ganymede=
        [
            propagation_setup.acceleration.spherical_harmonic_gravity(2, 2),
            propagation_setup.acceleration.aerodynamic()
        ],
        Jupiter=[
            propagation_setup.acceleration.spherical_harmonic_gravity(4, 0)
        ],
        Sun=[
            propagation_setup.acceleration.point_mass_gravity(),
            propagation_setup.acceleration.radiation_pressure()
        ],
        Saturn=[
            propagation_setup.acceleration.point_mass_gravity()
        ],
        Europa=[
            propagation_setup.acceleration.point_mass_gravity()
        ],
        Io=[
            propagation_setup.acceleration.point_mass_gravity()
        ],
        Callisto=[
            propagation_setup.acceleration.point_mass_gravity()
        ]
    )

    # Create global accelerations dictionary.
    acceleration_settings = {"JUICE": acceleration_settings_on_vehicle}

    # Create acceleration models.
    acceleration_models = propagation_setup.create_acceleration_models(
        bodies, acceleration_settings, bodies_to_propagate, central_bodies
    )

    # ----- Propagation Settings -----

    # Define start and stop epoch

    # # student number: 1244779 --> 1244ABC
    A = 1
    B = 5
    C = 4

    simulation_start_epoch = (
            35.4 * constants.JULIAN_YEAR
            + A * 7.0 * constants.JULIAN_DAY
            + B * constants.JULIAN_DAY
            + C * constants.JULIAN_DAY / 24.0
    )
    simulation_end_epoch = simulation_start_epoch + 344.0 * constants.JULIAN_DAY / 24.0

    # Define initial state.
    system_initial_state = spice.get_body_cartesian_state_at_epoch(
        target_body_name="JUICE",
        observer_body_name=cb,
        reference_frame_name="ECLIPJ2000",
        aberration_corrections="NONE",
        ephemeris_time=simulation_start_epoch,
    )

    # Define required outputs
    dependent_variables_to_save = [
        propagation_setup.dependent_variable.relative_position("Ganymede", "Jupiter"),
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

    # ----- Propagate the orbit -----

    # Create simulation object and propagate dynamics.
    dynamics_simulator = simulator.create_dynamics_simulator(
        bodies, propagator_settings
    )


    # Retrieve all data produced by simulation
    propagation_results = dynamics_simulator.propagation_results

    # Extract numerical solution for states and dependent variables
    state_history = propagation_results.state_history
    dependent_variables = propagation_results.dependent_variable_history

    dep_vars_history = dynamics_simulator.propagation_results.dependent_variable_history
    dep_vars_array = result2array(dep_vars_history)
    states_history = dynamics_simulator.propagation_results.state_history
    states_array = result2array(states_history)
    time = dependent_variables.keys()
    time_days = [
        t / constants.JULIAN_DAY - simulation_start_epoch / constants.JULIAN_DAY
        for t in time
    ]


    if cb == "Ganymede" :
        position = states_array[:, 1:4]
    else :
        r_Js = states_array[:, 1:4]
        r_JG = dep_vars_array[:, 1:4]
        position = r_Js - r_JG

    position = np.linalg.norm(position, axis=1)

    return position, time_days

total_position_21, time_days = run_scenario_q2("Ganymede")
total_position_22, time_days = run_scenario_q2("Jupiter")

delta_total_position_2 = total_position_21 - total_position_22

plt.figure(figsize=(10, 8))
plt.plot(time_days, np.linalg.norm(delta_total_position_2, axis=1))
plt.xlabel("Time (days)")
plt.ylabel("Distance (m)")
plt.title("Distance between JUICE and Ganymede vs Jupiter")
plt.grid(True)
plt.show()