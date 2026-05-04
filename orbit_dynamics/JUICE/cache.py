
time = results.state_history.keys()

time_days = [
    t / constants.JULIAN_DAY - simulation_start_epoch / constants.JULIAN_DAY
    for t in time
]

np.savetxt(f'./Data/{current_task}/epochs_days.dat', time_days)


####################################################################################################
# Post-Processing 
####################################################################################################

UJ_results =   all_case_results['Unperturbed_Jupiter_centered']
PJ_results =   all_case_results['Perturbed_Jupiter_centered']   
UJ_G_results = all_case_results['Unperturbed_Jupiter_centered_Ganymede']
PJ_G_results = all_case_results['Perturbed_Jupiter_centered_Ganymede']

UJ_cartesian_state_Gcentered = np.array(list(UJ_results.state_history.values()))[:, :3] - np.array(list(UJ_results.dependent_variable_history.values()))[:, 6:9]
PJ_cartesian_state_Gcentered = np.array(list(PJ_results.state_history.values()))[:, :3] - np.array(list(PJ_results.dependent_variable_history.values()))[:, 6:9]

np.savetxt(f'./Data/{current_task}/UJ_cartesian_state_Gcentered.dat', UJ_cartesian_state_Gcentered)
np.savetxt(f'./Data/{current_task}/PJ_cartesian_state_Gcentered.dat', PJ_cartesian_state_Gcentered)





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






# question 4 pp
    # U_JG_rel_pos = np.loadtxt(f'./Data/question4/JUICEPropagationHistory_DependentVariables_Q4_Unperturbed_Jupiter_centered.dat')[:, 7:10]
    # P_JG_rel_pos = np.loadtxt(f'./Data/question4/JUICEPropagationHistory_DependentVariables_Q4_Perturbed_Jupiter_centered.dat')[:, 7:10]

    # pt.plot_normal(time_days, np.linalg.norm(U_JG_rel_pos, axis=1), 'Time [days]', r'Relative Position Norm [m]',
    #                 "Relative Position between JUICE and Jupiter in Jupiter-centered Unperturbed Simulation", 
    #                 'rel_pos_UJ.pdf', None, 'question4', lwidth=2)  


    # pt.plot_two_scales(time_days, 
    #                    np.linalg.norm(UJ_pos, axis=1), 'Abs pos [m]', r'$||r_{JUICE}||$',
    #                    np.linalg.norm(U_JG_rel_pos, axis=1), 'Rel pos JG [m]', r'$||\Delta r_{JG}||$',
    #                    title="Relative Position between JUICE and Jupiter vs Absolute Position of JUICE in Jupiter-centered Unperturbed Simulation", 
    #                    filename='rel_pos_UJ_vs_abs_pos_UJ.pdf',
    #                    task='question4') 


    # pt.plot_normal(time_days, np.linalg.norm(UG_pos, axis=1), 'Time [days]', r'Absolute Position Norm [m]',
    #                 "Absolute Position of JUICE in Ganymede-centered Unperturbed Simulation", 
    #                 'abs_pos_UG.pdf', None, 'question4', lwidth=1)
    
    # pt.plot_normal(time_days, np.linalg.norm(UJ_pos - U_JG_rel_pos, axis=1), 'Time [days]', r'Absolute Position Norm [m]',
    #                 "Absolute Position of JUICE in Jupiter-centered Unperturbed Simulation \n - Relative Position Ganymede-Jupiter", 
    #                 'abs_pos_UJ-JG.pdf', None, 'question4', lwidth=1  )



    # # # --- 3D ORBIT PLOT ---
    # # # Load the raw Jupiter-centered state of JUICE (columns 1, 2, 3 correspond to x, y, z)
    # # UJ_state_Jcentered = np.loadtxt(f'./Data/question4/JUICEPropagationHistory_Q4_Unperturbed_Jupiter_centered.dat')
    # # UJ_pos_Jcentered = UJ_state_Jcentered[:, 1:4] 
    
    # # # U_JG_rel_pos is already loaded in your script (Ganymede's state wrt Jupiter)
    
    # # print("Generating 3D orbit comparison plot...")
    # # pt.plot_3d_orbits(
    # #     pos1=UJ_pos_Jcentered, label1="JUICE", 
    # #     pos2=U_JG_rel_pos, label2="Ganymede", 
    # #     title="3D Trajectories of JUICE and Ganymede \n (Unperturbed, Jupiter-Centered)", 
    # #     filename='3d_orbits_UJ.pdf', 
    # #     task='question4'
    # # )



    # # PJ_state_Jcentered = np.loadtxt(f'./Data/question4/JUICEPropagationHistory_Q4_Perturbed_Jupiter_centered.dat')
    # # PJ_pos_Jcentered = PJ_state_Jcentered[:, 1:4] 
    
    # # # U_JG_rel_pos is already loaded in your script (Ganymede's state wrt Jupiter)
    
    # # print("Generating 3D orbit comparison plot...")
    # # pt.plot_3d_orbits(
    # #     pos1=PJ_pos_Jcentered, label1="JUICE", 
    # #     pos2=P_JG_rel_pos, label2="Ganymede", 
    # #     title="3D Trajectories of JUICE and Ganymede \n (Perturbed, Jupiter-Centered)", 
    # #     filename='3d_orbits_PJ.pdf', 
    # #     task='question4'
    # # )


    # # --- 3D ORBIT PLOTS ---
    
    # # Load the raw Jupiter-centered state of JUICE (columns 1, 2, 3 correspond to x, y, z)
    # UJ_state_Jcentered = np.loadtxt(f'./Data/question4/JUICEPropagationHistory_Q4_Unperturbed_Jupiter_centered.dat')
    # UJ_pos_Jcentered = UJ_state_Jcentered[:, 1:4] 
    
    # print("Generating 3D orbit plots...")

    # # 1. Jupiter-Centered (Unperturbed)
    # # Shows Ganymede orbiting Jupiter while JUICE flies off into space
    # pt.plot_3d_orbits(
    #     trajectories={"JUICE": UJ_pos_Jcentered, "Ganymede": U_JG_rel_pos}, 
    #     origin_name="Jupiter",
    #     title="3D Trajectories of JUICE and Ganymede \n (Unperturbed, Jupiter-Centered)", 
    #     filename='3d_orbits_UJ_Jupiter_Centered.pdf', 
    #     task='question4'
    # )

    # # 2. Ganymede-Centered (Unperturbed)
    # # Shows the massive 7,000 km distortion of JUICE's orbit
    # pt.plot_3d_orbits(
    #     trajectories={"JUICE": UJ_pos}, # UJ_pos is already Ganymede-centered from your previous subtraction
    #     origin_name="Ganymede",
    #     title="3D Trajectory of JUICE \n (Unperturbed, Ganymede-Centered)", 
    #     filename='3d_orbits_UJ_Ganymede_Centered.pdf', 
    #     task='question4'
    # )
    
    # # 3. Ganymede-Centered (Perturbed)
    # # Shows the stable, correct orbit for comparison!
    # pt.plot_3d_orbits(
    #     trajectories={"JUICE": PJ_pos}, # PJ_pos is already Ganymede-centered
    #     origin_name="Ganymede",
    #     title="3D Trajectory of JUICE \n (Perturbed, Ganymede-Centered)", 
    #     filename='3d_orbits_PJ_Ganymede_Centered.pdf', 
    #     task='question4'
    # )


    # pt.plot_3d_orbits(
    #     trajectories={"JUICE": UJ_pos[:10000,:]-UG_pos[:10000, :]}, # PJ_pos is already Ganymede-centered
    #     origin_name="Ganymede",
    #     title="3D Trajectory error of JUICE \n (Unperturbed, (Jupiter-Ganymede)-Centered)", 
    #     filename='3d_orbits_UJ-UG_Ganymede_Centered.pdf', 
    #     task='question4'
    # )

    # UJ_G_pos = np.loadtxt(f'./Data/question4/JUICEPropagationHistory_Q4_Unperturbed_Jupiter_centered_Ganymede.dat')[:, 1:4]
    # UJ_G_G_pos = np.loadtxt(f'./Data/question4/JUICEPropagationHistory_DependentVariables_Q4_Unperturbed_Jupiter_centered_Ganymede.dat')[:, 7:10]   
    # UJ_G_pos_Gcentered = UJ_G_pos - UJ_G_G_pos  # Subtract Ganymede's position to get JUICE's position in Ganymede-centered frame
    # pt.plot_3d_orbits(
    #     trajectories={"JUICE": UJ_G_pos_Gcentered[:10000,:]-UG_pos[:10000, :]}, # PJ_pos is already Ganymede-centered
    #     origin_name="Ganymede",
    #     title="3D Trajectory error of JUICE \n (Unperturbed, (Jupiter-Ganymede)-Centered)", 
    #     filename='3d_orbits_UJ_G-UG_Ganymede_Centered.pdf', 
    #     task='question4'
    # )





import os
import numpy as np
from matplotlib import pyplot as plt

from tudatpy import constants
from tudatpy.data import save2txt
from tudatpy.interface import spice
from tudatpy.dynamics import environment_setup, propagation_setup, simulator
from tudatpy.dynamics.environment_setup import rotation_model

current_directory = os.getcwd()

A = 7
B = 5
C = 3

simulation_start_epoch = (
    35.4 * constants.JULIAN_YEAR
    + A * 7.0 * constants.JULIAN_DAY
    + B * constants.JULIAN_DAY
    + C * constants.JULIAN_DAY / 24.0
)
simulation_end_epoch = simulation_start_epoch + 344.0 * constants.JULIAN_DAY / 24.0

fixed_step_size = 10.0

###########################################################################
# SIMULATION FUNCTION SETUP ###############################################
###########################################################################

def run_q4_simulation(model="Q1", center="Ganymede", concurrent=False):
    """
    Funzione per eseguire i diversi casi richiesti dalla Q4.
    """
    spice.load_standard_kernels()
    spice.load_kernel(os.path.join(current_directory, "juice_mat_crema_5_1_150lb_v01.bsp"))

    bodies_to_create = ["Jupiter", "Ganymede", "Sun", "Saturn", "Io", "Europa", "Callisto"]
    global_frame_origin = center
    global_frame_orientation = "ECLIPJ2000"

    body_settings = environment_setup.get_default_body_settings(
        bodies_to_create, global_frame_origin, global_frame_orientation
    )

    if model == "Q2":
        rho0 = 2.0e-9
        H = 40.0e3
        body_settings.get("Ganymede").atmosphere_settings = environment_setup.atmosphere.exponential(H, rho0)
        
        L_sun = 3.828e26
        luminosity_settings = environment_setup.radiation_pressure.constant_luminosity(L_sun)
        body_settings.get("Sun").radiation_source_settings = environment_setup.radiation_pressure.isotropic_radiation_source(luminosity_settings)

    body_settings.add_empty_settings("JUICE")
    body_settings.get("JUICE").constant_mass = 2.0e3

    bodies = environment_setup.create_system_of_bodies(body_settings)

    if model == "Q2":
        reference_area = 100.0
        CD = 1.2
        Cr = 1.2
        aero_coefficient_settings = environment_setup.aerodynamic_coefficients.constant(reference_area, [CD, 0.0, 0.0])
        environment_setup.add_aerodynamic_coefficient_interface(bodies, "JUICE", aero_coefficient_settings)
        
        occulting_bodies_dict = {"Sun": ["Ganymede"]}
        juice_target_settings = environment_setup.radiation_pressure.cannonball_radiation_target(
            reference_area=reference_area, radiation_pressure_coefficient=Cr, per_source_occulting_bodies=occulting_bodies_dict
        )
        environment_setup.add_radiation_pressure_target_model(bodies, "JUICE", juice_target_settings)

    bodies_to_propagate = ["JUICE"]
    central_bodies = [center]

    if concurrent:
        bodies_to_propagate.append("Ganymede")
        central_bodies.append("Jupiter")

    acc_dict = {"JUICE": {}}
    if concurrent:
        acc_dict["Ganymede"] = {}

    if model == "Q1":
        acc_dict["JUICE"]["Ganymede"] = [propagation_setup.acceleration.point_mass_gravity()]
        
    elif model == "Q2":
        acc_dict["JUICE"]["Ganymede"] = [propagation_setup.acceleration.spherical_harmonic_gravity(2, 2), propagation_setup.acceleration.aerodynamic()]
        acc_dict["JUICE"]["Jupiter"] = [propagation_setup.acceleration.spherical_harmonic_gravity(4, 0)]
        acc_dict["JUICE"]["Sun"] = [propagation_setup.acceleration.point_mass_gravity(), propagation_setup.acceleration.radiation_pressure()]
        acc_dict["JUICE"]["Saturn"] = [propagation_setup.acceleration.point_mass_gravity()]
        acc_dict["JUICE"]["Io"] = [propagation_setup.acceleration.point_mass_gravity()]
        acc_dict["JUICE"]["Europa"] = [propagation_setup.acceleration.point_mass_gravity()]
        acc_dict["JUICE"]["Callisto"] = [propagation_setup.acceleration.point_mass_gravity()]
        
        if concurrent:
            acc_dict["Ganymede"]["Jupiter"] = [propagation_setup.acceleration.spherical_harmonic_gravity(4, 0)]
            acc_dict["Ganymede"]["Sun"] = [propagation_setup.acceleration.point_mass_gravity()]
            acc_dict["Ganymede"]["Saturn"] = [propagation_setup.acceleration.point_mass_gravity()]
            acc_dict["Ganymede"]["Io"] = [propagation_setup.acceleration.point_mass_gravity()]
            acc_dict["Ganymede"]["Europa"] = [propagation_setup.acceleration.point_mass_gravity()]
            acc_dict["Ganymede"]["Callisto"] = [propagation_setup.acceleration.point_mass_gravity()]

    acceleration_models = propagation_setup.create_acceleration_models(
        bodies, acc_dict, bodies_to_propagate, central_bodies
    )

    initial_state = spice.get_body_cartesian_state_at_epoch(
        "JUICE", center, global_frame_orientation, "NONE", simulation_start_epoch
    )
    
    if concurrent:
        initial_state_ganymede = spice.get_body_cartesian_state_at_epoch(
            "Ganymede", "Jupiter", global_frame_orientation, "NONE", simulation_start_epoch
        )
        initial_state = np.concatenate((initial_state, initial_state_ganymede))

    integrator_settings = propagation_setup.integrator.runge_kutta_fixed_step(
        fixed_step_size, propagation_setup.integrator.CoefficientSets.rk_4
    )
    termination_settings = propagation_setup.propagator.time_termination(simulation_end_epoch)

    propagator_settings = propagation_setup.propagator.translational(
        central_bodies, acceleration_models, bodies_to_propagate, initial_state,
        simulation_start_epoch, integrator_settings, termination_settings
    )

    dyn = simulator.create_dynamics_simulator(bodies, propagator_settings)
    prop_results = dyn.propagation_results

    # Data Extraction
    times = np.array(list(prop_results.state_history.keys()))
    state_history_array = np.vstack(list(prop_results.state_history.values()))
    
    if center == "Ganymede":
        states_wrt_ganymede = state_history_array[:, 0:6]
        
    elif center == "Jupiter":
        states_wrt_ganymede = np.zeros_like(state_history_array[:, 0:6])
        if concurrent:
            x_J_s = state_history_array[:, 0:6]
            x_J_G_prop = state_history_array[:, 6:12]
            states_wrt_ganymede = x_J_s - x_J_G_prop
        else:
            for i in range(len(times)):
                x_J_G_spice = spice.get_body_cartesian_state_at_epoch(
                    "Ganymede", "Jupiter", global_frame_orientation, "NONE", times[i])
                states_wrt_ganymede[i, :] = state_history_array[i, 0:6] - x_J_G_spice

    return times, states_wrt_ganymede, prop_results

###########################################################################
# PROPAGATION RUNNING (in all cases)#######################################
###########################################################################

t_i, s_i, _ = run_q4_simulation("Q1", "Ganymede", concurrent=False)
t_iii, s_iii, _ = run_q4_simulation("Q2", "Ganymede", concurrent=False)

t_ii_bad, s_ii_bad, _ = run_q4_simulation("Q1", "Jupiter", concurrent=False)
t_iv_bad, s_iv_bad, _ = run_q4_simulation("Q2", "Jupiter", concurrent=False)

t_ii_good, s_ii_good, _ = run_q4_simulation("Q1", "Jupiter", concurrent=True)
t_iv_good, s_iv_good, results_q4_nominal = run_q4_simulation("Q2", "Jupiter", concurrent=True)

###########################################################################
# DATA EXTRACTION AND ERRORS CALCULATION ##################################
###########################################################################

time_hours = (t_i - t_i[0]) / 3600.0

# Errors in i-ii case (first formulation)
err_q1_bad = np.linalg.norm(s_i[:, 0:3] - s_ii_bad[:, 0:3], axis=1)
err_q2_bad = np.linalg.norm(s_iii[:, 0:3] - s_iv_bad[:, 0:3], axis=1)

# Errors in iii-iv case (second formulation)
err_q1_good = np.linalg.norm(s_i[:, 0:3] - s_ii_good[:, 0:3], axis=1)
err_q2_good = np.linalg.norm(s_iii[:, 0:3] - s_iv_good[:, 0:3], axis=1)







































































# def plot_mixed_subplots(time, 
#                          data1_top, label1_top, ylabel1_top, 
#                          data2_top, label2_top, ylabel2_top,
#                          data_bot, label_bot, ylabel_bot, 
#                          lwidth_top=6,
#                          lwidth_bot=0.5,
#                          title="Generic Title", filename="mixed_subplots.pdf", 
#                          task="default"):
#     """
#     Creates a single figure with 2 vertical subplots. 
#     TOP: Two Y-axes (twinx) for comparing two scales.
#     BOTTOM: Normal single Y-axis plot.
#     """
#     task_dir = Path(f"Plots/{task}")
#     task_dir.mkdir(exist_ok=True, parents=True)
#     filepath = task_dir / filename

#     # Create figure with 2 rows, 1 column
#     fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(8, 10))
#     # hspace adds room between plots; top leaves room for suptitle/legend
#     fig.subplots_adjust(hspace=0.2, top=0.88) 

#     color1, color2 = 'tab:blue', 'tab:red'

#     # --- TOP SUBPLOT (Two Scales) ---
#     # Axis 1 (Left)
#     ax_top.plot(time, data1_top, color=color1, label=label1_top, linewidth=lwidth_top-lwidth_top/3, zorder=1)
#     ax_top.set_ylabel(ylabel1_top, color=color1, fontsize=LABEL_SIZE)
#     ax_top.tick_params(axis='y', labelcolor=color1, labelsize=TICK_SIZE)
#     ax_top.grid(True, linestyle=':', alpha=0.6)
    
#     # Axis 2 (Right)
#     ax_top_twin = ax_top.twinx()
#     ax_top_twin.plot(time, data2_top, color=color2, label=label2_top, linestyle=':', linewidth=lwidth_top)
#     ax_top_twin.set_ylabel(ylabel2_top, color=color2, fontsize=LABEL_SIZE)
#     ax_top_twin.tick_params(axis='y', labelcolor=color2, labelsize=TICK_SIZE)
    
#     # Legend for Top (Combined from both axes)
#     lines_t1, labels_t1 = ax_top.get_legend_handles_labels()
#     lines_t2, labels_t2 = ax_top_twin.get_legend_handles_labels()
#     ax_top.legend(lines_t1 + lines_t2, labels_t1 + labels_t2, loc='lower center', 
#                   bbox_to_anchor=(0.5, 1.02), ncol=2, fontsize=LEGEND_SIZE, frameon=True)

#     # --- BOTTOM SUBPLOT (Normal Plot) ---
#     ax_bot.plot(time, data_bot, color='blue', label=label_bot, linewidth=lwidth_bot)
#     ax_bot.set_xlabel("Time [days]", fontsize=LABEL_SIZE)
#     ax_bot.set_ylabel(ylabel_bot, fontsize=LABEL_SIZE)
#     ax_bot.tick_params(axis='both', labelsize=TICK_SIZE)
#     ax_bot.grid(True, linestyle=':', alpha=0.6)
    
#     # Legend for Bottom
#     ax_bot.legend(loc='lower center', bbox_to_anchor=(0.5, 1.02), 
#                   fontsize=LEGEND_SIZE, frameon=True)

#     # Global Title
#     fig.suptitle(title, fontsize=TITLE_SIZE, y=0.98, fontweight='bold')

#     # Save logic
#     plt.savefig(filepath, bbox_inches="tight", dpi=300)
#     plt.close(fig)
#     return str(filepath)


# def plot_three_vertical_subplots(time, 
#                                 data_top, ylabel_top,
#                                 data_mid, ylabel_mid,
#                                 data_bot, ylabel_bot, 
#                                 lwidth=1.5,
#                                 title="Generic Title", filename="three_compact_subplots.pdf", 
#                                 task="default"):
#     """
#     Creates a compact figure with 3 vertical subplots sharing the same X-axis.
#     No legends, shared x-label at the bottom.
#     """
#     task_dir = Path(f"Plots/{task}")
#     task_dir.mkdir(exist_ok=True, parents=True)
#     filepath = task_dir / filename

#     # Create figure with 3 rows, shared X axis
#     fig, (ax_top, ax_mid, ax_bot) = plt.subplots(3, 1, figsize=(8, 10), sharex=True)
    
#     # hspace=0.05 makes them very compact; top margin for the main title
#     fig.subplots_adjust(hspace=0.1, top=0.92, bottom=0.1) 

#     # --- TOP SUBPLOT ---
#     ax_top.plot(time, data_top, color='tab:blue', linewidth=lwidth)
#     ax_top.set_ylabel(ylabel_top, fontsize=LABEL_SIZE)
#     ax_top.tick_params(axis='both', labelsize=TICK_SIZE)
#     ax_top.grid(True, linestyle=':', alpha=0.6)

#     # --- MIDDLE SUBPLOT ---
#     ax_mid.plot(time, data_mid, color='tab:orange', linewidth=lwidth)
#     ax_mid.set_ylabel(ylabel_mid, fontsize=LABEL_SIZE)
#     ax_mid.tick_params(axis='both', labelsize=TICK_SIZE)
#     ax_mid.grid(True, linestyle=':', alpha=0.6)

#     # --- BOTTOM SUBPLOT ---
#     ax_bot.plot(time, data_bot, color='tab:green', linewidth=lwidth)
#     ax_bot.set_ylabel(ylabel_bot, fontsize=LABEL_SIZE)
#     ax_bot.set_xlabel("Time [days]", fontsize=LABEL_SIZE)
#     ax_bot.tick_params(axis='both', labelsize=TICK_SIZE)
#     ax_bot.grid(True, linestyle=':', alpha=0.6)

#     # Global Title
#     fig.suptitle(title, fontsize=TITLE_SIZE, y=0.97, fontweight='bold')

#     # Save logic
#     plt.savefig(filepath, bbox_inches="tight", dpi=300)
#     plt.close(fig)
#     return str(filepath)