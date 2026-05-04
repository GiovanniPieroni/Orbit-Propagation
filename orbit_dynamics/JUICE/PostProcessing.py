import numpy as np
import plotting as pt
from debug import pick

import subprocess
import sys

############################################################################################################################################
#                                                           Folders setup                                                                  # 
############################################################################################################################################

current_task = pick('task')
only_plotting = pick('only_plotting')

# Define the shared file name
filename = "cartesian results AE4868 2026 A1 6541151.txt"




all_scripts = [
    "question1.py",       
    "question2.py",
    "question3.py",
    "question4.py"
]

# Decide which task to execute
if current_task == "all":
    scripts_to_run = all_scripts
    print("All scripts selected to be executed.")
elif f"{current_task}.py" in all_scripts:
    scripts_to_run = [f"{current_task}.py"]
    print(f"Executing only the requested task: {current_task}.py")
else:
    print(f"Warning: no task found with name '{current_task}'.")
    scripts_to_run = []


# 4. Executing selected scripts - if not only plotting
if not only_plotting:
    # Initialize the file 
    print(f"Initializing {filename}...")
    with open(filename, 'w') as f:
        f.write('')

    for script in scripts_to_run:

        

        print(f"\n--- Running {script} ---")
        
        
        result = subprocess.run([sys.executable, script])
        
        # Check for errors
        if result.returncode != 0:
            print(f"WARNING: {script} encountered an error.")
        else:
            print(f"SUCCESS: {script} finished and appended its data.")

    print("\nAll simulations completed.")
else:
    print('\n Only plotting selected, using previous simulations data.\n')





############################################################################################################################################
#                                                           QUESTION 1                                                                     # 
############################################################################################################################################
if current_task == 'question1' or current_task == 'all':

    # Settings
    SC = pick('SpaceCraft', override_task='question1')
    CB = pick('CentralBody', override_task='question1')
    gModel = pick('GravityModel', override_task='question1')


    # Keplerian elements
    time_days = np.loadtxt(f'./Data/question1/epochs_days.txt')
    kep_state_num = np.loadtxt(f'./Data/question1/JUICEPropagationHistory_DependentVariables_Q1.txt')

    pt.plot_kep_state_separate(time_days, kep_state_num[:, 1:]-kep_state_num[0, 1:], SC, CB, gModel, 'question1')

    kep_state_num[:, 4] = np.unwrap(kep_state_num[:, 4])
    kep_state_num[:, 6] = np.unwrap(kep_state_num[:, 6])

    print(f'Starting to plot data for question 1')

    pt.plot_keplerian_state(time_days, kep_state_num[:, 1:]-kep_state_num[0, 1:], 'question1', 'kep_state.pdf')


    # Residuals between analytical and numerical solution
    kep_state_residual_PM = np.loadtxt(f'./Data/question1/keplerian_residual.txt')
    car_state_residual_PM = np.loadtxt(f'./Data/question1/cartesian_residual.txt')



    pt.plot_kep_state_error_separate(time_days, kep_state_residual_PM[:, 1:], SC, CB, gModel, None, 'question1')


    pt.plot_keplerian_state_error(time_days, kep_state_residual_PM[:, 1:], 'question1', 'kep_state_error.pdf', 
                                  'Residuals of Keplerian Elements between Numerical and Analytical Solution')

    pt.plot_car_state_error(time_days, car_state_residual_PM[:, 1:], SC, CB, gModel, None, 'question1')

    print(f'Completed to plot data for question 1')






############################################################################################################################################
#                                                           QUESTION 2                                                                     # 
############################################################################################################################################
if current_task == 'question2' or current_task == 'all':

    # Settings
    SC = pick('SpaceCraft', override_task='question2')
    CB = pick('CentralBody', override_task='question2')
    gModel = pick('GravityModel', override_task='question2')

    # (i)
    Io_grav_acc = np.loadtxt(f'./Data/question2/JUICEPropagationHistory_DependentVariables_Q2.txt')[:, 10:13]
    Sun_grav_acc = np.loadtxt(f'./Data/question2/JUICEPropagationHistory_DependentVariables_Q2.txt')[:, 13:16]
    
    time_days = np.loadtxt(f'./Data/question2/epochs_days.txt')

    print(f'Starting to plot data for question 2')

    
    pt.plot_acceleration_norm(time_days, Io_grav_acc, SC, 'Io', CB, gModel, None, 'question2')
    pt.plot_acceleration_norm(time_days, Sun_grav_acc, SC, 'Sun', CB, gModel, None, 'question2')

    print(f'Completed to plot data for question 2')
    


    # (ii)
    Io_freq_spctrum = np.loadtxt(f'./Data/question2/Io_frequency_spectrum.txt')
    Sun_freq_spectrum = np.loadtxt(f'./Data/question2/Sun_frequency_spectrum.txt')

    pt.plot_normal(1/Io_freq_spctrum[1:, 0]/3600, Io_freq_spctrum[1:, 1],'Period [hr]',r'Power Density $[m^2/s^4]$',
                    f"Io's Acceleration Frequency Spectrum", 'Io_freq_spect.pdf',None, 'question2', xlim=(0, 100))

    pt.plot_normal(1/Sun_freq_spectrum[1:, 0]/3600, Sun_freq_spectrum[1:, 1],'Period [hr]',r'Power Density $[m^2/s^4]$',
                    f"Sun's Acceleration Frequency Spectrum", 'Sun_freq_spect.pdf',None, 'question2', xlim=(0, 5))






############################################################################################################################################
#                                                           QUESTION 3                                                                     # 
############################################################################################################################################
if current_task == 'question3' or current_task == 'all':

    # Settings
    SC = pick('SpaceCraft', override_task='question3')
    CB = pick('CentralBody', override_task='question3')
    gModel = pick('GravityModel', override_task='question3')

# (a)
    # Aerodynamic acceleration numerical and analytical comparison
    aer_acc_num = np.loadtxt(f'./Data/question3/JUICEPropagationHistory_DependentVariables_Q3_complete.txt')[:, 7:10]
    aer_acc_analyt = np.loadtxt(f'./Data/question3/acc_aer_analyt.txt')
    aer_acc_num_mag = np.linalg.norm(aer_acc_num, axis=1)
    aer_acc_analyt_mag = np.linalg.norm(aer_acc_analyt, axis=1)

    time_days = np.loadtxt(f'./Data/question3/epochs_days.txt')

    print(f'Starting to plot data for question 3')

    pt.plot_mixed_subplots(
                            time_days, 
                            # Top: Two magnitudes to see if they track each other
                            aer_acc_num_mag, 'Numerical', r'$a [m/s^2]$',
                            aer_acc_analyt_mag, 'Manually computed', r'$a [m/s^2]$',
                            # Bottom: The absolute error or relative difference
                            (aer_acc_num_mag - aer_acc_analyt_mag), 'Residual', r'$\Delta a [m/s^2]$',
                            lwidth_top=2, lwidth_bot=0.5,
                            # Meta
                            title="Atmospheric Drag Acceleration Magnitude \n Comparison and Residuals",
                            filename="aer_acc_comparison.pdf",
                            task="question3")
    


    aer_acc_analyt_dir = aer_acc_analyt / np.linalg.norm(aer_acc_analyt, axis=1)[:, None]
    aer_acc_num_dir = aer_acc_num / np.linalg.norm(aer_acc_num, axis=1)[:, None]
    aer_acc_dir_residual = aer_acc_analyt_dir - aer_acc_num_dir




    pt.plot_three_vertical_subplots(time_days, 
                                aer_acc_dir_residual[:, 0], r'$\Delta a_x [m/s^2]$',
                                aer_acc_dir_residual[:, 1], r'$\Delta a_y [m/s^2]$',
                                aer_acc_dir_residual[:, 2], r'$\Delta a_z [m/s^2]$',
                                lwidth=0.5,
                                title="Atmospheric Drag Acceleration Direction Residuals", filename="aer_acc_direction_residuals.pdf", 
                                task="question3")
    


    # SRP acceleration numerical and analytical comparison
    srp_acc_num = np.loadtxt(f'./Data/question3/JUICEPropagationHistory_DependentVariables_Q3_complete.txt')[:, 10:13]
    srp_acc_analyt = np.loadtxt(f'./Data/question3/acc_rad_sun_analyt.txt')
    srp_acc_num_mag = np.linalg.norm(srp_acc_num, axis=1)
    srp_acc_analyt_mag = np.linalg.norm(srp_acc_analyt, axis=1) 



    pt.plot_mixed_subplots(
                            time_days, 
                            # Top: Two magnitudes to see if they track each other
                            srp_acc_num_mag, 'Numerical', r'$a [m/s^2]$',
                            srp_acc_analyt_mag, 'Manually computed', r'$a [m/s^2]$',
                            # Bottom: The absolute error or relative difference
                            (srp_acc_num_mag - srp_acc_analyt_mag), 'Residual', r'$\Delta a [m/s^2]$',
                            lwidth_top=6, lwidth_bot=0.5,
                            title="Solar Radiation Pressure Acceleration Magnitude \n Comparison and Residuals",
                            filename="srp_acc_comparison.pdf",
                            task="question3")



    srp_acc_analyt_dir = srp_acc_analyt / np.linalg.norm(srp_acc_analyt, axis=1)[:, None]
    srp_acc_num_dir = srp_acc_num / np.linalg.norm(srp_acc_num, axis=1)[:, None]
    srp_acc_dir_residual = srp_acc_analyt_dir - srp_acc_num_dir


    pt.plot_three_vertical_subplots(time_days, 
                                srp_acc_dir_residual[:, 0], r'$\Delta a_x [m/s^2]$',
                                srp_acc_dir_residual[:, 1], r'$\Delta a_y [m/s^2]$',
                                srp_acc_dir_residual[:, 2], r'$\Delta a_z [m/s^2]$',
                                lwidth=0.5,
                                title="Solar Radiation Pressure Acceleration Direction Residuals", filename="srp_acc_direction_residuals.pdf", 
                                task="question3")


# (b)
    acc_aer_rsw = np.loadtxt(f'./Data/question3/acc_aer_rsw.txt')

    pt.plot_three_vertical_subplots(time_days, 
                                acc_aer_rsw[:, 0], r'$a_R [m/s^2]$',
                                acc_aer_rsw[:, 1], r'$a_S [m/s^2]$',
                                acc_aer_rsw[:, 2], r'$a_W [m/s^2]$',
                                lwidth=1.5,
                                title="Atmospheric Drag Acceleration in RSW frame", filename="aer_acc_RSW.pdf", 
                                task="question3")


    acc_srp_rsw = np.loadtxt(f'./Data/question3/acc_srp_rsw.txt')

    pt.plot_three_vertical_subplots(time_days, 
                                acc_srp_rsw[:, 0], r'$a_R [m/s^2]$',
                                acc_srp_rsw[:, 1], r'$a_S [m/s^2]$',
                                acc_srp_rsw[:, 2], r'$a_W [m/s^2]$',
                                lwidth=1.5,
                                title="Solar Radiation Pressure Acceleration in RSW frame", filename="srp_acc_RSW.pdf", 
                                task="question3")


# (c) 
    
    no_drag_cartesian_state = np.loadtxt(f'./Data/question3/no_drag_cartesian_state.txt')
    no_srp_cartesian_state = np.loadtxt(f'./Data/question3/no_srp_cartesian_state.txt')
    complete_cartesian_state = np.loadtxt(f'./Data/question3/complete_cartesian_state.txt') 

    no_drag_pos_residual = complete_cartesian_state[:, :3] - no_drag_cartesian_state[:, :3]
    no_srp_pos_residual = complete_cartesian_state[:, :3] - no_srp_cartesian_state[:, :3]   

    pt.plot_normal(time_days, np.linalg.norm(no_drag_pos_residual, axis=1), 
                   'Time [days]', r'Position Residual Norm [m]',
                    "Position Residuals between Complete and No-Drag Simulations", 
                    'pos_residual_no_drag.pdf', None, 'question3', lwidth=2)                           

    pt.plot_normal(time_days, np.linalg.norm(no_srp_pos_residual, axis=1), 
                   'Time [days]', r'Position Residual Norm [m]',
                    "Position Residuals between Complete and No-SRP Simulations", 
                    'pos_residual_no_srp.pdf', None, 'question3', lwidth=2)                           
    

    print(f'Completed to plot data for question 3')






############################################################################################################################################
#                                                           QUESTION 4                                                                     # 
############################################################################################################################################
if current_task == 'question4' or current_task == 'all':

    # Settings
    SC = pick('SpaceCraft', override_task='question4')
    CB = pick('CentralBody', override_task='question4')

    time_days = np.loadtxt(f'./Data/question4/epochs_days.txt')


    print(f'Starting to plot data for question 4')

# (0) Ganymede not propagated

    UJ_pos = np.loadtxt(f'./Data/question4/UJ_cartesian_state_Gcentered.txt')[:, :3]
    PJ_pos = np.loadtxt(f'./Data/question4/PJ_cartesian_state_Gcentered.txt')[:, :3]
   
    UG_pos = np.loadtxt(f'./Data/question1/JUICEPropagationHistory_Q1.txt')[:, 1:4]
    PG_pos = np.loadtxt(f'./Data/question2/JUICEPropagationHistory_Q2.txt')[:, 1:4]


    pos_difference_UG_UJ = np.linalg.norm(UJ_pos - UG_pos, axis=1)
    pos_difference_PG_PJ = np.linalg.norm(PJ_pos - PG_pos, axis=1)

    pt.plot_normal(time_days, pos_difference_UG_UJ, 'Time [days]', r'$||\Delta r(t)||_{i, ii}$ [m]',
                    "Position Difference between Ganymede-centered and Jupiter-centered Unperturbed Simulations", 
                    'pos_diff_UG_UJ.pdf', None, 'question4', lwidth=2)

    pt.plot_normal(time_days, pos_difference_PG_PJ, 'Time [days]', r'$||\Delta r(t)||_{iii, iv}$ [m]',
                    "Position Difference between Ganymede-centered and Jupiter-centered Perturbed Simulations", 
                    'pos_diff_PG_PJ.pdf', None, 'question4', lwidth=2)
    

# (d) Ganymede propagated:

    UJ_Gpropagated_pos = np.loadtxt(f'./Data/question4/UJ_G_cartesian_state_Gcentered.txt')[:, :3]
    PJ_Gpropagated_pos = np.loadtxt(f'./Data/question4/PJ_G_cartesian_state_Gcentered.txt')[:, :3]

    pos_difference_UJ_G = np.linalg.norm(UJ_Gpropagated_pos - UG_pos, axis=1)
    pos_difference_PJ_G = np.linalg.norm(PJ_Gpropagated_pos - PG_pos, axis=1)

    pt.plot_normal(time_days, pos_difference_UJ_G, 'Time [days]', r'$||\Delta r(t)||_{i, ii}$ [m]',
                    "Position Difference between Ganymede- and Jupiter-centered Unperturbed Simulations \n with Ganymede's Propagated State", 
                    'pos_diff_UJ_G.pdf', None, 'question4', lwidth=2)

    pt.plot_normal(time_days, pos_difference_PJ_G, 'Time [days]', r'$||\Delta r(t)||_{iii, iv}$ [m]',
                    "Position Difference between Ganymede- and Jupiter-centered Perturbed Simulations \n with Ganymede's Propagated State", 
                    'pos_diff_PJ_G.pdf', None, 'question4', lwidth=2)


    pt.plot_normal(time_days, pos_difference_UG_UJ-pos_difference_UJ_G, 'Time [days]', r'Position Residual Norm [m]',
                    "Position Residuals between Ganymede-centered and Jupiter-centered Unperturbed Simulations \n with and without Ganymede's Propagated State", 
                    'pos_residual_UG_UJ_vs_UJ_G.pdf', None, 'question4', lwidth=2)  


    print(f'Completed to plot data for question 4')


