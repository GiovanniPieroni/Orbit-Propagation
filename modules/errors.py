from __future__ import annotations
import numpy as np
from debug import pick
from frame_transformation import cartesian_to_keplerian
from plots import plot_three_scales


# ----------------------------------------------------- A4 ----------------------------------------------------------- #
def compute_gpd(r_num: np.ndarray, r_sgp4: np.ndarray) -> np.ndarray:
    """
    Computes the Global Position Difference (GPD) between two orbits.

    Parameters:
        r_num  : Nx3 array of positions from numerical integrator (Euler, RK4, etc.)
        r_sgp4 : Nx3 array of positions from SGP4 propagation

    Returns:
        gpd : Nx1 array of L2 norms of the position residual at each epoch
    """

    if r_num.shape[0] != r_sgp4.shape[0]:
        raise RuntimeError("The two vectors have different rows!\n")
    else:
        if r_num.shape[1] != r_sgp4.shape[1]:
            raise RuntimeError("The two vectors have different columns!\n")

    residual = r_num - r_sgp4  # Nx3
    gpd = np.linalg.norm(residual, axis=1)  # Nx
    return gpd



from plots import mask_discontinuity, plot_comparison_residual
# import matplotlib.pyplot as plt



def kep_el_and_residuals(state_u, time_u,
                         state_d, time_d,
                         label_u, ylabel_u,
                         label_d, ylabel_d,
                         label_diff, ylabel_diff,
                         fig_title, fig_name, plot=True):
    # Retrieve Keplerian elements
    if not np.array_equal(time_u, time_d):
        raise RuntimeError("The two orbits have different time domains, residuals cannot be computed.")
    a_u, e_u, _, i_u, _, Om_u, _, om_u, _, th_u, u_u = cartesian_to_keplerian(state_u)
    a_d, e_d, _, i_d, _, Om_d, _, om_d, _, th_d, u_d = cartesian_to_keplerian(state_d)
    da = a_d - a_u
    de = e_d - e_u
    di = i_d - i_u
    dOm = Om_d - Om_u
    dom = om_d - om_u
    dth = th_d - th_u
    du = u_d - u_u
    # kep_el_u = np.column_stack((a_u, e_u, i_u, Om_u, om_u, th_u, u_u))
    # kep_el_d = np.column_stack((a_d, e_d, i_d, Om_d, om_d, th_d, u_d))
    #
    # kep_el_diff = np.column_stack((da, de, di, dOm, dom, dth, du))
    kep_el_u = np.column_stack((a_u, e_u, i_u, Om_u, u_u))
    kep_el_d = np.column_stack((a_d, e_d, i_d, Om_d, u_d))
    kep_el_diff = np.column_stack((da, de, di, dOm, du))
    generated_files = []

    if plot:
        for i in range(5):
            # Extract data
            u_data = kep_el_u[:, i]
            d_data = kep_el_d[:, i]
            diff_data = kep_el_diff[:, i]

            # Handle Angle Wrapping for Omega (idx 3) and omega (idx 4)
            if i in [3, 4]:
                u_data = mask_discontinuity(u_data)
                d_data = mask_discontinuity(d_data)
                # diff data usually doesn't need masking if it's small residuals,
                # but if it jumps due to wrapping, we might need it.
                # For residuals, we usually don't wrap them to 0-360, they are small deltas.
                # diff_data = mask_discontinuity(diff_data)

            # Handle filename list or single pattern
            # If 'filename' is a list, pick i-th. If string, append _i
            fname = fig_name[i] if isinstance(fig_name, list) else f"{fig_name}_{i}.png"
            tit = fig_title[i] if isinstance(fig_title, list) else f"{fig_title} - Element {i}"

            # Call the plotter
            out = plot_comparison_residual(
                time_u,
                u_data, label_u[i],
                d_data, label_d[i],
                diff_data, label_diff[i],
                ylabel_u[i], ylabel_diff[i],  # Main Y-label matches the elements
                title=tit,
                filename=fname
            )
            # generated_files.append(out)

    return kep_el_u, kep_el_d, kep_el_diff


def mag_dist(acc_v: np.ndarray):
    """
    :param acc_v: Time series of acceleration
    :return: magnitude of acceleration vector
    """
    # ensure exactly one dimension is 3
    if 3 not in acc_v.shape:
        raise ValueError("Input must have one dimension equal to 3")

    # orient as (N, 3)
    if acc_v.shape[1] == 3:
        v2 = acc_v
    else:  # shape must be (3, N)
        v2 = acc_v.T

    return np.linalg.norm(v2, axis=1)

from typing import Any

import numpy as np
from numpy import floating
import time as pytime
import frame_transformation
from debug import pick


def compute_GPD(r_num: np.ndarray, r_sgp4: np.ndarray) -> np.ndarray:
    """
    Computes the Global Position Difference (GPD) between two orbits.

    Parameters:
        r_num  : Nx3 array of positions from numerical integrator (Euler, RK4, etc.)
        r_sgp4 : Nx3 array of positions from SGP4 propagation

    Returns:
        gpd : Nx1 array of L2 norms of the position residual at each epoch
    """

    if r_num.shape[0] != r_sgp4.shape[0]:
        raise RuntimeError("The two vectors have different rows!\n")
    else:
        if r_num.shape[1] != r_sgp4.shape[1]:
            raise RuntimeError("The two vectors have different columns!\n")

    residual = r_num - r_sgp4  # Nx3
    gpd = np.linalg.norm(residual, axis=1)  # Nx
    return gpd


def compute_CGPD(time: np.ndarray, GPD: np.ndarray) -> float | Any:
    """
    Computes the Global Position Difference (GPD) between two orbits.

    Parameters:
        time  : 1xN array of time
        GPD : 1xN array of errors

    Returns:
        cgpd : Nx1 array of RMS norms of the GPD
    """

    if time.size != GPD.size:
        raise RuntimeError("The two vectors have different rows!\n")

    mask = ~np.isnan(GPD)
    if not np.any(mask):
        return np.nan

    GPD_valid = GPD[mask]
    n = GPD_valid.size
    cgpd = np.linalg.norm(GPD_valid) / np.sqrt(n)
    return cgpd

import numpy as np
from numpy import floating
import time as pytime
import frame_transformation
from debug import pick

# ... compute_GPD and compute_CGPD stay as in your file ...


def run_GPD_iteration(
    N_iter,
    ti,
    tf_GPD,
    y_i,
    gravity_model,
    integrator_func,
    reference_orbit_func,
    KeplerianIterationResults,
    GPDIterationResults,
):
    """
    Single GPD iteration for a given number of steps N_iter.

    integrator_func      : function(y0, dt, tf, accel) -> (N+1,6) state array
    reference_orbit_func : function(time_grid, dt, tf) -> (N+1,6) state array
    """
    # Time grid
    time_iter = np.linspace(ti, tf_GPD, N_iter + 1)
    dt_iter = float(time_iter[1] - time_iter[0])

    # Reference orbit (SGP4, Euler, Heun, ...)
    y_ref_iter = reference_orbit_func(time_iter, dt_iter, tf_GPD)

    # Integrated orbit with chosen integrator (Euler, Heun, RK4, ...)
    start = pytime.time()
    y_num_iter = integrator_func(y_i, dt_iter, tf_GPD, gravity_model.point_mass)
    elapsed = pytime.time() - start

    if elapsed >= pick("max_int_time"):
        print(f"Computation time: {elapsed:.3f} s\n")
        raise RuntimeError(
            f"Integration with chosen integrator takes longer than {pick('max_int_time')} seconds; increase dt.\n"
        )

    # Position magnitude of integrated orbit
    y_pos_mag_iter = np.linalg.norm(y_num_iter[:, 0:3], axis=1)

    #  GPD calculation: numerical vs reference
    gpd = compute_GPD(y_num_iter[:, 0:3], y_ref_iter[:, 0:3])

    # Keplerian elements for both
    (a_num, e_num, _,
     i_num, _, Om_num, _,
     om_num, _, th_num, u_num) = frame_transformation.cartesian_to_keplerian(y_num_iter)

    (a_ref, e_ref, _,
     i_ref, _, Om_ref, _,
     om_ref, _, th_ref, u_ref) = frame_transformation.cartesian_to_keplerian(y_ref_iter)

    keplerian_results = KeplerianIterationResults(
        a_iter_euler=a_num,
        e_iter_euler=e_num,
        i_iter_euler=i_num,
        Om_iter_euler=Om_num,
        om_iter_euler=om_num,
        th_iter_euler=th_num,
        u_iter_euler=u_num,
        a_iter_SGP4=a_ref,
        e_iter_SGP4=e_ref,
        i_iter_SGP4=i_ref,
        Om_iter_SGP4=Om_ref,
        om_iter_SGP4=om_ref,
        th_iter_SGP4=th_ref,
        u_iter_SGP4=u_ref,
    )

    label_dt = f"dt = {dt_iter:.3f} [s], iterations = {N_iter + 1}"

    return GPDIterationResults(
        time_iter=time_iter,
        dt_iter=dt_iter,
        y_iter=y_num_iter,
        y_SGP4_iter=y_ref_iter,     # kept name for compatibility with Main
        gpd=gpd,
        y_pos_mag_iter=y_pos_mag_iter,
        label_dt=label_dt,
        keplerian_results=keplerian_results,
    )


def run_all_GPD_iterations(
    N_GPD,
    ti,
    tf_GPD,
    y_i,
    gravity_model,
    integrator_func,
    reference_orbit_func,
    KeplerianIterationResults,
    GPDIterationResults,
    GPDAllIterationsResults,
):
    """
    Runs GPD iterations for a list of step counts N_GPD, using a generic
    integrator and a generic reference orbit.
    """
    max_len = N_GPD[-1] + 1
    time_GPD = np.full((N_GPD.size, max_len), np.nan)
    GPD_num_ref = np.full((N_GPD.size, max_len), np.nan)
    y_num_iter_all = np.full((max_len, 6 * N_GPD.size), np.nan)
    y_pos_mag = np.full((N_GPD.size, max_len), np.nan)

    dt_GPD = []
    GPD_label_list = []

    a_num_list, e_num_list, i_num_list = [], [], []
    Om_num_list, om_num_list, th_num_list, u_num_list = [], [], [], []
    a_ref_list, e_ref_list, i_ref_list = [], [], []
    Om_ref_list, om_ref_list, th_ref_list, u_ref_list = [], [], [], []

    for k, N_iter in enumerate(N_GPD):
        iteration_results = run_GPD_iteration(
            N_iter=N_iter,
            ti=ti,
            tf_GPD=tf_GPD,
            y_i=y_i,
            gravity_model=gravity_model,
            integrator_func=integrator_func,
            reference_orbit_func=reference_orbit_func,
            KeplerianIterationResults=KeplerianIterationResults,
            GPDIterationResults=GPDIterationResults,
        )

        time_iter = iteration_results.time_iter
        dt_iter = iteration_results.dt_iter
        y_num_iter = iteration_results.y_iter
        gpd = iteration_results.gpd
        y_pos_mag_iter = iteration_results.y_pos_mag_iter
        label_dt = iteration_results.label_dt
        keplerian = iteration_results.keplerian_results

        n = time_iter.size
        dt_GPD.append(dt_iter)

        time_GPD[k, :n] = time_iter
        y_num_iter_all[:n, 6 * k:6 * (k + 1)] = y_num_iter
        y_pos_mag[k, :n] = y_pos_mag_iter
        GPD_num_ref[k, :n] = gpd

        a_num_list.append(keplerian.a_iter_euler)
        e_num_list.append(keplerian.e_iter_euler)
        i_num_list.append(keplerian.i_iter_euler)
        Om_num_list.append(keplerian.Om_iter_euler)
        om_num_list.append(keplerian.om_iter_euler)
        th_num_list.append(keplerian.th_iter_euler)
        u_num_list.append(keplerian.u_iter_euler)

        a_ref_list.append(keplerian.a_iter_SGP4)
        e_ref_list.append(keplerian.e_iter_SGP4)
        i_ref_list.append(keplerian.i_iter_SGP4)
        Om_ref_list.append(keplerian.Om_iter_SGP4)
        om_ref_list.append(keplerian.om_iter_SGP4)
        th_ref_list.append(keplerian.th_iter_SGP4)
        u_ref_list.append(keplerian.u_iter_SGP4)

        GPD_label_list.append(label_dt)

    return GPDAllIterationsResults(
        dt_GPD=dt_GPD,
        time_GPD=time_GPD,
        GPD_SGP4_Euler=GPD_num_ref,      # kept name for Main compatibility
        y_euler_iter=y_num_iter_all,     # "numerical" orbit
        y_pos_mag=y_pos_mag,
        GPD_label_list=GPD_label_list,
        a_euler_list=a_num_list,
        e_euler_list=e_num_list,
        i_euler_list=i_num_list,
        Om_euler_list=Om_num_list,
        om_euler_list=om_num_list,
        th_euler_list=th_num_list,
        u_euler_list=u_num_list,
        a_SGP4_list=a_ref_list,
        e_SGP4_list=e_ref_list,
        i_SGP4_list=i_ref_list,
        Om_SGP4_list=Om_ref_list,
        om_SGP4_list=om_ref_list,
        th_SGP4_list=th_ref_list,
        u_SGP4_list=u_ref_list,
    )


def run_CGPD_iteration(
    N_iter,
    ti,
    tf_GPD,
    y_i,
    gravity_model,
    integrator_func,
    reference_orbit_func,
    CGPDIterationResults,
):
    """
    Single CGPD iteration for a given N_iter using a generic integrator
    and a generic reference orbit.
    """
    time_iter = np.linspace(ti, tf_GPD, N_iter + 1)
    dt_iter = float(time_iter[1] - time_iter[0])

    # Reference orbit
    y_ref_iter = reference_orbit_func(time_iter, dt_iter, tf_GPD)

    # Numerical orbit
    start = pytime.time()
    y_num_iter = integrator_func(y_i, dt_iter, tf_GPD, gravity_model.point_mass)
    elapsed = pytime.time() - start

    if elapsed >= pick("max_int_time"):
        print(f"Computation time: {elapsed:.3f} s\n")
        raise RuntimeError(
            f"Integration with chosen integrator takes longer than {pick('max_int_time')} seconds; increase dt.\n"
        )

    y_pos_mag_iter = np.linalg.norm(y_num_iter[:, 0:3], axis=1)

    gpd = compute_GPD(y_num_iter[:, 0:3], y_ref_iter[:, 0:3])
    cgpd = compute_CGPD(time_iter, gpd)

    label_dt = f"dt = {dt_iter:.3f} [s], iterations = {N_iter + 1}"

    return CGPDIterationResults(
        time_iter=time_iter,
        dt_iter=dt_iter,
        y_iter=y_num_iter,
        y_SGP4_iter=y_ref_iter,   # kept name for compatibility
        cgpd=cgpd,
        y_pos_mag_iter=y_pos_mag_iter,
        label_dt=label_dt,
    )


def run_all_CGPD_iterations(
    N_CGPD,
    ti,
    tf_GPD,
    y_i,
    gravity_model,
    integrator_func,
    reference_orbit_func,
    CGPDIterationResults,
    CGPDAllIterationsResults,
):
    """
    Runs CGPD iterations for a list of N_CGPD, using a generic integrator
    and generic reference orbit.
    """
    max_len = N_CGPD[-1] + 1
    time_CGPD = np.full((N_CGPD.size, max_len), np.nan)
    CGPD_num_ref = np.zeros(N_CGPD.size)
    y_num_iter_all = np.full((max_len, 6 * N_CGPD.size), np.nan)
    y_pos_mag = np.full((N_CGPD.size, max_len), np.nan)

    dt_CGPD = []
    CGPD_label_list = []

    for k, N_iter in enumerate(N_CGPD):
        iteration_results = run_CGPD_iteration(
            N_iter=N_iter,
            ti=ti,
            tf_GPD=tf_GPD,
            y_i=y_i,
            gravity_model=gravity_model,
            integrator_func=integrator_func,
            reference_orbit_func=reference_orbit_func,
            CGPDIterationResults=CGPDIterationResults,
        )

        time_iter = iteration_results.time_iter
        dt_iter = iteration_results.dt_iter
        y_num_iter = iteration_results.y_iter
        cgpd = iteration_results.cgpd
        y_pos_mag_iter = iteration_results.y_pos_mag_iter
        label_dt = iteration_results.label_dt

        n = time_iter.size
        dt_CGPD.append(dt_iter)

        time_CGPD[k, :n] = time_iter
        y_num_iter_all[:n, 6 * k:6 * (k + 1)] = y_num_iter
        y_pos_mag[k, :n] = y_pos_mag_iter
        CGPD_num_ref[k] = cgpd

        CGPD_label_list.append(label_dt)

    return CGPDAllIterationsResults(
        dt_CGPD=dt_CGPD,
        time_CGPD=time_CGPD,
        CGPD_SGP4_Euler=CGPD_num_ref,   # kept name
        y_euler_iter=y_num_iter_all,    # "numerical" orbit
        y_pos_mag=y_pos_mag,
        CGPD_label_list=CGPD_label_list,
    )

# def run_GPD_iteration(N_iter, ti, tf_GPD, y_i, sat, gravity_model, integrator, sgp4_propagator,
#                       KeplerianIterationResults, GPDIterationResults):
#     # Time grid
#     time_iter = np.linspace(ti, tf_GPD, N_iter + 1)
#     dt_iter = float(time_iter[1] - time_iter[0])
#
#     #  SGP4 orbit
#     r_SGP4_E, v_SGP4_E = sgp4_propagator.propagated_position_velocity(sat, time_iter)
#     y_SGP4_iter = np.concatenate((r_SGP4_E, v_SGP4_E), axis=1)
#
#     #  Euler orbit
#     start = pytime.time()
#     y_iter = integrator.Euler(y_i, dt_iter, tf_GPD, gravity_model.point_mass)
#     elapsed = pytime.time() - start
#
#     if elapsed >= pick("max_int_time"):
#         print(f"Computation time: {elapsed:.3f} s\n")
#         raise RuntimeError(f"Integration with Euler takes longer than {pick("max_int_time")} seconds; increase dt.\n")
#
#     # position magnitude
#     y_pos_mag_iter = np.linalg.norm(y_iter[:, 0:3], axis=1)
#
#     #  GPD calculation
#     gpd = compute_GPD(y_iter[:, 0:3], r_SGP4_E)
#
#     # Keplerian elements calculation:
#     (a_iter_euler, e_iter_euler, _,
#      i_iter_euler, _, Om_iter_euler, _,
#      om_iter_euler, _, th_iter_euler, u_iter_euler) = frame_transformation.cartesian_to_keplerian(y_iter)
#
#     (a_iter_SGP4, e_iter_SGP4, _,
#      i_iter_SGP4, _, Om_iter_SGP4, _,
#      om_iter_SGP4, _, th_iter_SGP4, u_iter_SGP4) = frame_transformation.cartesian_to_keplerian(y_SGP4_iter)
#
#     keplerian_results = KeplerianIterationResults(
#         a_iter_euler=a_iter_euler,
#         e_iter_euler=e_iter_euler,
#         i_iter_euler=i_iter_euler,
#         Om_iter_euler=Om_iter_euler,
#         om_iter_euler=om_iter_euler,
#         th_iter_euler=th_iter_euler,
#         u_iter_euler=u_iter_euler,
#         a_iter_SGP4=a_iter_SGP4,
#         e_iter_SGP4=e_iter_SGP4,
#         i_iter_SGP4=i_iter_SGP4,
#         Om_iter_SGP4=Om_iter_SGP4,
#         om_iter_SGP4=om_iter_SGP4,
#         th_iter_SGP4=th_iter_SGP4,
#         u_iter_SGP4=u_iter_SGP4,
#     )
#
#     label_dt = f"dt = {dt_iter:.3f} [s], iterations = {N_iter + 1}"
#
#     return GPDIterationResults(
#         time_iter=time_iter,
#         dt_iter=dt_iter,
#         y_iter=y_iter,
#         y_SGP4_iter=y_SGP4_iter,
#         gpd=gpd,
#         y_pos_mag_iter=y_pos_mag_iter,
#         label_dt=label_dt,
#         keplerian_results=keplerian_results,
#     )
#
#
# def run_all_GPD_iterations(N_GPD, ti, tf_GPD, y_i, sat, gravity_model, integrator, sgp4_propagator,
#                            KeplerianIterationResults, GPDIterationResults, GPDAllIterationsResults):
#     # allocate matrices with max_len = max N_iter + 1
#     max_len = N_GPD[-1] + 1
#     time_GPD = np.full((N_GPD.size, max_len), np.nan)
#     GPD_SGP4_Euler = np.full((N_GPD.size, max_len), np.nan)
#     y_euler_iter = np.full((max_len, 6 * N_GPD.size), np.nan)
#     y_pos_mag = np.full((N_GPD.size, max_len), np.nan)
#
#     dt_GPD = []
#     GPD_label_list = []
#
#     a_euler_list, e_euler_list, i_euler_list, Om_euler_list, om_euler_list, th_euler_list, u_euler_list = [], [], [], [], [], [], []
#     a_SGP4_list, e_SGP4_list, i_SGP4_list, Om_SGP4_list, om_SGP4_list, th_SGP4_list, u_SGP4_list = [], [], [], [], [], [], []
#
#     # internal loop, hidden from main
#     for k, N_iter in enumerate(N_GPD):
#         iteration_results = run_GPD_iteration(
#             N_iter=N_iter,
#             ti=ti,
#             tf_GPD=tf_GPD,
#             y_i=y_i,
#             sat=sat,
#             gravity_model=gravity_model,
#             integrator=integrator,
#             sgp4_propagator=sgp4_propagator,
#             KeplerianIterationResults=KeplerianIterationResults,
#             GPDIterationResults=GPDIterationResults,
#         )
#
#         time_iter = iteration_results.time_iter
#         dt_iter = iteration_results.dt_iter
#         y_iter = iteration_results.y_iter
#         gpd = iteration_results.gpd
#         y_pos_mag_iter = iteration_results.y_pos_mag_iter
#         label_dt = iteration_results.label_dt
#         keplerian = iteration_results.keplerian_results
#
#         n = time_iter.size
#         dt_GPD.append(dt_iter)
#
#         # save time grid
#         time_GPD[k, :n] = time_iter
#
#         # store full Euler state
#         y_euler_iter[:n, 6 * k:6 * (k + 1)] = y_iter
#
#         # position magnitude (always length n)
#         y_pos_mag[k, :n] = y_pos_mag_iter
#
#         # GPD
#         GPD_SGP4_Euler[k, :n] = gpd
#
#         # Keplerian lists
#         a_euler_list.append(keplerian.a_iter_euler)
#         e_euler_list.append(keplerian.e_iter_euler)
#         i_euler_list.append(keplerian.i_iter_euler)
#         Om_euler_list.append(keplerian.Om_iter_euler)
#         om_euler_list.append(keplerian.om_iter_euler)
#         th_euler_list.append(keplerian.th_iter_euler)
#         u_euler_list.append(keplerian.u_iter_euler)
#
#         a_SGP4_list.append(keplerian.a_iter_SGP4)
#         e_SGP4_list.append(keplerian.e_iter_SGP4)
#         i_SGP4_list.append(keplerian.i_iter_SGP4)
#         Om_SGP4_list.append(keplerian.Om_iter_SGP4)
#         om_SGP4_list.append(keplerian.om_iter_SGP4)
#         th_SGP4_list.append(keplerian.th_iter_SGP4)
#         u_SGP4_list.append(keplerian.u_iter_SGP4)
#
#         GPD_label_list.append(label_dt)
#
#     return GPDAllIterationsResults(
#         dt_GPD=dt_GPD,
#         time_GPD=time_GPD,
#         GPD_SGP4_Euler=GPD_SGP4_Euler,
#         y_euler_iter=y_euler_iter,
#         y_pos_mag=y_pos_mag,
#         GPD_label_list=GPD_label_list,
#         a_euler_list=a_euler_list,
#         e_euler_list=e_euler_list,
#         i_euler_list=i_euler_list,
#         Om_euler_list=Om_euler_list,
#         om_euler_list=om_euler_list,
#         th_euler_list=th_euler_list,
#         u_euler_list=u_euler_list,
#         a_SGP4_list=a_SGP4_list,
#         e_SGP4_list=e_SGP4_list,
#         i_SGP4_list=i_SGP4_list,
#         Om_SGP4_list=Om_SGP4_list,
#         om_SGP4_list=om_SGP4_list,
#         th_SGP4_list=th_SGP4_list,
#         u_SGP4_list=u_SGP4_list,
#     )
#
# def run_GPD_iteration(
#     N_iter,
#     ti,
#     tf_GPD,
#     y_i,
#     gravity_model,
#     integrator_func,
#     reference_orbit_func,
#     KeplerianIterationResults,
#     GPDIterationResults,
# ):
#     # Time grid
#     time_iter = np.linspace(ti, tf_GPD, N_iter + 1)
#     dt_iter = float(time_iter[1] - time_iter[0])
#
#     # Reference orbit (e.g. SGP4, Euler, etc.):
#     # must return a (N+1, 6) array [rx ry rz vx vy vz]
#     y_SGP4_iter = reference_orbit_func(time_iter, dt_iter, tf_GPD)
#
#     #  Integrated orbit with chosen integrator (Euler, Heun, ...)
#     start = pytime.time()
#     y_iter = integrator_func(y_i, dt_iter, tf_GPD, gravity_model.point_mass)
#     elapsed = pytime.time() - start
#
#     if elapsed >= pick("max_int_time"):
#         print(f"Computation time: {elapsed:.3f} s\n")
#         raise RuntimeError(
#             f"Integration with chosen integrator takes longer than {pick('max_int_time')} seconds; increase dt.\n"
#         )
# #
#     # Position magnitude of integrated orbit
#     y_pos_mag_iter = np.linalg.norm(y_iter[:, 0:3], axis=1)
#
#     #  GPD calculation: integrated vs reference
#     gpd = compute_GPD(y_iter[:, 0:3], y_SGP4_iter[:, 0:3])
#
#     # Keplerian elements calculation:
#     (a_iter_euler, e_iter_euler, _,
#      i_iter_euler, _, Om_iter_euler, _,
#      om_iter_euler, _, th_iter_euler, u_iter_euler) = frame_transformation.cartesian_to_keplerian(y_iter)
#
#     (a_iter_SGP4, e_iter_SGP4, _,
#      i_iter_SGP4, _, Om_iter_SGP4, _,
#      om_iter_SGP4, _, th_iter_SGP4, u_iter_SGP4) = frame_transformation.cartesian_to_keplerian(y_SGP4_iter)
#
#     keplerian_results = KeplerianIterationResults(
#         a_iter_euler=a_iter_euler,
#         e_iter_euler=e_iter_euler,
#         i_iter_euler=i_iter_euler,
#         Om_iter_euler=Om_iter_euler,
#         om_iter_euler=om_iter_euler,
#         th_iter_euler=th_iter_euler,
#         u_iter_euler=u_iter_euler,
#         a_iter_SGP4=a_iter_SGP4,
#         e_iter_SGP4=e_iter_SGP4,
#         i_iter_SGP4=i_iter_SGP4,
#         Om_iter_SGP4=Om_iter_SGP4,
#         om_iter_SGP4=om_iter_SGP4,
#         th_iter_SGP4=th_iter_SGP4,
#         u_iter_SGP4=u_iter_SGP4,
#     )
#
#     label_dt = f"dt = {dt_iter:.3f} [s], iterations = {N_iter + 1}"
#
#     return GPDIterationResults(
#         time_iter=time_iter,
#         dt_iter=dt_iter,
#         y_iter=y_iter,
#         y_SGP4_iter=y_SGP4_iter,
#         gpd=gpd,
#         y_pos_mag_iter=y_pos_mag_iter,
#         label_dt=label_dt,
#         keplerian_results=keplerian_results,
#     )
#
# def run_all_GPD_iterations(
#     N_GPD,
#     ti,
#     tf_GPD,
#     y_i,
#     gravity_model,
#     integrator_func,
#     reference_orbit_func,
#     KeplerianIterationResults,
#     GPDIterationResults,
#     GPDAllIterationsResults,
# ):
#     # allocate matrices with max_len = max N_iter + 1
#     max_len = N_GPD[-1] + 1
#     time_GPD = np.full((N_GPD.size, max_len), np.nan)
#     GPD_SGP4_Euler = np.full((N_GPD.size, max_len), np.nan)
#     y_euler_iter = np.full((max_len, 6 * N_GPD.size), np.nan)
#     y_pos_mag = np.full((N_GPD.size, max_len), np.nan)
#
#     dt_GPD = []
#     GPD_label_list = []
#
#     a_euler_list, e_euler_list, i_euler_list, Om_euler_list, om_euler_list, th_euler_list, u_euler_list = [], [], [], [], [], [], []
#     a_SGP4_list, e_SGP4_list, i_SGP4_list, Om_SGP4_list, om_SGP4_list, th_SGP4_list, u_SGP4_list = [], [], [], [], [], [], []
#
#     # internal loop, hidden from main
#     for k, N_iter in enumerate(N_GPD):
#         iteration_results = run_GPD_iteration(
#             N_iter=N_iter,
#             ti=ti,
#             tf_GPD=tf_GPD,
#             y_i=y_i,
#             gravity_model=gravity_model,
#             integrator_func=integrator_func,
#             reference_orbit_func=reference_orbit_func,
#             KeplerianIterationResults=KeplerianIterationResults,
#             GPDIterationResults=GPDIterationResults,
#         )
#
#         time_iter = iteration_results.time_iter
#         dt_iter = iteration_results.dt_iter
#         y_iter = iteration_results.y_iter
#         gpd = iteration_results.gpd
#         y_pos_mag_iter = iteration_results.y_pos_mag_iter
#         label_dt = iteration_results.label_dt
#         keplerian = iteration_results.keplerian_results
#
#         n = time_iter.size
#         dt_GPD.append(dt_iter)
#
#         # save time grid
#         time_GPD[k, :n] = time_iter
#
#         # store full integrated state (still called "euler" for compatibility)
#         y_euler_iter[:n, 6 * k:6 * (k + 1)] = y_iter
#
#         # position magnitude (always length n)
#         y_pos_mag[k, :n] = y_pos_mag_iter
#
#         # GPD
#         GPD_SGP4_Euler[k, :n] = gpd
#
#         # Keplerian lists
#         a_euler_list.append(keplerian.a_iter_euler)
#         e_euler_list.append(keplerian.e_iter_euler)
#         i_euler_list.append(keplerian.i_iter_euler)
#         Om_euler_list.append(keplerian.Om_iter_euler)
#         om_euler_list.append(keplerian.om_iter_euler)
#         th_euler_list.append(keplerian.th_iter_euler)
#         u_euler_list.append(keplerian.u_iter_euler)
#
#         a_SGP4_list.append(keplerian.a_iter_SGP4)
#         e_SGP4_list.append(keplerian.e_iter_SGP4)
#         i_SGP4_list.append(keplerian.i_iter_SGP4)
#         Om_SGP4_list.append(keplerian.Om_iter_SGP4)
#         om_SGP4_list.append(keplerian.om_iter_SGP4)
#         th_SGP4_list.append(keplerian.th_iter_SGP4)
#         u_SGP4_list.append(keplerian.u_iter_SGP4)
#
#         GPD_label_list.append(label_dt)
#
#     return GPDAllIterationsResults(
#         dt_GPD=dt_GPD,
#         time_GPD=time_GPD,
#         GPD_SGP4_Euler=GPD_SGP4_Euler,
#         y_euler_iter=y_euler_iter,
#         y_pos_mag=y_pos_mag,
#         GPD_label_list=GPD_label_list,
#         a_euler_list=a_euler_list,
#         e_euler_list=e_euler_list,
#         i_euler_list=i_euler_list,
#         Om_euler_list=Om_euler_list,
#         om_euler_list=om_euler_list,
#         th_euler_list=th_euler_list,
#         u_euler_list=u_euler_list,
#         a_SGP4_list=a_SGP4_list,
#         e_SGP4_list=e_SGP4_list,
#         i_SGP4_list=i_SGP4_list,
#         Om_SGP4_list=Om_SGP4_list,
#         om_SGP4_list=om_SGP4_list,
#         th_SGP4_list=th_SGP4_list,
#         u_SGP4_list=u_SGP4_list,
#     )
#
#
#
#
# def run_all_CGPD_iterations(N_CGPD, ti, tf_GPD, y_i, sat, gravity_model, integrator, sgp4_propagator,
#                             CGPDIterationResults, CGPDAllIterationsResults):
#     # allocate matrices with max_len = max N_iter + 1
#     max_len = N_CGPD[-1] + 1
#     time_CGPD = np.full((N_CGPD.size, max_len), np.nan)
#     CGPD_SGP4_Euler = np.zeros(N_CGPD.size)
#     y_euler_iter = np.full((max_len, 6 * N_CGPD.size), np.nan)
#     y_pos_mag = np.full((N_CGPD.size, max_len), np.nan)
#
#     dt_CGPD = []
#     CGPD_label_list = []
#
#     # internal loop
#     for k, N_iter in enumerate(N_CGPD):
#         iteration_results = run_CGPD_iteration(
#             N_iter=N_iter,
#             ti=ti,
#             tf_GPD=tf_GPD,
#             y_i=y_i,
#             sat=sat,
#             gravity_model=gravity_model,
#             integrator=integrator,
#             sgp4_propagator=sgp4_propagator,
#             CGPDIterationResults=CGPDIterationResults,
#         )
#
#         time_iter = iteration_results.time_iter
#         dt_iter = iteration_results.dt_iter
#         y_iter = iteration_results.y_iter
#         cgpd = iteration_results.cgpd
#         y_pos_mag_iter = iteration_results.y_pos_mag_iter
#         label_dt = iteration_results.label_dt
#
#         n = time_iter.size
#         dt_CGPD.append(dt_iter)
#
#         # save time grid
#         time_CGPD[k, :n] = time_iter
#
#         # store full Euler state
#         y_euler_iter[:n, 6 * k:6 * (k + 1)] = y_iter
#
#         # position magnitude (always length n)
#         y_pos_mag[k, :n] = y_pos_mag_iter
#
#         # GPD
#         CGPD_SGP4_Euler[k] = cgpd
#
#         CGPD_label_list.append(label_dt)
#
#     return CGPDAllIterationsResults(
#         dt_CGPD=dt_CGPD,
#         time_CGPD=time_CGPD,
#         CGPD_SGP4_Euler=CGPD_SGP4_Euler,
#         y_euler_iter=y_euler_iter,
#         y_pos_mag=y_pos_mag,
#         CGPD_label_list=CGPD_label_list,
#     )
#
#
# def run_CGPD_iteration(N_iter, ti, tf_GPD, y_i, sat, gravity_model, integrator, sgp4_propagator, CGPDIterationResults):
#     # Time grid
#     time_iter = np.linspace(ti, tf_GPD, N_iter + 1)
#     dt_iter = float(time_iter[1] - time_iter[0])
#
#     #  SGP4 orbit
#     r_SGP4_E, v_SGP4_E = sgp4_propagator.propagated_position_velocity(sat, time_iter)
#     y_SGP4_iter = np.concatenate((r_SGP4_E, v_SGP4_E), axis=1)
#
#     #  Euler orbit
#     start = pytime.time()
#     y_iter = integrator.Euler(y_i, dt_iter, tf_GPD, gravity_model.point_mass)
#     elapsed = pytime.time() - start
#
#     if elapsed >= pick("max_int_time"):
#         print(f"Computation time: {elapsed:.3f} s\n")
#         raise RuntimeError(f"Integration with Euler takes longer than {pick("max_int_time")} seconds; increase dt.\n")
#
#     # position magnitude
#     y_pos_mag_iter = np.linalg.norm(y_iter[:, 0:3], axis=1)
#
#     #  GPD calculation
#     gpd = compute_GPD(y_iter[:, 0:3], r_SGP4_E)
#     cgpd = compute_CGPD(time_iter, gpd)
#
#     label_dt = f"dt = {dt_iter:.3f} [s], iterations = {N_iter + 1}"
#
#     return CGPDIterationResults(
#         time_iter=time_iter,
#         dt_iter=dt_iter,
#         y_iter=y_iter,
#         y_SGP4_iter=y_SGP4_iter,
#         cgpd=cgpd,
#         y_pos_mag_iter=y_pos_mag_iter,
#         label_dt=label_dt,
#     )
