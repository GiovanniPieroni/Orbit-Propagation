import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as mticker
from matplotlib.ticker import LogLocator, NullFormatter
from modules import errors
import utils
from matplotlib.lines import Line2D  # <--- Make sure this is imported at the top

from debug import pick

formatter = mticker.ScalarFormatter(useMathText=True)
formatter.set_powerlimits((-3, 2))

# Constants for consistent plotting styles
LABEL_SIZE = 16
TITLE_SIZE = 18
TICK_SIZE = 14
LEGEND_SIZE = 'large'

with open("TLE_GONETSM24.txt", "r") as f:
    TLE_data = f.read().splitlines()

sat = utils.create_satrec(TLE_data)

mean_motion = sat.no_kozai / 60  # [rad/s]
period = 2 * np.pi / mean_motion  # Orbital period


def mask_discontinuity(data, threshold=179):
    """
    Replaces values with NaN where the jump between points exceeds threshold.
    Useful for removing vertical lines in -180/180 wrapping angles.
    """
    # Create a copy so we don't modify the original data array
    clean_data = data.copy()

    # Calculate difference between consecutive points
    # prepend=data[0] keeps the shape consistent
    diff = np.abs(np.diff(clean_data, prepend=clean_data[0]))

    # Find indices where the jump is too big (wrapping)
    # Threshold 300 deg is safe for standard orbits
    jump_indices = diff > threshold

    # Set those points to NaN to break the line
    clean_data[jump_indices] = np.nan
    return clean_data


def plot_two_scales(time, data1, label1, ylabel1, data2, label2, ylabel2,
                    title="Orbital Elements", filename="plot_2_axes.png"):
    """
    Plots 2 quantities with different scales on a single plot using 2 Y-axes.
    """

    if "Latitude" in title or "Argument" in title:
        data2 = mask_discontinuity(data2)

    fig, ax1 = plt.subplots(figsize=(10, 8))  # Taller figure
    # Adjust top margin significantly to separate Title, Legend, and 1T labels
    fig.subplots_adjust(top=0.80)

    color1, color2 = 'tab:blue', 'tab:red'

    # --- Axis 1 ---
    ax1.plot(time, data1, color=color1, label=label1, linestyle='-')
    ax1.set_xlabel("Time [s]", fontsize=LABEL_SIZE)
    ax1.set_ylabel(ylabel1, color=color1, fontsize=LABEL_SIZE)
    ax1.tick_params(axis='y', labelcolor=color1, labelsize=TICK_SIZE)
    ax1.tick_params(axis='x', labelsize=TICK_SIZE)
    ax1.grid(True, linestyle=':', alpha=0.6)

    # --- Axis 2 ---
    ax2 = ax1.twinx()
    ax2.plot(time, data2, color=color2, label=label2, linestyle='-')
    ax2.set_ylabel(ylabel2, color=color2, fontsize=LABEL_SIZE)
    ax2.tick_params(axis='y', labelcolor=color2, labelsize=TICK_SIZE)

    # --- Legend & Title ---
    lines = [ax1.get_lines()[0], ax2.get_lines()[0]]
    labels = [l.get_label() for l in lines]

    # --- ADD ORBITAL PERIOD LABEL ---
    # lines.append(Line2D([0], [0], color='none', label='1T = one orbital period'))
    # labels.append('1T = one orbital period')

    # Place legend ABOVE the plot area but BELOW the title
    # y=1.08 puts it clear of the x-axis top frame where "1T" labels sit
    ax1.legend(lines, labels, loc='lower center', bbox_to_anchor=(0.5, 1.08),
               ncol=2, borderaxespad=0., frameon=True, fontsize=LEGEND_SIZE)

    # Push title way up
    ax1.set_title(title, y=1.20, fontsize=TITLE_SIZE)

    fig.savefig(filename, bbox_inches="tight")
    plt.close(fig)
    return str(filename)


def plot_three_scales(time,
                      data1, label1, ylabel1,
                      data2, label2, ylabel2,
                      data3, label3, ylabel3,
                      title="Orbital Elements",
                      filename="plot_3_axes.png",
                      ax1_big=False):
    # --- Setup & Figure Size ---
    max_time = np.nanmax(np.array(time)) if len(time) > 0 else 0
    orbital_lines = []

    # Calculate lines
    if period is not None and max_time > period:
        orbital_lines = np.arange(period, max_time + 1, period, dtype=float)
        fig, ax1 = plt.subplots(figsize=(10 + len(orbital_lines), 9))
    else:
        fig, ax1 = plt.subplots(figsize=(10, 8))  # Increased height

    fig.subplots_adjust(right=0.75, top=0.82)  # Adjust top margin

    try:
        # --- Masking Logic ---
        if "Latitude" in title or "Argument" in title:
            data1 = mask_discontinuity(data1)
            data2 = mask_discontinuity(data2)

        color1, color2, color3 = 'tab:blue', 'tab:green', 'tab:orange'

        # --- Axis 2 (Right) ---
        if "Latitude" in title or "Argument" in title:
            ax2 = ax1.twinx()
            ax2.plot(time, data2, color=color2, label=label2, linestyle='--', linewidth=4, zorder=0)
            ax2.set_ylabel(ylabel2, color=color2, fontsize=LABEL_SIZE)
            ax2.yaxis.set_minor_locator(mticker.AutoMinorLocator())
            ax2.tick_params(axis='y', labelcolor=color2, labelsize=TICK_SIZE)
            ax2.yaxis.get_offset_text().set_x(1)
        else:
            ax2 = ax1.twinx()
            if ax1_big:
                ax2.plot(time, data2, color=color2, label=label2, linestyle='-', linewidth=6, zorder=0, alpha=0.5)
            else:
                ax2.plot(time, data2, color=color2, label=label2, linestyle='-', linewidth=6, zorder=0)

            ax2.set_ylabel(ylabel2, color=color2, fontsize=LABEL_SIZE)
            ax2.yaxis.set_minor_locator(mticker.AutoMinorLocator())
            ax2.tick_params(axis='y', labelcolor=color2, labelsize=TICK_SIZE)
            ax2.yaxis.get_offset_text().set_x(1)

        # --- Axis 1 (Left) ---
        if "Latitude" in title or "Argument" in title:
            ax1.plot(time, data1, color=color1, label=label1, linestyle='-', linewidth=6, zorder=1)
        else:
            if ax1_big:
                ax1.plot(time, data1, color="tab:blue", label=label1, linestyle='--', linewidth=6, zorder=1)
            else:
                ax1.plot(time, data1, color=color1, label=label1, linestyle='-', linewidth=2, zorder=1)

        ax1.set_xlabel("Time [s]", fontsize=LABEL_SIZE)

        # FIX 1: Massive padding to make room for the offset text on the left
        ax1.set_ylabel(ylabel1, color=color1, fontsize=LABEL_SIZE)

        # Ticks
        ax1.xaxis.set_minor_locator(mticker.AutoMinorLocator())
        ax1.yaxis.set_minor_locator(mticker.AutoMinorLocator())
        ax1.tick_params(axis='both', which='major', length=7, width=1.2, labelsize=TICK_SIZE)
        ax1.tick_params(axis='both', which='minor', length=4, width=0.8)
        ax1.tick_params(axis='y', labelcolor=color1)
        ax1.grid(True, which='major', linestyle=':', alpha=0.6)

        # FIX 2: Move Scientific Notation (Offset Text) FAR LEFT
        # (-0.25 means 25% of the plot width to the left of the y-axis)
        ax1.yaxis.get_offset_text().set_horizontalalignment('right')
        ax1.yaxis.get_offset_text().set_position((-0.25, 1.0))

        # --- Vertical Lines (Orbital Period) ---
        if period is not None and max_time > period:
            for line_x in orbital_lines:
                ax1.axvline(x=line_x, color='black', linestyle='-.', alpha=0.4, zorder=0)

                # FIX 3: Place text using Axis Coordinates (y=1.01) so it sits just above the frame
                # This separates it vertically from the scientific notation if they are close
                label_text = rf"{int(line_x / period)} T $\approx$ {int(period * (line_x / period))} [s]"
                ax1.text(line_x, 1.04, label_text,
                         transform=ax1.get_xaxis_transform(),  # <--- Key fix: uses axis y-coords (0 to 1)
                         ha='center', va='bottom', fontsize=12, color='gray')

        # --- Axis 3 (Right Offset) ---
        ax3 = ax1.twinx()
        ax3.spines.right.set_position(("axes", 1.15))
        ax3.set_frame_on(True)
        ax3.patch.set_visible(False)
        ax3.plot(time, data3, color=color3, label=label3, linestyle='--', linewidth=4)
        ax3.set_ylabel(ylabel3, color=color3, fontsize=LABEL_SIZE)
        ax3.yaxis.set_minor_locator(mticker.AutoMinorLocator())
        ax3.tick_params(axis='y', labelcolor=color3, labelsize=TICK_SIZE)

        # Formatter
        ax3.yaxis.set_major_formatter(formatter)
        ax3.yaxis.get_offset_text().set_x(1.18)

        # --- Legend & Title ---
        lines = [ax1.get_lines()[0], ax2.get_lines()[0], ax3.get_lines()[0]]
        labels = [l.get_label() for l in lines]

        # --- ADD ORBITAL PERIOD LABEL ---
        lines.append(Line2D([0], [0], color='none', label='1T = one orbital period'))
        labels.append('1T = one orbital period')

        # Place legend ABOVE the plot, UNDER the title
        ax1.legend(lines, labels, loc='lower center', bbox_to_anchor=(0.5, 1.08),
                   ncol=3, borderaxespad=0., frameon=True, fontsize=LEGEND_SIZE)

        # Push title higher
        ax1.set_title(title, y=1.22, fontsize=TITLE_SIZE)

        # --- Save ---
        try:
            fig.savefig(filename, bbox_inches="tight")
        except ValueError:
            print(f"Warning: Infinite data detected in {filename}. Saving without 'tight' layout.")
            fig.savefig(filename)

    finally:
        plt.close(fig)

    return str(filename)


def plot_three_scales_mat(time,
                          data1, label1, ylabel1,
                          data2, label2, ylabel2,
                          data3, label3, ylabel3,
                          title="Orbital Elements", filename="plot_3_axes.png"):
    # --- 1. Data Prep ---
    time = np.asarray(time)
    data1 = np.asarray(data1)
    data2 = np.asarray(data2)
    data3 = np.asarray(data3)

    if time.ndim != 2:
        raise ValueError("time must be a 2D array (n_runs, N)")

    n_runs = time.shape[0]

    # --- 2. Mask Discontinuities ---
    if "Latitude" in title or "Argument" in title:
        for k in range(n_runs):
            data1[k] = mask_discontinuity(data1[k])
            data2[k] = mask_discontinuity(data2[k])
            data3[k] = mask_discontinuity(data3[k])

    # --- 3. Figure Setup ---
    max_time = np.nanmax(time)
    orbital_lines = []

    # Check for global 'period' variable, otherwise skip lines
    if 'period' in globals() and period is not None and max_time > period:
        nT = int(max_time // period)
        orbital_lines = np.arange(1, nT + 1) * period
        fig, ax1 = plt.subplots(figsize=(10 + nT, 9))
    else:
        fig, ax1 = plt.subplots(figsize=(10, 8))

    fig.subplots_adjust(right=0.75, top=0.82, bottom=0.15)  # Adjusted margins

    # --- 4. Style: Colors & Widths ---
    # Gradient: Light -> Dark
    c_indices = np.linspace(0.4, 1.0, n_runs)
    colors1 = [plt.cm.Blues(i) for i in c_indices]
    colors2 = [plt.cm.Greens(i) for i in c_indices]
    colors3 = [plt.cm.Oranges(i) for i in c_indices]

    # --- 5. DT Logic (Reverted to simple calculation) ---
    # We calculate dt for all runs upfront to determine plotting order
    dt_values = time[:, 1] - time[:, 0]

    # Sort: Largest dt (Thick/Light) -> Smallest dt (Thin/Dark)
    # This ensures the high-precision (small dt) run is on top.
    order = np.argsort(dt_values)[::-1]

    # Linewidths: Thick -> Thin
    line_widths = np.linspace(6.0, 1.5, n_runs)

    # Lists for Legend (We will append in loop order)
    legend_handles = []
    legend_labels = []

    # --- Axis 1 (Left, Blue) ---
    ax1.set_xlabel("Time [s]", fontsize=LABEL_SIZE)
    ax1.set_ylabel(ylabel1, color='tab:blue', fontsize=LABEL_SIZE)
    ax1.tick_params(axis='y', labelcolor='tab:blue', labelsize=TICK_SIZE)
    ax1.tick_params(axis='x', labelsize=TICK_SIZE)
    ax1.grid(True, linestyle=':', alpha=0.6)

    ax1.yaxis.get_offset_text().set_horizontalalignment('right')
    ax1.yaxis.get_offset_text().set_position((-0.25, 1.0))
    ax1.xaxis.set_minor_locator(mticker.AutoMinorLocator())
    ax1.yaxis.set_minor_locator(mticker.AutoMinorLocator())

    for i, k in enumerate(order):
        dt = dt_values[k]  # Use the simple dt value
        lbl = f"{label1}, dt = {dt:.3f} s"

        h, = ax1.plot(time[k], data1[k],
                      color=colors1[i],
                      linewidth=line_widths[i],
                      label=lbl,
                      zorder=10 + i)
        legend_handles.append(h)
        legend_labels.append(h.get_label())

    # --- Axis 2 (Right, Green) ---
    ax2 = ax1.twinx()
    ax2.set_ylabel(ylabel2, color='tab:green', fontsize=LABEL_SIZE)
    ax2.tick_params(axis='y', labelcolor='tab:green', labelsize=TICK_SIZE)
    ax2.yaxis.get_offset_text().set_x(1.1)
    ax2.yaxis.set_minor_locator(mticker.AutoMinorLocator())

    for i, k in enumerate(order):
        dt = dt_values[k]
        lbl = f"{label2}, dt = {dt:.3f} s"

        h, = ax2.plot(time[k], data2[k],
                      color=colors2[i],
                      linewidth=line_widths[i],
                      linestyle='-',
                      label=lbl,
                      zorder=10 + i)
        legend_handles.append(h)
        legend_labels.append(h.get_label())

    # --- Axis 3 (Right Offset, Orange) ---
    ax3 = ax1.twinx()
    ax3.spines.right.set_position(("axes", 1.15))
    ax3.set_ylabel(ylabel3, color='tab:orange', fontsize=LABEL_SIZE)
    ax3.tick_params(axis='y', labelcolor='tab:orange', labelsize=TICK_SIZE)
    ax3.set_frame_on(True)
    ax3.patch.set_visible(False)
    ax3.yaxis.get_offset_text().set_x(1.18)

    ax3.yaxis.set_major_formatter(formatter)
    ax3.yaxis.get_offset_text().set_x(1.18)
    ax3.yaxis.set_minor_locator(mticker.AutoMinorLocator())  # Added minor ticks for ax3

    for i, k in enumerate(order):
        dt = dt_values[k]
        lbl = f"{label3}, dt = {dt:.3f} s"

        h, = ax3.plot(time[k], data3[k],
                      color=colors3[i],
                      linewidth=line_widths[i],
                      linestyle='--',
                      label=lbl,
                      zorder=10 + i)
        legend_handles.append(h)
        legend_labels.append(h.get_label())

    # --- Orbital Period Lines ---
    if 'period' in globals() and period is not None and len(orbital_lines) > 0:
        for lx in orbital_lines:
            ax1.axvline(lx, color='black', linestyle='-.', alpha=0.3, zorder=0)
            ax1.text(lx, 1.05, f"{int(lx / period)} T ≈ {int(lx)} s",
                     transform=ax1.get_xaxis_transform(),
                     ha='center', va='bottom', fontsize=12, color='gray')

    # --- Legend ---
    # --- ADD ORBITAL PERIOD LABEL ---
    legend_handles.append(Line2D([0], [0], color='none', label='1T = one orbital period'))
    legend_labels.append('1T = one orbital period')

    # Place legend ABOVE the plot, UNDER the title
    ax1.legend(legend_handles, legend_labels,
               loc='lower center', bbox_to_anchor=(0.5, 1.08),  # Adjusted higher to clear labels
               fontsize=LEGEND_SIZE, ncol=3, frameon=True)

    ax1.set_title(title, y=1.22, fontsize=TITLE_SIZE)  # Push title higher

    try:
        fig.savefig(filename, bbox_inches="tight")
    except ValueError:
        print(f"Warning: Infinite data detected in {filename}. Saving without 'tight' layout.")
        fig.savefig(filename)

    plt.close(fig)
    return str(filename)



def plot_normal(x: np.ndarray, y: np.ndarray,
                xlabel: str, ylabel: str, label_list: list,
                title="Generic title", filename="Generic_filename.png"):
    """
    Args:
        x: (N,) or (#args, N) array
        y: (N,) or (#args, N) array
        label_list: List of label strings
        filename: Figure filename
    """
    # --- 3. Figure Setup ---
    max_time = np.nanmax(x)
    orbital_lines = []
    # Calculate lines
    if period is not None and max_time > period:
        orbital_lines = np.arange(period, max_time + 1, period, dtype=float)
        fig, ax = plt.subplots(figsize=(7 + len(orbital_lines), 8))  # Increased height
    else:
        fig, ax = plt.subplots(figsize=(7, 8))  # Increased height

    # Increase top margin to give room for labels, legend, and title stack
    fig.subplots_adjust(top=0.80)

    # --- Force numpy arrays ---
    x = np.asarray(x)
    y = np.asarray(y)

    # --- Ensure at least 2D ---
    if x.ndim == 1:
        x = x[None, :]  # (1, N)
    elif x.ndim != 2:
        raise ValueError("x must be 1D or 2D")

    if y.ndim == 1:
        y = y[None, :]  # (1, N)
    elif y.ndim != 2:
        raise ValueError("y must be 1D or 2D")

    # --- Ensure shape is (n_args, N) ---
    # If passed as (N, n_args), transpose
    if x.shape[0] != y.shape[0] and x.shape[1] == y.shape[0]:
        x = x.T
        y = y.T

    n_args = min(x.shape[0], y.shape[0])

    if len(label_list) < n_args:
        raise ValueError("label_list shorter than number of curves")

    # --- Plot ---
    for k in range(n_args):
        plt.plot(x[k], y[k], linewidth=2, label=label_list[k])

    # --- Orbital Period Lines ---
    if 'period' in globals() and period is not None and len(orbital_lines) > 0:
        for lx in orbital_lines:
            ax.axvline(lx, color='black', linestyle='-.', alpha=0.3, zorder=0)
            ax.text(lx, 1.05, f"{int(lx / period)} T ≈ {int(lx)} s",
                    transform=ax.get_xaxis_transform(),
                    ha='center', va='bottom', fontsize=12, color='gray')

    plt.xlabel(xlabel, fontsize=LABEL_SIZE)
    plt.ylabel(ylabel, fontsize=LABEL_SIZE)

    # Push title way up
    plt.title(title, fontsize=TITLE_SIZE, y=1.20)
    plt.grid(True, linestyle=':', alpha=0.6)

    ax.tick_params(labelsize=TICK_SIZE)


    # --- ADD ORBITAL PERIOD LABEL ---
    handles, labels = ax.get_legend_handles_labels()
    handles.append(Line2D([0], [0], color='none', label='1T = one orbital period'))
    labels.append('1T = one orbital period')
    # Legend ABOVE axis, sitting between title and plot
    ax.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, 1.08),
              ncol=2, fontsize=LEGEND_SIZE)

    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close(fig)

    return str(filename)


def plot_log(x: np.ndarray, y: np.ndarray, xlabel: str, ylabel: str, label_list: list,
             title="Generic title", filename="Generic_filename.png", scatter=False, x_inverted=False):
    """
    Generic semi-log y plot for multiple datasets.

    Args:
        x: (#args, N) array (x-vector)
        y: (#args, N) array (y-vector)
        label_list: List of label strings
        filename: Figure filename
    """
    fig, ax = plt.subplots(figsize=(8, 7))  # Increased height
    fig.subplots_adjust(top=0.80)  # Room for legend

    # Handle cases where input might be 1D or 2D
    if x.ndim == 1:
        x = x.reshape(1, -1)
    if y.ndim == 1:
        y = y.reshape(1, -1)

    n_args = min(y.shape[0], y.shape[1])

    # If shape is (N, args), transpose to (args, N)
    if y.shape[0] > y.shape[1]:
        x = x.T
        y = y.T
        n_args = y.shape[0]

    for k in range(n_args):
        if scatter:
            if k & 1:
                plt.plot(x[k, :], y[k, :], '-*', linewidth=1.2, label=label_list[k])
            else:
                plt.plot(x[k, :], y[k, :], '--*', linewidth=1.5, label=label_list[k])
        else:
            if k & 1:
                plt.plot(x[k, :], y[k, :], '-', linewidth=1.2, label=label_list[k])
            else:
                plt.plot(x[k, :], y[k, :], '--', linewidth=1.5, label=label_list[k])

    plt.xlabel(xlabel, fontsize=LABEL_SIZE)
    plt.ylabel(ylabel, fontsize=LABEL_SIZE)
    plt.title(title, fontsize=TITLE_SIZE, y=1.20)

    # Minor ticks between powers of 10
    ax.set_yscale('log')
    ax.set_xscale('log')
    # Configure major and minor log ticks
    ax.yaxis.set_major_locator(LogLocator(base=10))
    ax.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10)))

    # Draw major and minor grids
    ax.grid(True, which='major', linewidth=0.8)
    ax.grid(True, which='minor', linestyle=':', linewidth=0.5)
    plt.gca().yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10)))
    plt.gca().xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10)))

    # plt.grid(True, which='major', linewidth=0.8)
    # plt.grid(True, which='minor', linestyle=':', linewidth=0.5)
    plt.minorticks_on()
    # ax.set_yscale('log')
    if x_inverted:
        plt.gca().invert_xaxis()

    ax.tick_params(labelsize=TICK_SIZE)

    # --- ADD ORBITAL PERIOD LABEL ---
    handles, labels = ax.get_legend_handles_labels()
    handles.append(Line2D([0], [0], color='none', label='1T = one orbital period'))
    labels.append('1T = one orbital period')

    # Legend ABOVE axis
    ax.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, 1.08),
              ncol=2, fontsize=LEGEND_SIZE)

    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close()
    return str(filename)


# ----------------------------------------------------- A4 ----------------------------------------------------------- #


# def plot_kep_el_and_diff(kep_el_u, kep_el_d, kep_el_diff,
#                          time,
#                          label1, ylabel1,
#                          label2, ylabel2,
#                          label3, ylabel3,
#                          title="Orbital Elements and difference", filename="plot_3_axes.png"):
#     for i in range(0, 5):
#         plot_three_scales(time,
#                           kep_el_u[:, i], label1[i], ylabel1[i],
#                           kep_el_d[:, i], label2[i], ylabel2[i],
#                           kep_el_diff[:, i], label3[i], ylabel3[i],
#                           title[i], filename[i]
#                           )
#     return 0
#


def plot_comparison_residual(time,
                             data_u, label_u,
                             data_d, label_d,
                             data_diff, label_diff,
                             ylabel_main, ylabel_diff,
                             title="Comparison", filename="comp_plot.png"):
    """
    Plots Disturbed and Undisturbed on Left Axis, Residual on Right Axis.
    """
    # --- Setup & Figure Size ---
    max_time = np.nanmax(np.array(time)) if len(time) > 0 else 0
    orbital_lines = []

    # Calculate lines
    if period is not None and max_time > period:
        orbital_lines = np.arange(period, max_time + 1, period, dtype=float)
        fig, ax1 = plt.subplots(figsize=(8 + len(orbital_lines), 8))
    else:
        fig, ax1 = plt.subplots(figsize=(8, 7))

    fig.subplots_adjust(top=0.82)

    # Colors matching plot_three_scales
    color_u = 'tab:blue'
    color_d = 'tab:green'
    color_diff = 'tab:orange'  # Changed to orange to match 3rd axis of plot_three_scales

    # --- Axis 1 (Left) - Main Values (Unperturbed & Perturbed) ---
    # Plot Undisturbed
    l1, = ax1.plot(time, data_u, color=color_u, label=label_u, linestyle='--', linewidth=5, alpha=1)
    # Plot Disturbed
    l2, = ax1.plot(time, data_d, color=color_d, label=label_d, linestyle='-', linewidth=7, alpha=0.5)

    ax1.set_xlabel("Time [s]", fontsize=LABEL_SIZE)
    ax1.set_ylabel(ylabel_main, fontsize=LABEL_SIZE)
    ax1.tick_params(axis='y', labelsize=TICK_SIZE)
    ax1.tick_params(axis='x', labelsize=TICK_SIZE)
    ax1.grid(True, which='major', linestyle=':', alpha=0.6)

    # Move Scientific Notation (Offset Text) FAR LEFT
    ax1.yaxis.get_offset_text().set_horizontalalignment('right')
    ax1.yaxis.get_offset_text().set_position((-0.15, 1.0))
    ax1.xaxis.set_minor_locator(mticker.AutoMinorLocator())
    ax1.yaxis.set_minor_locator(mticker.AutoMinorLocator())

    # --- Axis 2 (Right) - Residual ---
    ax2 = ax1.twinx()
    l3, = ax2.plot(time, data_diff, color=color_diff, label=label_diff, linestyle='-.', linewidth=6)

    ax2.set_ylabel(ylabel_diff, color=color_diff, fontsize=LABEL_SIZE)
    ax2.tick_params(axis='y', labelcolor=color_diff, labelsize=TICK_SIZE)

    # Configure Right Axis (Scientific notation if needed)
    formatter = mticker.ScalarFormatter(useMathText=True)
    formatter.set_powerlimits((-3, 2))
    ax2.yaxis.set_major_formatter(formatter)
    ax2.yaxis.get_offset_text().set_x(1.1)
    ax2.yaxis.set_minor_locator(mticker.AutoMinorLocator())

    # --- Orbital Period Lines ---
    if 'period' in globals() and period is not None and len(orbital_lines) > 0:
        for lx in orbital_lines:
            ax1.axvline(lx, color='black', linestyle='-.', alpha=0.4, zorder=0)
            ax1.text(lx, 1.02, f"{int(lx / period)} T ≈ {int(period * (lx / period))} [s]",
                     transform=ax1.get_xaxis_transform(),
                     ha='center', va='bottom', fontsize=12, color='gray')

    # --- Legend & Title ---
    lines = [l1, l2, l3]
    labels = [l.get_label() for l in lines]

    # --- ADD ORBITAL PERIOD LABEL ---
    lines.append(Line2D([0], [0], color='none', label='1T = one orbital period'))
    labels.append('1T = one orbital period')

    # Legend at Top
    ax1.legend(lines, labels, loc='lower center', bbox_to_anchor=(0.5, 1.05),
               ncol=2, borderaxespad=0., frameon=True, fontsize=LEGEND_SIZE)

    ax1.set_title(title, y=1.18, fontsize=TITLE_SIZE)

    # --- Save ---
    try:
        fig.savefig(filename, bbox_inches="tight")
    except ValueError:
        print(f"Warning: Infinite data detected in {filename}. Saving without 'tight' layout.")
        fig.savefig(filename)

    plt.close(fig)
    return str(filename)



# def plot_kep_el_and_diff(kep_el_u, kep_el_d, kep_el_diff,
#                          time,
#                          label1, ylabel1,
#                          label2, ylabel2,
#                          label3, ylabel3,
#                          title, filename):
#     """
#     Generates 5 plots (one for each Keplerian element).
#     Inputs are arrays where columns correspond to [a, e, i, Omega, omega].
#     """
#
#     generated_files = []
#
#     # Indices map: 0:a, 1:e, 2:i, 3:Omega, 4:omega
#
#     for i in range(5):
#         # Extract data for the i-th element
#         u_data = kep_el_u[:, i]
#         d_data = kep_el_d[:, i]
#         diff_data = kep_el_diff[:, i]
#
#         # Handle Angle Wrapping for Omega (idx 3) and omega (idx 4)
#         if i ==4 :
#             u_data = mask_discontinuity(u_data)
#             d_data = mask_discontinuity(d_data)
#             diff_data = mask_discontinuity(diff_data)
#
#         # Handle filename list or single pattern
#         # If 'filename' is a list, pick i-th. If string, append _i
#         fname = filename[i] if isinstance(filename, list) else f"{filename}_{i}.png"
#         tit = title[i] if isinstance(title, list) else f"{title} - Element {i}"
#
#         # Call the plotter
#         out = plot_comparison_residual(
#             time,
#             u_data, label1[i],
#             d_data, label2[i],
#             diff_data, label3[i],
#             ylabel1[i], ylabel3[i],  # Main Y-label matches the elements
#             title=tit,
#             filename=fname
#         )
#         generated_files.append(out)
#
#     return generated_files



def plot_mag_acc(time: np.ndarray, acc_v: np.ndarray, xlabel: str, ylabel: str, label_str: str,
                 title="Generic title", filename="Generic_filename.png"):
    """
    Args:
        :param filename:
        :param title:
        :param label_str:
        :param ylabel:
        :param xlabel:
        :param acc_v:
        :param time:
    :returns: str(filename)
    """
    # --- 3. Figure Setup ---
    max_time = np.nanmax(time)
    orbital_lines = []
    # Calculate lines
    if period is not None and max_time > period:
        orbital_lines = np.arange(period, max_time + 1, period, dtype=float)
        fig, ax = plt.subplots(figsize=(8 + len(orbital_lines), 8))  # Increased height for legend
    else:
        fig, ax = plt.subplots(figsize=(8, 7))  # Increased height

    fig.subplots_adjust(top=0.82)  # Room for legend

    acc = errors.mag_dist(acc_v)

    plt.plot(time, acc, linewidth=2, label=label_str)  # Added label here for legend

    plt.xlabel(xlabel, fontsize=LABEL_SIZE)
    plt.ylabel(ylabel, fontsize=LABEL_SIZE)
    ax.set_title(title, fontsize=TITLE_SIZE, y=1.20)
    plt.grid(True, linestyle=':', alpha=0.6)

    ax.tick_params(labelsize=TICK_SIZE)
    # Turn on minor ticks
    ax.xaxis.set_minor_locator(mticker.AutoMinorLocator())
    ax.yaxis.set_minor_locator(mticker.AutoMinorLocator())
    # ax.legend(loc='best', fontsize=LEGEND_SIZE)

    # --- Orbital Period Lines ---
    if 'period' in globals() and period is not None and len(orbital_lines) > 0:
        for lx in orbital_lines:
            ax.axvline(lx, color='black', linestyle='-.', alpha=0.3, zorder=0)
            ax.text(lx, 1.05, f"{int(lx / period)} T ≈ {int(lx)} s",
                    transform=ax.get_xaxis_transform(),
                    ha='center', va='bottom', fontsize=12, color='gray')

    # --- ADD ORBITAL PERIOD LABEL ---
    handles, labels = ax.get_legend_handles_labels()
    handles.append(Line2D([0], [0], color='none', label='1T = one orbital period'))
    labels.append('1T = one orbital period')

    # Legend ABOVE axis
    ax.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, 1.08),
              ncol=2, fontsize=LEGEND_SIZE)

    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return str(filename)