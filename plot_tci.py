import os
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt

def get_cd(s10_max):
    cd = np.zeros_like(s10_max)

    cd1 = 1.0e-3
    cd2 = 2.4e-3
    s1 = 5.0
    s2 = 25.0

    for i, s10 in enumerate(s10_max):
        cd[i] = s10 * (cd2-cd1)/(s2-s1) + cd1
        cd[i] = min(cd[i], cd2)
        cd[i] = max(cd[i], cd1)

    return cd

def get_inertial_stability_at_rmax(v10, r, s10max_index):
    Im1 = np.zeros(v10.shape[0])
    Im2 = np.zeros(v10.shape[0])
    for i in range(v10.shape[0]):
        v10_i = v10[i,:]
        imax = s10max_index[i]
        Im1[i] = np.sqrt(f + 2*v10_i[imax]/r[imax])
        Im2[i] = np.sqrt(f + v10_i[imax]/r[imax] + (v10_i[imax+1] - v10_i[imax-1])/(r[imax+1] - r[imax-1]))
    
    Im = Im1 * Im2

    return Im, Im1, Im2

def get_efficiency_from_ir(s10max, cdmax, smpi, time):
    # Calculate the efficiency from intensification rate (decomposed)
    ir = (s10max[2:] - s10max[:-2])/(time[2:] - time[:-2])

    e_ir1 = ((H/cdmax[1:-1]) * ir)/(smpi**2)
    e_ir2 = (s10max[1:-1]**2)/(smpi**2)
    e_ir = e_ir1 + e_ir2

    return e_ir, e_ir1, e_ir2

def get_efficiency_from_is(v10, r, s10max_index, smpi_index):
    # Calculate the efficiency from inertial stability (decomposed)
    Im, Im1, Im2 = get_inertial_stability_at_rmax(v10, r, s10max_index)
    e1 = np.power((Im1/Im1[smpi_index]),n)
    e2 = np.power((Im2/Im2[smpi_index]),n)
    e = e1 * e2

    return e, e1, e2

def get_decompose_lc(v10, r, s10max_index, smpi_index):
    ratio_v10 = np.zeros(v10.shape[0])
    ratio_r = np.zeros(v10.shape[0])
    v10_mpi = v10[smpi_index,s10max_index[smpi_index]]
    r_mpi = r[s10max_index[smpi_index]]
    for i in range(v10.shape[0]):
        v10_i = v10[i,:]
        imax = s10max_index[i]
        ratio_v10[i] = v10_i[imax]/v10_mpi
        ratio_r[i] = r[imax]/r_mpi

    return ratio_v10, ratio_r, np.sqrt(ratio_v10/ratio_r)

def get_tcf(s10, r):
    tcf = np.zeros(s10.shape[0])
    for i in range(s10.shape[0]):
        s10_i = s10[i,:]
        imax = np.argmax(s10_i)
        rmw = r[imax]
        r17 = 0
        if s10_i[imax] < 17:
            tcf[i] = np.nan
            continue
        else:
            for j in range(imax, len(s10_i) - 1):
                if s10_i[j] >= 17 > s10_i[j + 1]:
                    r17 = r[j] + (r[j + 1] - r[j]) * (17 - s10_i[j]) / (s10_i[j + 1] - s10_i[j])
                    break

        tcf[i] = 1 - rmw / r17

    return tcf

def calculate_results(v10_list, s10_list, r, time):
    # s10 max for all members
    s10_max_list = np.max(s10_list, axis=2)

    # s10 max for ensemble mean
    s10 = np.average(s10_list, axis=0)
    s10_max = np.max(s10, axis=1)

    # e from intensification rate
    smpi = np.max(s10_max)
    cdmax = get_cd(s10_max)
    e_ir, e_ir1, e_ir2 = get_efficiency_from_ir(s10_max, cdmax, smpi, time)

    # e from inertial stability
    v10 = np.average(v10_list, axis=0)
    s10max_index = np.argmax(s10,axis=1)    # radial index
    smpi_index = np.argmax(s10_max)         # time index
    e_is, e_is1, e_is2 = get_efficiency_from_is(v10, r, s10max_index, smpi_index)

    # local coriolis term (without f)
    ratio_v10, ratio_r, lc = get_decompose_lc(v10, r, s10max_index, smpi_index)

    # TCF
    tcf = get_tcf(s10, r)

    return s10_max_list, s10_max, e_ir, e_is, e_ir1, e_ir2, e_is1, e_is2, ratio_v10, ratio_r, lc, tcf

def plot_results(s10_max_list, s10_max, e_ir, e_is, e_ir1, e_ir2, e_is1, e_is2, ratio_v10, ratio_r, lc, tcf, time, fig_dir):
    plt.figure()
    for s10_max_line in s10_max_list:
        plt.plot(time, s10_max_line, color='#D3D3D3', linewidth=0.5)
    plt.plot(time, s10_max, linewidth=2)
    plt.xlabel('Time(days)')
    plt.ylabel('10m wind speed (m/s)')
    plt.grid(True)
    fig_path_all = os.path.join(fig_dir, 's10_max_all.png')
    plt.savefig(fig_path_all)

    plt.figure()
    plt.plot(time[1:-1], e_ir, label='e_ir')
    plt.plot(time, e_is, label='e_is')
    plt.xlabel('Time (days)')
    plt.ylabel('Efficiency')
    plt.legend()
    plt.grid(True)
    fig_path_efficiency = os.path.join(fig_dir, 'efficiency.png')
    plt.savefig(fig_path_efficiency)

    plt.figure()
    plt.plot(time[1:-1], e_ir1, label='Intensification Rate') 
    plt.plot(time[1:-1], e_ir2, label='Vmax')
    plt.xlabel('Time (days)')
    plt.ylabel('Efficiency')
    plt.legend()
    plt.grid(True)
    plt.title('Decomposition of e_ir')
    fig_path_efficiency = os.path.join(fig_dir, 'e_ir_decomposition.png')
    plt.savefig(fig_path_efficiency)

    plt.figure()
    plt.plot(time, e_is1, label='Local Coriolis')
    plt.plot(time, e_is2, label='Vorticity')
    plt.xlabel('Time (days)')
    plt.ylabel('Efficiency')
    plt.legend()
    plt.grid(True)
    plt.title('Decomoosition of e_is')
    fig_path_efficiency = os.path.join(fig_dir, 'e_is_decomposition.png')
    plt.savefig(fig_path_efficiency)

    fig, ax1 = plt.subplots()
    color = 'tab:red'
    ax1.set_xlabel('Time (days)')
    ax1.set_ylabel('Ratio', color=color)
    ax1.plot(time, ratio_v10, label='v_max/v_mpi', color=color, linestyle='-')
    ax1.plot(time, ratio_r, label='r_max/r_mpi', color=color, linestyle='--')
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.legend(loc='upper left')
    ax1.grid(True)
    ax2 = ax1.twinx()
    color = 'tab:blue'
    ax2.set_ylabel('Local Coriolis', color=color)
    ax2.plot(time, lc, label='LC approximation', color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.legend(loc='upper right')
    fig.tight_layout()
    plt.title('Decomposition of LC')
    plt.gcf().subplots_adjust(left=0.10,top=0.90,bottom=0.10) 
    fig_path_lc = os.path.join(fig_dir, 'lc_decomposition.png')
    plt.savefig(fig_path_lc)

    plt.figure()
    plt.plot(time, tcf)
    plt.xlabel('Time (days)')
    plt.ylabel('TCF')
    plt.grid(True)
    plt.title('TCF')
    fig_path_tcf = os.path.join(fig_dir, 'tcf.png')
    plt.savefig(fig_path_tcf)

def main():
    # Default parameters
    global H, n, f
    H = 2000.0
    n = 1.0
    f = 5e-5

    ouput_file = 'ensemble_output.nc'
    fig_dir = './figures'

    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)

    with Dataset(ouput_file, 'r') as nc_file:
        v10_list = np.squeeze(nc_file.variables['v10_list'][:])
        s10_list = np.squeeze(nc_file.variables['s10_list'][:])
        r = np.squeeze(nc_file.variables['r'][:])
        time = np.squeeze(nc_file.variables['time'][:])

    results = calculate_results(v10_list, s10_list, r, time)

    plot_results(*results, time/86400, fig_dir)


if __name__ == "__main__":
    main()