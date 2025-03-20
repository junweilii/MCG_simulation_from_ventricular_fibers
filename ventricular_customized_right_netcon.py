# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 10:12:50 2024

@author: lijun
"""
from neuron import h
from neuron.units import ms, mV
import matplotlib.pyplot as plt
import scipy.io as sio
import os
import numpy as np
import time
import csv
import pickle
from numpy import savetxt
import biot_savart as bs
from scipy.ndimage import gaussian_filter1d

current_dir = os.getcwd()
h.nrn_load_dll(current_dir + "\\beelerReuter\\nrnmech.dll")

data = sio.loadmat(current_dir + "\\data\\customized_right_ventricular_fiber_half_final.mat")
x_spiral_interp = 1000 * data['x_spiral_interp']  # convert mm to um
y_spiral_interp = 1000 * data['y_spiral_interp']
z_spiral_interp = 1000 * data['z_spiral_interp']
x_spiral_key = 1000 * data['x_spiral_key_nan']
y_spiral_key = 1000 * data['y_spiral_key_nan']
z_spiral_key = 1000 * data['z_spiral_key_nan']
zero_angle_ind = data['zero_angle_ind']
layer_num = np.size(x_spiral_interp, 0)
latitude_num = np.size(x_spiral_key[0][0], 1)
fiber_num = np.size(x_spiral_key[0][0], 0)

cell_length = 100  # um
gap_length = 0.008  # um, gap width 80 A, neglected
cell_diameter = 16  # um
Ra = 1 / 0.005  # myoplasm resistivity, Ω.cm
gap_G = 0.01  # unit: uS

node = [None] * layer_num
syn_fiber_left = [None] * layer_num
nc_fiber_left = [None] * layer_num
syn_fiber_right = [None] * layer_num
nc_fiber_right = [None] * layer_num
node_ind = 0
for layer_ind in range(layer_num):
    node[layer_ind] = [None] * fiber_num
    syn_fiber_left[layer_ind] = [None] * fiber_num
    nc_fiber_left[layer_ind] = [None] * fiber_num
    syn_fiber_right[layer_ind] = [None] * fiber_num
    nc_fiber_right[layer_ind] = [None] * fiber_num
    for fiber_i in range(fiber_num):
        # a fiber consists of 5~1000 cells of 100um length
        cell_num = np.size(x_spiral_interp[layer_ind][0][fiber_i][0], 1) - 1
        node[layer_ind][fiber_i] = [None] * cell_num

        for cell_ind in range(cell_num):
            if ~np.isnan(x_spiral_interp[layer_ind][0][fiber_i][0][0, cell_ind]):
                node_name = "node" + str(node_ind)
                node_ind = node_ind + 1
                node[layer_ind][fiber_i][cell_ind] = h.Section(name=node_name)  # Construct a Section.
                node[layer_ind][fiber_i][cell_ind].nseg = 3  # Set the number of segments.
                node[layer_ind][fiber_i][cell_ind].cm = 1  # membrane capacitance, uF/cm2
                node[layer_ind][fiber_i][cell_ind].Ra = Ra
                node[layer_ind][fiber_i][cell_ind].insert("Cadynam")  # Insert a mechanism.
                node[layer_ind][fiber_i][cell_ind].insert("IK1")
                node[layer_ind][fiber_i][cell_ind].insert("INa")
                node[layer_ind][fiber_i][cell_ind].insert("Is")
                node[layer_ind][fiber_i][cell_ind].insert("IKx1")
                node[layer_ind][fiber_i][cell_ind].insert("K_acc")
                node[layer_ind][fiber_i][cell_ind].insert("Na_acc")
                node[layer_ind][fiber_i][cell_ind].ena = 51.8691
                node[layer_ind][fiber_i][cell_ind].ek = -87.2343
                node[layer_ind][fiber_i][cell_ind].gsbar_Is = 0.00005  # very important setting
                if layer_ind == 0:
                    node[layer_ind][fiber_i][cell_ind].gK1_IK1 = 0.000345
                elif layer_ind == 1:
                    node[layer_ind][fiber_i][cell_ind].gK1_IK1 = 0.00035
                elif layer_ind == 2:
                    node[layer_ind][fiber_i][cell_ind].gK1_IK1 = 0.000355
                elif layer_ind == 3:
                    node[layer_ind][fiber_i][cell_ind].gK1_IK1 = 0.00036
                elif layer_ind == 4:
                    node[layer_ind][fiber_i][cell_ind].gK1_IK1 = 0.000365
                elif layer_ind == 5:
                    node[layer_ind][fiber_i][cell_ind].gK1_IK1 = 0.00037
                elif layer_ind == 6:
                    node[layer_ind][fiber_i][cell_ind].gK1_IK1 = 0.000375
                elif layer_ind == 7:
                    node[layer_ind][fiber_i][cell_ind].gK1_IK1 = 0.00038
                elif layer_ind == 8:
                    node[layer_ind][fiber_i][cell_ind].gK1_IK1 = 0.000385
                elif layer_ind == 9:
                    node[layer_ind][fiber_i][cell_ind].gK1_IK1 = 0.00039
                elif layer_ind == 10:
                    node[layer_ind][fiber_i][cell_ind].gK1_IK1 = 0.000395
                elif layer_ind == 11:
                    node[layer_ind][fiber_i][cell_ind].gK1_IK1 = 0.00040
                elif layer_ind == 12:
                    node[layer_ind][fiber_i][cell_ind].gK1_IK1 = 0.000405
                elif layer_ind == 13:
                    node[layer_ind][fiber_i][cell_ind].gK1_IK1 = 0.00041
                elif layer_ind == 14:
                    node[layer_ind][fiber_i][cell_ind].gK1_IK1 = 0.000415

                # node[cell_ind].psection()
                node[layer_ind][fiber_i][cell_ind].pt3dclear()
                node[layer_ind][fiber_i][cell_ind].pt3dadd(x_spiral_interp[layer_ind][0][fiber_i][0][0, cell_ind],
                                                           y_spiral_interp[layer_ind][0][fiber_i][0][0, cell_ind],
                                                           z_spiral_interp[layer_ind][0][fiber_i][0][0, cell_ind],
                                                           cell_diameter)
                node[layer_ind][fiber_i][cell_ind].pt3dadd(x_spiral_interp[layer_ind][0][fiber_i][0][0, cell_ind+1],
                                                           y_spiral_interp[layer_ind][0][fiber_i][0][0, cell_ind+1],
                                                           z_spiral_interp[layer_ind][0][fiber_i][0][0, cell_ind+1],
                                                           cell_diameter)

                # # watch the length of every cell
                # cell_real_L = node[layer_ind][fiber_i][cell_ind].L

        # # seg.area() returns the surface area of the cell segment instead of cross-sectional area
        # # ri() is calculated by Ra(), the segment length and segment cross-sectional area
        # for seg in node[0][0][0].allseg():
        #     print(seg, seg.area(), seg.ri())

        # ps = h.Shape(True)
        # xmin = min(x_spiral_interp[layer_ind][0][fiber_i, :])
        # xmax = max(x_spiral_interp[layer_ind][0][fiber_i, :])
        # ymin = min(y_spiral_interp[layer_ind][0][fiber_i, :])
        # ymax = max(y_spiral_interp[layer_ind][0][fiber_i, :])
        # ps.size(xmin, xmax, ymin, ymax)
        # ps.show(0)

        syn_fiber_left[layer_ind][fiber_i] = [None] * (cell_num - 1)
        nc_fiber_left[layer_ind][fiber_i] = [None] * (cell_num - 1)
        syn_fiber_right[layer_ind][fiber_i] = [None] * (cell_num - 1)
        nc_fiber_right[layer_ind][fiber_i] = [None] * (cell_num - 1)
        for cell_j in range(cell_num - 1):
            if (~np.isnan(x_spiral_interp[layer_ind][0][fiber_i][0][0, cell_j]) and
                    ~np.isnan(x_spiral_interp[layer_ind][0][fiber_i][0][0, cell_j+1])):
                # define two symmetric synapse to mimic gap junction
                syn_fiber_left[layer_ind][fiber_i][cell_j] = h.ExpSyn(node[layer_ind][fiber_i][cell_j + 1](0))
                syn_fiber_left[layer_ind][fiber_i][cell_j].tau = 3  # ms, decay time constant
                syn_fiber_left[layer_ind][fiber_i][cell_j].e = 0
                nc_fiber_left[layer_ind][fiber_i][cell_j] = h.NetCon(node[layer_ind][fiber_i][cell_j](1)._ref_v, syn_fiber_left[layer_ind][fiber_i][cell_j], sec=node[layer_ind][fiber_i][cell_j])
                nc_fiber_left[layer_ind][fiber_i][cell_j].weight[0] = 1
                nc_fiber_left[layer_ind][fiber_i][cell_j].delay = 0  # set here to adjust node delay!!!!
                nc_fiber_left[layer_ind][fiber_i][cell_j].threshold = -30
                # symmetric
                syn_fiber_right[layer_ind][fiber_i][cell_j] = h.ExpSyn(node[layer_ind][fiber_i][cell_j](1))
                syn_fiber_right[layer_ind][fiber_i][cell_j].tau = 3  # ms, decay time constant
                syn_fiber_right[layer_ind][fiber_i][cell_j].e = 0
                nc_fiber_right[layer_ind][fiber_i][cell_j] = h.NetCon(node[layer_ind][fiber_i][cell_j+1](0)._ref_v,
                                                                     syn_fiber_right[layer_ind][fiber_i][cell_j],
                                                                     sec=node[layer_ind][fiber_i][cell_j+1])
                nc_fiber_right[layer_ind][fiber_i][cell_j].weight[0] = 1
                nc_fiber_right[layer_ind][fiber_i][cell_j].delay = 0  # set here to adjust node delay!!!!
                nc_fiber_right[layer_ind][fiber_i][cell_j].threshold = -30
    print(layer_ind)

    # h.topology()
# fig = plt.figure(figsize=(8, 8))
# ax = plt.axes(projection='3d')
# for layer_ind in range(layer_num):
#     for fiber_i in range(fiber_num):
#         ax.plot3D(x_spiral_key[layer_ind][0][fiber_i, :],
#                   y_spiral_key[layer_ind][0][fiber_i, :],
#                   z_spiral_key[layer_ind][0][fiber_i, :], linewidth=0.5, c='green', alpha=0.3)

nc_in = [None] * fiber_num
syn_in = [None] * fiber_num
nc_out = [None] * fiber_num
syn_out = [None] * fiber_num
gap_fiber_ind = [None] * fiber_num
gap_cell_ind = [None] * fiber_num
for fiber_i in range(fiber_num):
    nc_in[fiber_i] = [None] * (latitude_num - 1)
    syn_in[fiber_i] = [None] * (latitude_num - 1)
    nc_out[fiber_i] = [None] * (latitude_num - 1)
    syn_out[fiber_i] = [None] * (latitude_num - 1)
    gap_fiber_ind[fiber_i] = [None] * (latitude_num - 1)
    gap_cell_ind[fiber_i] = [None] * (latitude_num - 1)
    # because gap is in cell center, so the last latitude cell do not have gap
    for latitude_ind in range(latitude_num - 1):
        nc_in[fiber_i][latitude_ind] = [None] * layer_num
        syn_in[fiber_i][latitude_ind] = [None] * layer_num
        nc_out[fiber_i][latitude_ind] = [None] * layer_num
        syn_out[fiber_i][latitude_ind] = [None] * layer_num
        gap_fiber_ind[fiber_i][latitude_ind] = [None] * layer_num
        gap_cell_ind[fiber_i][latitude_ind] = [None] * layer_num
        gap_point = np.zeros((layer_num, 3))
        for layer_ind in range(layer_num):
            if zero_angle_ind[layer_ind, latitude_ind] + fiber_i + 1 <= fiber_num:
                gap_fiber_ind[fiber_i][latitude_ind][layer_ind] = zero_angle_ind[layer_ind, latitude_ind] + fiber_i
            else:
                gap_fiber_ind[fiber_i][latitude_ind][layer_ind] = zero_angle_ind[layer_ind, latitude_ind] + fiber_i - fiber_num

            if ~np.isnan(x_spiral_key[layer_ind][0][gap_fiber_ind[fiber_i][latitude_ind][layer_ind], latitude_ind]):
                gap_point[layer_ind, :] = [x_spiral_key[layer_ind][0]
                                           [gap_fiber_ind[fiber_i][latitude_ind][layer_ind], latitude_ind],
                                           y_spiral_key[layer_ind][0]
                                           [gap_fiber_ind[fiber_i][latitude_ind][layer_ind], latitude_ind],
                                           z_spiral_key[layer_ind][0]
                                           [gap_fiber_ind[fiber_i][latitude_ind][layer_ind], latitude_ind]]
                gap_point_z = gap_point[layer_ind, 2]
                # find the position of z_spiral_key point (unique) on z_spiral_interp
                temp = abs(z_spiral_interp[layer_ind][0][gap_fiber_ind[fiber_i][latitude_ind][layer_ind]][0][0, :] - gap_point_z)
                gap_cell_ind[fiber_i][latitude_ind][layer_ind] = np.where(temp == np.nanmin(temp))
                gap_cell_ind[fiber_i][latitude_ind][layer_ind] = gap_cell_ind[fiber_i][latitude_ind][layer_ind][0][0]
                if gap_cell_ind[fiber_i][latitude_ind][layer_ind] == len(temp)-1:
                    gap_cell_ind[fiber_i][latitude_ind][layer_ind] = gap_cell_ind[fiber_i][latitude_ind][layer_ind]-1

        # ax.plot3D(gap_point[:, 0], gap_point[:, 1], gap_point[:, 2], linewidth=0.5, c='red')
        for layer_ind in range(layer_num - 1):
            if (~np.isnan(x_spiral_key[layer_ind][0][gap_fiber_ind[fiber_i][latitude_ind][layer_ind], latitude_ind]) and
                    ~np.isnan(x_spiral_key[layer_ind+1][0][gap_fiber_ind[fiber_i][latitude_ind][layer_ind+1], latitude_ind])):
                syn_in[fiber_i][latitude_ind][layer_ind] = h.ExpSyn(
                    node[layer_ind + 1]
                    [gap_fiber_ind[fiber_i][latitude_ind][layer_ind + 1]]
                    [gap_cell_ind[fiber_i][latitude_ind][layer_ind + 1]](0.5))
                syn_in[fiber_i][latitude_ind][layer_ind].tau = 3
                syn_in[fiber_i][latitude_ind][layer_ind].e = 0

                nc_in[fiber_i][latitude_ind][layer_ind] = h.NetCon(
                    node[layer_ind]
                    [gap_fiber_ind[fiber_i][latitude_ind][layer_ind]]
                    [gap_cell_ind[fiber_i][latitude_ind][layer_ind]](0.5)._ref_v,
                    syn_in[fiber_i][latitude_ind][layer_ind], sec=node[layer_ind]
                    [gap_fiber_ind[fiber_i][latitude_ind][layer_ind]]
                    [gap_cell_ind[fiber_i][latitude_ind][layer_ind]])
                nc_in[fiber_i][latitude_ind][layer_ind].weight[0] = 1
                nc_in[fiber_i][latitude_ind][layer_ind].delay = 0
                nc_in[fiber_i][latitude_ind][layer_ind].threshold = -30

                syn_out[fiber_i][latitude_ind][layer_ind] = h.ExpSyn(
                    node[layer_ind]
                    [gap_fiber_ind[fiber_i][latitude_ind][layer_ind]]
                    [gap_cell_ind[fiber_i][latitude_ind][layer_ind]](0.5))
                syn_out[fiber_i][latitude_ind][layer_ind].tau = 3
                syn_out[fiber_i][latitude_ind][layer_ind].e = 0

                nc_out[fiber_i][latitude_ind][layer_ind] = h.NetCon(
                    node[layer_ind+1]
                    [gap_fiber_ind[fiber_i][latitude_ind][layer_ind+1]]
                    [gap_cell_ind[fiber_i][latitude_ind][layer_ind+1]](0.5)._ref_v,
                    syn_out[fiber_i][latitude_ind][layer_ind], sec=node[layer_ind+1]
                    [gap_fiber_ind[fiber_i][latitude_ind][layer_ind+1]]
                    [gap_cell_ind[fiber_i][latitude_ind][layer_ind+1]])
                nc_out[fiber_i][latitude_ind][layer_ind].weight[0] = 1
                nc_out[fiber_i][latitude_ind][layer_ind].delay = 0
                nc_out[fiber_i][latitude_ind][layer_ind].threshold = -30

# plt.show()

# define neural one spike
ns = h.NetStim()
ns.interval = 1000
ns.number = 1
unit_length = 50
# Track time with Vector t_vec.
t = h.Vector().record(h._ref_t)  # Time stamp vector
v050_node = [None] * layer_num
for layer_ind in range(layer_num):
    v050_node[layer_ind] = [None] * fiber_num
    if layer_ind == 0:
        syn = [None] * fiber_num
        nc = [None] * fiber_num
        stim = [None] * fiber_num
    for fiber_i in range(fiber_num):
        cell_num = np.size(x_spiral_interp[layer_ind][0][fiber_i][0][0, :], 0) - 1
        v050_node[layer_ind][fiber_i] = [None] * cell_num
        for cell_ind in range(cell_num):
            if ~np.isnan(x_spiral_interp[layer_ind][0][fiber_i][0][0, cell_ind]):
                v050_node[layer_ind][fiber_i][cell_ind] = (
                    h.Vector().record(node[layer_ind][fiber_i][cell_ind](0.5)._ref_v))

        if layer_ind == 0:  # only add stimulus on inner layer
            break_ind = list(range(0, cell_num // unit_length))  # every 50 cells are connected by gap junction
            break_ind = [i * unit_length + 25 for i in break_ind]
            unit_num = np.size(break_ind)
            syn[fiber_i] = [None] * unit_num
            nc[fiber_i] = [None] * unit_num
            stim[fiber_i] = [None] * unit_num
            for i in range(unit_num):
                syn_ind = break_ind[i]
                if ~np.isnan(x_spiral_interp[layer_ind][0][fiber_i][0][0, syn_ind]):
                    syn[fiber_i][i] = h.ExpSyn(node[layer_ind][fiber_i][syn_ind](0.25))
                    syn[fiber_i][i].tau = 3  # ms, decay time constant
                    syn[fiber_i][i].e = 0  # mV, reversal potential
                    nc[fiber_i][i] = h.NetCon(ns, syn[fiber_i][i])
                    nc[fiber_i][i].weight[0] = 1
                    nc[fiber_i][i].delay = 5
    print(layer_ind)
# Run simulation.
h.load_file("stdrun.hoc")
h.finitialize(-75 * mV)  # Initialize simulation with a voltage of -65 mV.
h.dt = 1  # 1ms or 0.5ms is not enough for cardiac cell, will lead to abnormal action potential
t_total = 500
visual_dt = 1
steps = int(t_total/visual_dt)
# 800ms simulation needs 10s computation, and a cycle of 800ms correspond to heart rate of 75 times/minute.
T1 = time.time()
for i in range(steps):
    h.continuerun(i * visual_dt * ms)  # advance to time t
    print("t=" + str(i * visual_dt) + "ms")

T2 = time.time()
print('程序运行时间:%s毫秒' % ((T2 - T1) * 1000))

# plt.figure()
# plt.plot(t, v025_node[0], 'b')
# plt.plot(t, v025_node[10], 'r')
# plt.plot(t, v025_node[cell_num-1], 'k')
# plt.title("Action potential")
# plt.xlabel("t (ms)")
# plt.ylabel("voltage (mV)")
# plt.show()

sigma = 2.5  # Standard deviation of the Gaussian kernel, 2.5 is proper, 5 is too smooth
radial_current = [None] * fiber_num
for fiber_i in range(fiber_num):
    radial_current[fiber_i] = [None] * (latitude_num - 1)
    # because gap is in cell center, so the last latitude cell do not have gap
    for latitude_ind in range(latitude_num - 1):
        radial_current[fiber_i][latitude_ind] = [None] * (layer_num-1)
        for layer_ind in range(layer_num-1):
            if (~np.isnan(x_spiral_key[layer_ind][0][gap_fiber_ind[fiber_i][latitude_ind][layer_ind], latitude_ind]) and
                    ~np.isnan(x_spiral_key[layer_ind + 1][0][gap_fiber_ind[fiber_i][latitude_ind][layer_ind + 1], latitude_ind])):
                radial_current[fiber_i][latitude_ind][layer_ind] = (
                    list(v050_node[layer_ind]
                    [gap_fiber_ind[fiber_i][latitude_ind][layer_ind]]
                    [gap_cell_ind[fiber_i][latitude_ind][layer_ind]]
                          - v050_node[layer_ind + 1]
                    [gap_fiber_ind[fiber_i][latitude_ind][layer_ind + 1]]
                    [gap_cell_ind[fiber_i][latitude_ind][layer_ind + 1]]))
                radial_current[fiber_i][latitude_ind][layer_ind] = 50 * 1e-3 * gaussian_filter1d(
                    radial_current[fiber_i][latitude_ind][layer_ind], sigma)
    print(fiber_i)

with open('radial_current_right_ventricle.pickle', 'wb') as handle:
    pickle.dump(radial_current, handle, protocol=pickle.HIGHEST_PROTOCOL)

fiber_current = [None] * layer_num
for layer_ind in range(layer_num):
    fiber_current[layer_ind] = [None] * fiber_num
    # because gap is in cell center, so the last latitude cell do not have gap
    for fiber_i in range(fiber_num):
        cell_num = np.size(x_spiral_interp[layer_ind][0][fiber_i][0][0, :], 0) - 1
        fiber_current[layer_ind][fiber_i] = [None] * (cell_num-1)
        for cell_ind in range(cell_num-1):
            if (~np.isnan(x_spiral_interp[layer_ind][0][fiber_i][0][0, cell_ind]) and
                    ~np.isnan(x_spiral_interp[layer_ind][0][fiber_i][0][0, cell_ind+1])):
                fiber_current[layer_ind][fiber_i][cell_ind] = (
                    list(v050_node[layer_ind][fiber_i][cell_ind]- v050_node[layer_ind][fiber_i][cell_ind+1]))
                fiber_current[layer_ind][fiber_i][cell_ind] = 50 * 1e-3 * gaussian_filter1d(
                    fiber_current[layer_ind][fiber_i][cell_ind], sigma)
    print(layer_ind)

with open('fiber_current_right_ventricle.pickle', 'wb') as handle:
    pickle.dump(fiber_current, handle, protocol=pickle.HIGHEST_PROTOCOL)

T1 = time.time()
t_num = len(t)
b_radial_sum = np.zeros((t_num, 3))
px_interest = 0
py_interest = -5
pz_interest = -2
for fiber_i in range(fiber_num):
    for latitude_ind in range(latitude_num - 1):
        heart_point = np.zeros((layer_num-1, 3))
        current = np.zeros((layer_num - 1, t_num))
        for layer_ind in range(layer_num-1):
            if (~np.isnan(x_spiral_key[layer_ind][0][gap_fiber_ind[fiber_i][latitude_ind][layer_ind], latitude_ind]) and
                    ~np.isnan(x_spiral_key[layer_ind + 1][0][gap_fiber_ind[fiber_i][latitude_ind][layer_ind + 1], latitude_ind])):
                x_gapJ = 0.5*(node[layer_ind]
                    [gap_fiber_ind[fiber_i][latitude_ind][layer_ind]]
                    [gap_cell_ind[fiber_i][latitude_ind][layer_ind]].x3d(0)+node[layer_ind]
                    [gap_fiber_ind[fiber_i][latitude_ind][layer_ind]]
                    [gap_cell_ind[fiber_i][latitude_ind][layer_ind]].x3d(1))
                y_gapJ = 0.5*(node[layer_ind]
                    [gap_fiber_ind[fiber_i][latitude_ind][layer_ind]]
                    [gap_cell_ind[fiber_i][latitude_ind][layer_ind]].y3d(0)+node[layer_ind]
                    [gap_fiber_ind[fiber_i][latitude_ind][layer_ind]]
                    [gap_cell_ind[fiber_i][latitude_ind][layer_ind]].y3d(1))
                z_gapJ = 0.5*(node[layer_ind]
                    [gap_fiber_ind[fiber_i][latitude_ind][layer_ind]]
                    [gap_cell_ind[fiber_i][latitude_ind][layer_ind]].z3d(0) + node[layer_ind]
                    [gap_fiber_ind[fiber_i][latitude_ind][layer_ind]]
                    [gap_cell_ind[fiber_i][latitude_ind][layer_ind]].z3d(1))
                heart_point[layer_ind, :] = [x_gapJ/10000, y_gapJ/10000, z_gapJ/10000]
                current[layer_ind, :] = radial_current[fiber_i][latitude_ind][layer_ind]
                heart_point = np.nan_to_num(heart_point)
                current = np.nan_to_num(current)
        b = []
        for t_ind in range(t_num):
            current_t = current[:, t_ind]
            fields = bs.produce_target_volume(heart_point, current_t, (0, 0, 0),
                                              (px_interest, py_interest, pz_interest), 1)
            b.append(np.reshape(fields[0][0][0], 3))

        b_radial_sum = b_radial_sum + np.array(b)
    print(fiber_i)

b_fiber_sum = np.zeros((t_num, 3))
for layer_ind in range(layer_num):
    for fiber_i in range(fiber_num):
        cell_num = np.size(x_spiral_interp[layer_ind][0][fiber_i][0][0, :], 0) - 1
        heart_point = np.zeros((cell_num - 1, 3))
        current = np.zeros((cell_num - 1, t_num))
        for cell_ind in range(cell_num-1):
            if hasattr(node[layer_ind][fiber_i][cell_ind], 'x3d'):
                x_gapJ = 0.5 * (node[layer_ind][fiber_i][cell_ind].x3d(0) + node[layer_ind][fiber_i][cell_ind].x3d(1))
                y_gapJ = 0.5 * (node[layer_ind][fiber_i][cell_ind].y3d(0) + node[layer_ind][fiber_i][cell_ind].y3d(1))
                z_gapJ = 0.5 * (node[layer_ind][fiber_i][cell_ind].z3d(0) + node[layer_ind][fiber_i][cell_ind].z3d(1))
                heart_point[cell_ind, :] = [x_gapJ / 10000, y_gapJ / 10000, z_gapJ / 10000]    # um to cm
                current[cell_ind, :] = fiber_current[layer_ind][fiber_i][cell_ind]
            heart_point = np.nan_to_num(heart_point)
            current = np.nan_to_num(current)
        b = []
        for t_ind in range(t_num):
            current_t = current[:, t_ind]
            fields = bs.produce_target_volume(heart_point, current_t, (0, 0, 0),
                                              (px_interest, py_interest, pz_interest), 1)
            b.append(np.reshape(fields[0][0][0], 3))
        b_fiber_sum = b_fiber_sum + np.array(b)
    print(layer_ind)

T2 = time.time()
print('程序运行时间:%s毫秒' % ((T2 - T1) * 1000))

plt.figure()
plt.plot(b_radial_sum[:, 0])
plt.plot(b_fiber_sum[:, 0])

from scipy.io import savemat
savemat('b_sum.mat',{'b_radial_sum':b_radial_sum,'b_fiber_sum':b_fiber_sum})

plt.figure()
plt.plot(v050_node[0][0][75])
plt.plot(v050_node[1]
         [gap_fiber_ind[0][1][1]]
         [gap_cell_ind[0][1][1]])

plt.figure()
for layer_ind in range(6):
    plt.plot(v050_node[layer_ind]
             [gap_fiber_ind[0][1][layer_ind]]
             [gap_cell_ind[0][1][layer_ind]],label='layer '+str(layer_ind))
