# Based on https://www.neuron.yale.edu/phpBB/viewtopic.php?f=2&t=2750

from neuron import *
import numpy as np
import sys
import matplotlib.pyplot as plt

duration = 500
dt = 0.01


def get_random_stim(rate):
    stimNc = h.NetStim()
    stimNc.noise = 1
    stimNc.start = 0
    stimNc.number = 1e9
    stimNc.interval = 1000.0 / rate
    return stimNc


def get_timed_stim():
    stimNc = h.NetStim()
    stimNc.noise = 0
    stimNc.start = 100
    stimNc.number = 3
    stimNc.interval = 30
    return stimNc


soma = h.Section()
soma.L = 17.841242
soma.diam = 17.841242
soma.push()
soma.insert("pas")
soma.g_pas = 0.0003
syn = h.ExpSyn(0.5, sec=soma)

stim1 = h.IClamp(0.5, sec=soma)
stim1.delay = 50.0
stim1.dur = 5.0
stim1.amp = 0.4

stim2 = h.IClamp(0.5, sec=soma)
stim2.delay = 140.0
stim2.dur = 5.0
stim2.amp = 0.4


def get_base_syn():
    syn = h.Syn4P(0.5, sec=soma)
    syn.A_LTD_pre = 0
    syn.A_LTP_pre = 0
    syn.A_LTD_post = 0
    syn.A_LTP_post = 0

    return syn


def get_ampa_syn():
    syn = get_base_syn()
    syn.s_ampa = 1
    syn.s_nmda = 0
    return syn


def get_preLTD_syn():
    syn = get_base_syn()
    syn.s_ampa = 1
    syn.s_nmda = 0
    syn.A_LTD_pre = 0.1
    return syn


def get_preLTP_syn():
    syn = get_base_syn()
    syn.s_ampa = 1
    syn.s_nmda = 0
    syn.A_LTP_pre = 0.1
    return syn


def get_postLTD_syn():
    syn = get_base_syn()
    syn.s_ampa = 1
    syn.s_nmda = 0
    syn.A_LTD_post = 0.1
    return syn


def get_nmda_syn():
    syn = get_base_syn()
    syn.s_ampa = 0
    syn.s_nmda = 1
    return syn


syn = get_base_syn()
syn = get_ampa_syn()
syn = get_preLTP_syn()
syn = get_postLTD_syn()


def run_sim(rate=10):
    print("Running simulation of %s with %s Hz input" % (duration, rate))

    stimNc = get_random_stim(rate)
    stimNc = get_timed_stim()

    vec_nc = h.Vector()

    nc = h.NetCon(stimNc, syn)
    nc.weight[0] = 0.001
    nc.delay = 0

    nc.record(vec_nc)

    vec = {}
    all_states = ["G", "Z", "E", "C", "P", "X", "u_bar", "T", "flag_D"]
    states = ["u_bar", "E", "T", "flag_D"]
    # states = all_states
    A_vals = ["A_LTD_pre", "A_LTP_pre", "A_LTD_post", "A_LTP_post"]
    A_vals = ["A_LTD_pre"]
    N_vals = ["N_alpha_bar", "N_alpha", "N_beta_bar", "N_beta", "N", "Z", "X"]
    postLTP_vals = ["G", "P", "C"]

    for var in (
        ["v", "t", "g", "g_ampa", "g_nmda", "w_pre", "w_post", "w"]
        + states
        + A_vals
        + N_vals
        + postLTP_vals
    ):
        vec[var] = h.Vector()
        if var != "v" and var != "t":
            exec("print('Recording:  %s')" % var)
            exec("vec ['%s'].record(syn._ref_%s)" % (var, var))

    # record the membrane potentials and
    # synaptic currents
    vec["v"].record(soma(0.5)._ref_v)
    vec["t"].record(h._ref_t)

    # run the simulation
    h.load_file("stdrun.hoc")
    h.init()
    h.tstop = duration
    h.dt = dt
    h.run()

    spikes = []
    isis = []
    lastSpike = None
    for t in vec_nc:
        spikes.append(t)
        # print(t)
        if lastSpike:
            isis.append(t - lastSpike)
        lastSpike = t

    hz = 1000 / (h.tstop / len(spikes))
    print("Spike times: %s" % ["%.3f" % t for t in vec_nc])
    print(
        "Num spikes: %s; avg rate: %s Hz; avg isi: %s ms; std isi: %s ms"
        % (len(spikes), hz, np.average(isis), np.std(isis))
    )
    # assert abs((hz-rate)/rate)<0.01

    scales = {
        "v": 0.001,
        "g": 1e-6,
        "u_bar": 1,
        "T": 1,
        "N_alpha_bar": 1,
        "N_beta_bar": 1,
        "N_alpha": 1,
        "N_beta": 1,
        "N": 1,
        "Z": 1,
        "X": 1,
        "w_pre": 1,
        "w_post": 1,
        "w": 1,
    }
    for v in postLTP_vals:
        scales[v] = 1
    for var in scales:
        var_file = open("%s.dat" % var, "w")
        for i in range(len(vec["t"])):
            var_file.write("%s\t%s\n" % (vec["t"][i] / 1000, vec[var][i] * scales[var]))
        var_file.close()

    if not "-nogui" in sys.argv:
        # plot the results
        plt.figure()
        plt.title("Membrane potential")
        plt.plot(vec["t"], vec["v"])

        plt.figure()
        plt.title("Conductance")
        plt.plot(vec["t"], vec["g"], label="g")
        plt.plot(vec["t"], vec["g_ampa"], label="g_ampa")
        plt.plot(vec["t"], vec["g_nmda"], label="g_nmda")
        plt.legend()

        plt.figure()
        plt.title("Weights")
        # plt.plot(vec['t'],vec['w'],label='w')
        plt.plot(vec["t"], vec["w_pre"], label="w_pre")
        plt.plot(vec["t"], vec["w_post"], label="w_post")
        plt.legend()

        plt.figure()
        plt.title("States")
        for s in states:
            plt.plot(vec["t"], vec[s], label=s)
        plt.legend()

        plt.figure()
        plt.title("A vals")
        for s in A_vals:
            plt.plot(vec["t"], vec[s], label=s)
        plt.legend()

        plt.figure()
        plt.title("N vals")
        for s in N_vals:
            plt.plot(vec["t"], vec[s], label=s)
        plt.legend()

        plt.figure()
        plt.title("postLTP vals")
        for s in postLTP_vals:
            plt.plot(vec["t"], vec[s], label=s)
        plt.legend()


run_sim()

if not "-nogui" in sys.argv:
    plt.show()
