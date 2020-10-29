import numpy as np
import matplotlib.pyplot as plt
import cycler
from unsorted import quicksummary, analyse_performance


def extract_performance(run_ids,run_names):
    for run_namelist,run_id in zip(run_names,run_ids):
        for run_name in run_namelist:
            print("summarising performance:",run_id,run_name)
            analyse_performance.dump_performance(run_id, run_name)

def plot_performance(run_ids,run_names,colour_cycle):
    plt.figure()
    plt.rc('axes',prop_cycle=colour_cycle)

    for run_namelist,run_id in zip(run_names,run_ids):
        for run_name in run_namelist:
            print("plotting performance:",run_id,run_name)
            perf_data = np.loadtxt("data/info_plotable_"+run_id+run_name+".dat")
#             perf_data = np.random.random((100,2))
            plt.plot(perf_data[:,0],perf_data[:,1],label=run_name)

    plt.legend()
    plt.savefig("../figures/SN_performance_summary.png",dpi=200)
    plt.close()

def extract_summary(run_ids,run_names,ordered):
    for run_namelist,run_id,dumpsOrdered in zip(run_names,run_ids,ordered):
        for run_name in run_namelist:
            print("summarising energy:",run_id,run_name)
            quicksummary.summarise_and_dump(run_id, run_name, dumpsOrdered=dumpsOrdered)

def plot_energy(run_ids,run_names,colour_cycle):
    plt.figure(figsize=(9,9))
    plt.rc('axes',prop_cycle=colour_cycle)

    for run_namelist,run_id in zip(run_names,run_ids):
        for run_name in run_namelist:
            print("plotting energy:",run_id,run_name)
            energy_data = np.loadtxt("data/quicksummary"+run_id+run_name+".dat",skiprows=1)
#             perf_data = np.random.random((100,2))
            plt.plot(energy_data[:,0],energy_data[:,3]*energy_data[:,1],label=run_name)

    plt.legend()
    plt.savefig("../figures/SN_energy_summary.png",dpi=200)
    plt.close()

def compare_energy_time(run_ids,run_names):
    print("run_name dEnergy_dt dtime_dstep")
    for run_namelist,run_id in zip(run_names,run_ids):
        for run_name in run_namelist:
            energy_data = np.loadtxt("data/quicksummary"+run_id+run_name+".dat",skiprows=1)
            perf_data = np.loadtxt("data/info_plotable_"+run_id+run_name+".dat")
            dEnergy_dt = (energy_data[-1,3]-energy_data[0,3])/(energy_data[-1,0])
            dtime_dstep = (perf_data[-1,3]-perf_data[0,1])/(perf_data[-1,0])
            print(run_name,dEnergy_dt,dtime_dstep)

run_names = [
    ["SN_test_m1","SN_test_m2","SN_test_m3","SN_test_m3_thinner","SN_test_m3_thin","SN_test_T5_half",
    "SN_test_T5","SN_test_T6_half","SN_test_T6","SN_test_T7_half","SN_test_T7"],
    ["SN_test_m2_2","SN_test_m2_5","SN_test_m2_7","SN_test_m2_9","SN_test_m3_slowSN"]
    ]

run_ids = ["2030","3030"]

ordered=[True,False]

cmap = plt.cm.get_cmap('jet')

nlines = np.sum([len(x) for x in run_names])
colour_cycle = cycler.cycler('color',cmap(np.linspace(0.,1.,nlines)))


# extract_performance(run_ids,run_names)

# plot_performance(run_ids,run_names,colour_cycle)

# extract_summary(run_ids,run_names,ordered)

# plot_energy(run_ids,run_names,colour_cycle)

compare_energy_time(run_ids,run_names)




