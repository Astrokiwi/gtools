import time_series_plot
import extract_sfr
import quickangle

# run_names = ["aniso_4","aniso_4_bigrad","aniso_4_biggerrad"]
run_names = ["superhr","hyperhr"]
run_ids = ['2020']*len(run_names)
outfile = "../figures/hr_time_series_summary.png"

extract_sfr.extract_sfr(run_ids,run_names)
quickangle.dump_full_evolution(run_ids,run_names)
time_series_plot.plot(run_ids,run_names,outfile)
