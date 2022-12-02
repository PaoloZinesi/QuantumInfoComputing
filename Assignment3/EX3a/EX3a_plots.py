import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit

types = ['naive', 'opt', 'builtin']

def pow_func(x, ln_const, pow):
    return ln_const + pow * x


##################################
##### CREATION OF DATAFRAMES #####
##################################
time_dfs = pd.DataFrame({"size":[], "type":[], "time":[]})

# fill datasets
for type_ in types:
    file = "results/"+type_+"_performances.dat"

    # add to existing dataframe
    df = pd.read_csv(file, header=None, names=['rowA','colA','rowB','colB','time'], skiprows=0, delimiter=";")\
        .astype({'rowA':'int32','colA':'int32','rowB':'int32','colB':'int32'})
    df['size'] = df.apply(lambda x: (x.rowA*x.colA*x.rowB*x.colB)**0.25, axis=1)
    df['type'] = type_
    time_dfs = pd.concat([time_dfs, df], axis=0, ignore_index=True)[["size","type", "time"]]



# groupby repetitions of the same configuration
# then take the average, std, and rename columns 
agg_times = time_dfs.groupby(["size","type"]).agg(["mean", "std"])
agg_times.columns = ['_'.join(col) for col in agg_times.columns.to_flat_index().values]
agg_times.reset_index(inplace=True)




##################################
##### PLOTS AND FITS #############
##################################

# plot styling
col_dict = {"builtin":'#377eb8', "opt":'#ff7f00', "naive":'#4daf4a'}
label_dict = {"builtin":"Intrinsic", "opt":"Optimized", "naive":"Naive"}


fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 8))

# plot all the desider plots
for type_ in types:
    plot_df = agg_times[(agg_times["type"]==type_) & (agg_times["size"]>=100)]
    ax.errorbar(x="size", y="time_mean", yerr="time_std", data=plot_df,
                fmt="o", label=f"{label_dict[type_]} function with -O3 flag",
                capsize=5, c=col_dict[type_])

    # fit
    if type_=='naive':
        # restrict domain
        popt, pcov = curve_fit(f=pow_func, xdata=np.log(plot_df.loc[plot_df['size']>=800, 'size']),
                               ydata=np.log(plot_df.loc[plot_df['size']>=800, 'time_mean']))
    else:
        popt, pcov = curve_fit(f=pow_func, xdata=np.log(plot_df['size']), ydata=np.log(plot_df['time_mean']))
    ln_const, power = popt

    # plot
    x0 = np.linspace(min(plot_df['size']), max(plot_df['size']), 1000).reshape(-1,1)
    ax.plot(x0, np.exp(pow_func(np.log(x0), ln_const=ln_const, pow=power)), ls='dotted', c=col_dict[type_],
            label=f"ln(const) = {ln_const:.2f}, power = {power:.2f}")
    

# set plot options
ax.set_title("Performances of matrix multiplication algorithms", fontsize=30)
ax.set_xlabel("Dimension of the square matrices", fontsize=15)
ax.set_ylabel("Execution time [seconds]", fontsize=15)
ax.tick_params(axis="both", which="major", labelsize=15)
ax.grid(True, which="both", linestyle=":")
ax.set_yscale("log")
ax.set_xscale("log")
ax.legend(fontsize=12)

# print formula
ax.text(0.6, 0.2, s=r'time = const $\cdot x^{pow}$', transform=ax.transAxes, fontsize=25, bbox=dict(facecolor="w"))

# save image
fig.savefig("results/fit_matmul.png", dpi=400, format='png')