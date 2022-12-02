import os, re, warnings
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
warnings.filterwarnings('ignore', category=RuntimeWarning)

def spacing_dist(x, a, alpha, b, beta):
    return a * x**alpha * np.exp(b * x**beta)

types = [re.findall(r"(\d+)_([A-Za-z]+)", file) for file in os.listdir("results/")]
types = [t[0] for t in types if len(t)>0]
types = {t[0]+"_"+t[1]:{"size":int(t[0]), "method":t[1]} for t in types}


##################################
##### CREATION OF ARRAYS #########
##################################
# fill datasets
hist_data = {}

for type_ in types:
    filename = "results/norm_spac_"+type_+".dat"
    hist_data[type_] = np.loadtxt(filename)


##################################
##### PLOTS AND FITS #############
##################################

nbins = 100

# plot styling
col_dict = {k:'#377eb8' if val['method']=='HS' else "#ff7f00" for k, val in types.items()}
label_dict = {k:'Hermitian, size '+str(val['size'])+', '+str(nbins)+' bins' if val['method']=='HS'
                else 'Diagonal real, size '+str(val['size'])+', '+str(nbins)+' bins' for k, val in types.items()}


fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 7))

# plot all the desider plots
for type_ in types:

    data = hist_data[type_][hist_data[type_] < 0.5*np.max(hist_data[type_])]
    probs, bin_edges, _ = ax.hist(data, density=True, bins=nbins, alpha=.5, color=col_dict[type_], edgecolor="black", label=label_dict[type_])
    bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])

    # fit
    popt, pcov = curve_fit(f=spacing_dist, xdata=bin_centers, ydata=probs, p0=[1.0, 1.0, -1.0, 1.0])
    print(popt)
    print()
    print(pcov)

    # plot
    x0 = np.linspace(min(data), max(data), 1000)
    ax.plot(x0, spacing_dist(x0, *popt), ls='solid', c=col_dict[type_], lw=1.5,
            label=rf"a = {popt[0]:.2f}, $\alpha$= {popt[1]:.2f}, b = {popt[2]:.2f}, $\beta$ = {popt[3]:.2f}")


    ax.set_title("Distribution of normalized spacings \nfor different random matrix types ", fontsize=30)
    ax.set_xlabel(r"Eigenvalue spacing ($s$)", fontsize=15)
    ax.set_ylabel(r"Probability density $P(s)$", fontsize=15)
    ax.tick_params(axis="both", which="major", labelsize=15)
    ax.grid(True, which="both", linestyle=":")
    ax.legend(fontsize=15)

# print formula
ax.text(0.5, 0.5, s=r'$P(s)=a\, s^{\alpha}\, \exp{(b s^{\beta})}$', transform=ax.transAxes, fontsize=25, bbox=dict(facecolor="w"))


# save image
fig.savefig("results/fit_norm_spacing.png", dpi=400, format='png')

