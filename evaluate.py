# -*- coding: utf-8 -*-
"""
Created on Fri Nov 23 15:04:17 2018

@author: Saskia Philippen
"""

# import needed packages:
from __future__ import division
from numpy import *

import matplotlib
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'
matplotlib.rcParams['xtick.top'] = True
matplotlib.rcParams['ytick.right'] = True
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.ticker as ticker
from scipy.interpolate import UnivariateSpline, RectBivariateSpline
import glob
import pickle as pickle  

from scipy.stats import chi2, norm

def combinePartialResults(infile):
    dic = pickle.load(open("%s_0_mediumEnergyCut.pickle"%(infile.split(".pickle")[0]), 'rb'), encoding='latin1')
    logL = dic["logL"].T
    nSignal = dic["nSignal"].T
    logLBackground = [dic["logLBackground"]]
    for i in range(1,8):
        dic = pickle.load(open("%s_%d_mediumEnergyCut.pickle"%(infile.split(".pickle")[0], i), 'rb'), encoding='latin1')
        logL = concatenate((LLH, dic["logL"].T))
        nSignal = concatenate((nSignal, dic["nSignal"].T))
        logLBackground = concatenate((logLBackground, dic["logLBackground"].T))
    print(logLBackground)
    logLBackground = min(logLBackground)
    return nSignal, logL, logLBackground
#%%%    


def readResult(infile):
    if (isinstance(infile, list) and len(infile) == 1) or isinstance(infile, list)==False:
        print("Read in result from a single analysis:")
        if "*" in infile:
            print("combine partial sets to full maps...")
            nSignal, logL, logLBackground = combinePartialResultssort(glob.glob("%s.pickle"%(infile.split(".pickle")[0])))
        else:
            print("read in a full maps...")        
    elif isinstance(infile, list) and len(infile) != 1:
        
        
        if isinstance(infile, list): file = infile[0]
        elif isinstance(infile, str): file = infile
        else: 
            raise Exception("invalid input, infile has to be a single string or a list of strings")
        dic = pickle.load(open("%s.pickle"%(file.split(".pickle")[0]), 'rb'), encoding='latin1')
        logLBackground = dic["logLBackground"]      # background PDF Bi
        logL = dic["logL"].T
        nSignal = dic["nSignal"].T
        del dic
    return nSignal, logL, logLBackground


"""  
def readResult(infile):
    if (isinstance(infile, list) and len(infile) != 1):
        combinePartialResults(infile)   # TODO
        
    else:
        if isinstance(infile, list): file = infile[0]
        elif isinstance(infile, str): file = infile
        else: 
            raise Exception("invalid input, infile has to be a single string or a list of strings")
        dic = pickle.load(open("%s.pickle"%(file.split(".pickle")[0]), 'rb'), encoding='latin1')
        logLBackground = dic["logLBackground"]      # background PDF Bi
        logL = dic["logL"].T
        nSignal = dic["nSignal"].T
        del dic
    return nSignal, logL, logLBackground
"""
#%%%





def evaluate(infile, azRange=3, zenRange=3, plot=True, figSize=[4,3], 
             scanPoints=1601, moonPictureLocation="C:/Users/Saskia/Desktop/Uni/Master/Masterarbeit/MondFreigestellt.png"):
    img = plt.imread(moonPictureLocation)
    nSignal, logL, logLBackground = readResult(infile)
    
    # interpolate maps:
    gridPointsZen, gridPointsAz = shape(nSignal)
    azGrid = linspace(-azRange, azRange, gridPointsAz+1)
    azGrid = (azGrid[1:] + azGrid[:-1])/2
    zenGrid = linspace(-zenRange, zenRange, gridPointsZen+1)
    zenGrid = (zenGrid[1:] + zenGrid[:-1])/2
    azInterp = linspace(-azRange, azRange, scanPoints)
    zenInterp = linspace(-zenRange, zenRange, scanPoints)
    
    logL = RectBivariateSpline(azGrid, zenGrid, logL, s=0)
    logL = logL(azInterp, zenInterp)
    
    nSignal = RectBivariateSpline(azGrid, zenGrid, nSignal, s=0)
    nSignal = nSignal(azInterp, zenInterp)
    
    ############################
    ### CALCULATE QUANTITIES ###
    ############################   

    
    ### LLH ###
    
    
    
    # LLH minimum
    xlogL = (argmin(logL)/scanPoints - int(argmin(logL)/scanPoints))*scanPoints
    ylogL = int(argmin(logL)/scanPoints)
    xlogL = 2*azRange/scanPoints*xlogL - azRange
    ylogL = 2*zenRange/scanPoints*ylogL - zenRange
    
    ### Delta LLHs ###
    dLLH_ns0 = logL - logLBackground
    dLLH_min = logL - logL.min()
    
    ### nSignal ###
    
    # extrema nSignal
    xnSignal = (argmin(nSignal)/intersize - int(argmin(nSignal)/intersize))*intersize
    ynSignal = int(argmin(nSignal)/intersize)
    xnSignal = 6/intersize*xnSignal - 3
    ynSignal = 6/intersize*ynSignal - 3
    maximum = abs(nSignal).max()
    
    ### significance ###
    
    # significance 1: sigma(LLH - LLH(nSignal=0))
    dLLH_sig_ns0 = dLLH_ns0.copy()
    dLLH_sig_ns0[nSignal>0] = dLLH_sig_ns0[nSignal>0]*(-1)
    #sig_ns0 = empty_like(dLLH_ns0)
    #sig_ns0[dLLH_sig_ns0 < 0] = sqrt(-2.*dLLH_sig_ns0[dLLH_sig_ns0 < 0])
    #sig_ns0[dLLH_sig_ns0 > 0] = -sqrt(2.*dLLH_sig_ns0[dLLH_sig_ns0 > 0])
    sig_ns0 = empty_like(dLLH_ns0)
    sig_ns0[dLLH_sig_ns0 < 0] = norm.isf(chi2(1).sf(-dLLH_ns0[dLLH_sig_ns0 < 0]*2.)*0.5 )
    sig_ns0[dLLH_sig_ns0 > 0] = -norm.isf(chi2(1).sf(-dLLH_ns0[dLLH_sig_ns0 > 0]*2.)*0.5 )
    
    dLLHmin = dLLH_ns0 + logLBackground - (logL).min()
    sig_min = norm.isf(chi2(2).sf(dLLH_min*2.)*0.5 )
    
    
    #%%%
    ################
    ### PLOTTING ###
    ################
    
    ###########
    ### LLH ###
    ###########
    ms = 5
    fig = plt.figure(figsize=(n,m))
    plt.pcolormesh(interp, interp, logL, cmap="PuBu", vmin=(logL).min(), vmax=(logL).max())
    plt.plot(xlogL,ylogL, "+", label=r"$\log(\mathcal{L})_\mathrm{max}=%.1g$"%logL.max(), markersize=ms, color="k")
    plt.xlabel(r"azimuth $(\phi_\mathrm{event} - \phi_\mathrm{moon})\cdot \mathrm{sin}(\theta_\mathrm{event})$")
    plt.ylabel(r"zenith  $\theta_\mathrm{event} - \theta_\mathrm{moon}$")
    plt.gca().yaxis.set_major_formatter(ticker.FormatStrFormatter(u'%d°'))
    plt.gca().xaxis.set_major_formatter(ticker.FormatStrFormatter(u'%d°'))
    cbar = plt.colorbar()
    cbar.set_label(r"$\mathrm{log}(\mathcal{L})$")
    #plt.legend()
    plt.tight_layout()
    plt.grid(linestyle=(0,(5, 10)),  linewidth=0.5)
    fig.savefig("%s/%s_LLH_interpolated.png"%(pathSave, fileName), dpi=256)
    #%%%
    #################
    ### Delta LLH ###
    #################
    
    # 1. LLH - LLH(nSignal=0)
    fig = plt.figure(figsize=(n,m))
    plt.pcolormesh(interp, interp, dLLH_ns0, cmap="PuBu")
    plt.plot(xlogL,ylogL, "+", label="minimum", markersize=ms, color="k")
    plt.xlabel(r"azimuth $(\phi_\mathrm{event} - \phi_\mathrm{moon})\cdot \mathrm{sin}(\theta_\mathrm{event})$")
    plt.ylabel(r"zenith  $\theta_\mathrm{event} - \theta_\mathrm{moon}$")
    plt.gca().yaxis.set_major_formatter(ticker.FormatStrFormatter(u'%d°'))
    plt.gca().xaxis.set_major_formatter(ticker.FormatStrFormatter(u'%d°'))
    cbar = plt.colorbar()
    #cbar.set_label(r"$\Delta\mathrm{log}(\mathcal{L})$")
    cbar.set_label(r"$-\mathrm{log}(\mathcal{L})+\mathrm{log}(\mathcal{L}(n_\mathrm{s}=0))$")
    #plt.legend()
    plt.text(-2.8, 2.5, r"$\Delta\log(\mathcal{L})_\mathrm{max}$", color="w")
    plt.text(-2.8, 2.0, r"$\Delta\log(\mathcal{L})_\mathrm{min}$", color="w")
    plt.text(-.9, 2.5, r"$=$", color="w")
    plt.text(-.9, 2.0, r"$=$", color="w")
    plt.text(1.3, 2.48, (r"$%.0e"%dLLH_ns0.max()).replace('e', r'\cdot 10^{') + r'}$', color="w",  ha="right")
    plt.text(1.3, 2.0, r"$%.1f$"%dLLH_ns0.min(), color="w",  ha="right")
    plt.tight_layout()
    plt.grid(linestyle=(0,(5, 10)),  linewidth=0.5)
    fig.savefig("%s/%s_dLLH1_interpolated.png"%(pathSave, fileName), dpi=256)
    
    # 2. LLH - min(LLH)
    fig = plt.figure(figsize=(n,m))
    plt.pcolormesh(interp, interp, dLLH_min, cmap="PuBu")
    plt.plot(xlogL,ylogL, "+", label="minimum", markersize=ms, color="k")
    plt.xlabel(r"azimuth $(\phi_\mathrm{event} - \phi_\mathrm{moon})\cdot \mathrm{sin}(\theta_\mathrm{event})$")
    plt.ylabel(r"zenith  $\theta_\mathrm{event} - \theta_\mathrm{moon}$")
    plt.gca().yaxis.set_major_formatter(ticker.FormatStrFormatter(u'%d°'))
    plt.gca().xaxis.set_major_formatter(ticker.FormatStrFormatter(u'%d°'))
    cbar = plt.colorbar()
    cbar.set_label(r"$-\mathrm{log}(\mathcal{L})+\mathrm{log}(\mathcal{L}(n_\mathrm{s, min}))$")
    plt.text(-2.8, 2.5, r"$\Delta\log(\mathcal{L})_\mathrm{max}=%.1f$"%dLLH_min.max(), color="w")
    #plt.legend()
    plt.tight_layout()
    plt.grid(linestyle=(0,(5, 10)),  linewidth=0.5)
    fig.savefig("%s/%s_dLLH2_interpolated.png"%(pathSave, fileName), dpi=256)
    #%%%
    
    ################
    ### nSignal ####
    ################
    
    fig = plt.figure(figsize=(n,m))
    plt.pcolormesh(interp,interp, nSignal, cmap="seismic_r", vmin=-int(maximum)-1, vmax=int(maximum)+1)
    plt.plot(xnSignal,ynSignal, "+", label="minimum", markersize=ms, color="w")
    plt.xlabel(r"azimuth $(\phi_\mathrm{event} - \phi_\mathrm{moon})\cdot \mathrm{sin}(\theta_\mathrm{event})$")
    plt.ylabel(r"zenith  $\theta_\mathrm{event} - \theta_\mathrm{moon}$")
    plt.gca().yaxis.set_major_formatter(ticker.FormatStrFormatter(u'%d°'))
    plt.gca().xaxis.set_major_formatter(ticker.FormatStrFormatter(u'%d°'))
    cbar = plt.colorbar()
    cbar.set_label(r"$n_\mathrm{signal}$")
    plt.text(-2.8, 2.5, r"$n_\mathrm{signal,min}$")
    plt.text(-0.9, 2.0, r"$=$")
    plt.text(1.3, 2.5, r"$%.1f$"%nSignal.min(),  ha="right")
    plt.text(-2.8, 2.0, r"$n_\mathrm{signal,max}$")
    plt.text(-0.9, 2.5, r"$=$")
    plt.text(1.3, 2.0, r"$%.1f$"%nSignal.max(), ha="right")
    #plt.legend(framealpha=1)
    plt.grid(linestyle=(0,(5, 10)),  linewidth=0.5)
    plt.tight_layout()
    fig.savefig("%s/%s_nSignal_interpolated.png"%(pathSave, fileName), dpi=256)
    print(nSignal.min(), nSignal.max())
    #%%%
    
    ####################
    ### significance ###
    ####################
    
    
    fig = plt.figure(figsize=(n,m))
    cmap = plt.cm.get_cmap("seismic", 2*int(sig_ns0.max())+2)
    plt.pcolormesh(interp, interp, sig_ns0, cmap=cmap, vmin=-int(sig_ns0.max())-1, vmax=int(sig_ns0.max())+1)
    plt.grid(linestyle=(0,(5, 10)),  linewidth=0.5)
    plt.xlabel(r"azimuth $(\phi_\mathrm{event} - \phi_\mathrm{moon})\cdot \mathrm{sin}(\theta_\mathrm{event})$")
    plt.ylabel(r"zenith  $\theta_\mathrm{event} - \theta_\mathrm{moon}$")
    plt.gca().yaxis.set_major_formatter(ticker.FormatStrFormatter(u'%d°'))
    plt.gca().xaxis.set_major_formatter(ticker.FormatStrFormatter(u'%d°'))
    cbar = plt.colorbar()
    cbar.set_label(r"significance $\sigma(n_\mathrm{s}=0)$")
    plt.tight_layout()
    fig.savefig("%s/%s_sig_ns0_interpolated.png"%(pathSave, fileName), dpi=256)
    
    #%%%
    
    # significance 2: sigma(LLH - min(LLH))
    
    mask = sig_min < 1
    fig = plt.figure(figsize=(n,m))
    color = array(["w", "silver", "dimgrey", "k", "yellow", "gold", "darkorange", "red", "firebrick", "lime","mediumseagreen",
             "darkgreen", "lightcyan", "cyan", "dodgerblue", "navy", "lavenderblush", "pink", "deeppink"])[::-1]
    color = array(["w", "silver", "dimgrey", "k", "yellow", "gold", "darkorange", "red", "firebrick", "lime","mediumseagreen",
             "darkgreen", "lightcyan", "cyan", "dodgerblue", "navy", "lavenderblush", "pink", "deeppink", "silver", "dimgrey", "k", "yellow", "gold", "darkorange", "red", "firebrick", "lime","mediumseagreen",
             "darkgreen", "lightcyan", "cyan", "dodgerblue", "navy", "lavenderblush", "pink", "deeppink"])[::-1]
    cmap = colors.ListedColormap(color[len(color)-int(sig_min.max())-1:])
    plt.pcolormesh(interp, interp, sig_min, cmap=cmap, vmin=0, vmax=int(sig_min.max())+1)
    
    plt.grid(linestyle=(0,(5, 10)),  linewidth=0.5)
    plt.xlabel(r"azimuth $(\phi_\mathrm{event} - \phi_\mathrm{moon})\cdot \mathrm{sin}(\theta_\mathrm{event})$")
    plt.ylabel(r"zenith  $\theta_\mathrm{event} - \theta_\mathrm{moon}$")
    plt.gca().yaxis.set_major_formatter(ticker.FormatStrFormatter(u'%d°'))
    plt.gca().xaxis.set_major_formatter(ticker.FormatStrFormatter(u'%d°'))
    cbar = plt.colorbar()
    cbar.set_label(r"significance $\sigma(n_\mathrm{s,min}, n_\mathrm{dof}=2)$")
    #cbar.set_label(r"significance $\sigma(\Delta\mathrm{log}(\mathcal{L})(n_\mathrm{s,min}))$")
    plt.tight_layout()
    fig.savefig("%s/%s_sig_mindof2_interpolated.png"%(pathSave, fileName), dpi=256)
    
    #%%%
    # CONSITENCY
    ms = 4
    rMoon = 0.26
    
    manual_locations = [(-1, 1.4),(-0.1, -0.2)]
    
    pos = where(sig_ns0 == sig_ns0.max())
    yP = interp[pos[0][0]]
    xP = interp[pos[1][0]]
    
    fig, ax = plt.subplots(figsize=(n,m))
    ax.imshow(img, extent=[-rMoon,rMoon,-rMoon,rMoon])
    cs = ax.contour(interp, interp, sig_ns0, levels=[int(sig_ns0.max())-1,int(sig_ns0.max())], colors=["fireBrick","fireBrick"], linestyles=["dotted", "solid"])#colors=["dodgerblue","navy"])#
    #ax.clabel(cs, inline=0.001, fontsize=10, inline_spacing=-0.1, manual=manual_locations, fmt=r"$%d\sigma$")
    plt.grid(linestyle=(0,(5, 10)),  linewidth=0.5)
    plt.plot(xP, yP, "+", color="fireBrick")
    #plt.plot(1,1,color="navy",label="SplineMPE")
    #plt.plot(1,1,color="dodgerblue",label="MPE")
    plt.xlim(-1,1)
    plt.ylim(-1,1)
    plt.text(-0.9, 0.85, r"significance$=%.2f\sigma$"%(sig_ns0.max()))
    plt.xlabel(r"azimuth $(\phi_\mathrm{event} - \phi_\mathrm{moon})\cdot \mathrm{sin}(\theta_\mathrm{event})$")
    plt.ylabel(r"zenith  $\theta_\mathrm{event} - \theta_\mathrm{moon}$")
    plt.gca().yaxis.set_major_formatter(ticker.FormatStrFormatter(u'%.1f°'))
    plt.gca().xaxis.set_major_formatter(ticker.FormatStrFormatter(u'%.1f°'))
    plt.tight_layout()
    plt.plot(-10,-10, color="fireBrick", linestyle="dotted", label=r"$%d\sigma$"%(int(sig_ns0.max())-1))
    plt.plot(-10,-10, color="firebrick", linestyle="solid", label=r"$%d\sigma$"%(int(sig_ns0.max())))
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    fig.savefig("%s/%s_contour_ns0_interpolated.png"%(pathSave, fileName), dpi=256, bbox_inches='tight')
    #%%%
    # POINTING
    manual_locations = [(-1, 1.4),(-0.1, -0.2)]
    
    
    pos = where(sig_min == sig_min.min())
    yP = interp[pos[0][0]]
    xP = interp[pos[1][0]]
    err = where(sig_min<1)
    errAzMax = interp[err[1].max()] -xP
    errAzMin = xP - interp[err[1].min()]
    errZeMax = interp[err[0].max()] -yP
    errZeMin = yP - interp[err[0].min()]
    
    
    fig, ax = plt.subplots(figsize=(n,m))
    ax.imshow(img, extent=[-rMoon,rMoon,-rMoon,rMoon])
    cs = ax.contour(interp,interp,sig_min, levels=[1,3], colors=["fireBrick"], linestyles=["solid", "dotted"])#colors=["dodgerblue","navy"])
    plt.grid(linestyle=(0,(5, 10)),  linewidth=0.5)
    plt.plot(1,1,color="fireBrick", linestyle="solid", label=r"$1\sigma$")
    plt.plot(1,1,color="fireBrick", linestyle="dotted", label=r"$3\sigma$")
    plt.plot(xP, yP, "+", color="fireBrick")
    plt.xlim(-1,1)
    plt.ylim(-1,1)
    plt.text(-0.9, 0.86, r"Best Fit:")
    plt.text(-0.9, 0.71, r"$\Delta\phi\sin(\theta)$")
    plt.text(-0.9, 0.55, r"$\Delta\theta$")
    plt.text(-0.4, 0.71, r"$=$")
    plt.text(-0.4, 0.55, r"$=$")
    if xP >= 0:
        plt.text(0.8, 0.71, r"$(%.3f\pm^{%.3f}_{%.3f})^\circ$"%(abs(xP),errAzMax,errAzMin),  ha="right")
    else: 
        plt.text(0.8, 0.71, r"$-(%.3f\pm^{%.3f}_{%.3f})^\circ$"%(abs(xP),errAzMax,errAzMin),  ha="right")
    if yP >= 0:
        plt.text(0.8, 0.55, r"$(%.3f\pm^{%.3f}_{%.3f})^\circ$"%(abs(yP),errZeMax,errZeMin),  ha="right")
    else:
        plt.text(0.8, 0.55, r"$-(%.3f\pm^{%.3f}_{%.3f})^\circ$"%(abs(yP),errZeMax,errZeMin),  ha="right")
    plt.xlabel(r"azimuth $(\phi_\mathrm{event} - \phi_\mathrm{moon})\cdot \mathrm{sin}(\theta_\mathrm{event})$")
    plt.ylabel(r"zenith  $\theta_\mathrm{event} - \theta_\mathrm{moon}$")
    plt.gca().yaxis.set_major_formatter(ticker.FormatStrFormatter(u'%.1f°'))
    plt.gca().xaxis.set_major_formatter(ticker.FormatStrFormatter(u'%.1f°'))
    plt.tight_layout()
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    fig.savefig("%s/%s_contour_nsmin_interpolated.png"%(pathSave, fileName), dpi=256, bbox_inches='tight')
     #%%%
    
    print(sig_ns0.max(), interp[pos[1]], interp[pos[0]])



#%%%
dic = {"sig_min":sig_min, "sig_ns0":sig_ns0}
pickle.dump(dic, open("%s/%s_interp_eval.pickle"%(pathData,fileName), "wb"))


#low = [0, 36, 72, 108, 144, 180, 216, 252, 288, 324]
#up = [30, 66, 102,  138, 174, 210, 246, 282, 318, 354]


# fig size
n = 4       # x- size
m = 3       # y-size


"""
#####################
### 1D Histogramm ###
#####################

fig = plt.figure(figsize=(6,4))
plt.axhline(0, color="k", linewidth=1, zorder=-1)
plt.errorbar((b[1:]+b[:-1])/2, (h - hBGD/34)/(hBGD/34), xerr =(b[1:]-b[:-1])/2, yerr=sqrt(h)/hBGD*34 , fmt=".", color="r")
plt.xlabel(r"distance between event and moon $\Delta\Psi=\sqrt{\Delta\phi^2 + \Delta\theta^2}$")
plt.ylabel(r"$\frac{N_\mathrm{source}}{N_\mathrm{background}}- 1$", fontsize=14)
plt.grid(linestyle=(0,(5, 10)),  linewidth=0.5, zorder=-10)
plt.gca().xaxis.set_major_formatter(ticker.FormatStrFormatter(u'%d°'))
plt.tight_layout()
fig.savefig("%s/%s_1dhistdeficit.png"%(pathSave, fileName), dpi=256)
"""


