import os, numpy, math, random#, csv, fnmatch, , math, seaborn as sns #glob #pandas, re, glob, sys, 
import matplotlib.pyplot as plt
# from scipy.stats import norm, stats
from scipy.stats import mannwhitneyu, wilcoxon, shapiro, normaltest, anderson, t
import matplotlib.patches as mpatches
# from timeit import default_timer as timer
# import matplotlib.axes as Axes
# from time import time
# from mpl_toolkits.mplot3d import Axes3D
from Blinds_Core_Functions_0030 import findtitles, histogram, read_all_data, append_value, samplecomparison_histogram, bootstrap_confidence_intervals, dictionarystats, uncertainty_plotter, mean_and_uncertainty_plot, read_bendstrain #readstats_readlastframe, remove_chosen_outliers
# import Core_Functions_0001 #Imports functions and the modules. 

from Dataset_file_locations import get_dataset_info 
from statistics import NormalDist

# first_start=timer()
# start = timer()

tablefilename = 'fixed_zero_stats_04.txt'   #The .txt file name to read roughness data from
txtfilepaths=[]
dataset_names=[]
root_1='F:/3.SiC Blinds/1.ILC/1.Section 1/'
root_2='F:/3.SiC Blinds/3.ILC Repeat/1.Section 1/'

##Get txtfilepaths for the different samples, the graph saving file root location and the dataset name.
"""dataset_number=3 corresponds to the "1st Dataset" 
   dataset_number=4 corresponds to the "2nd Dataset"
   For purposes of naming, the 1st and 2nd datasets are 
   called 3rd and 4th respectively in the code due to obsolete 
   datasets being previously used. 
   """

txtfilepathways, graphsave_root, dataset_name, roughness_par_names, roughness_par_labels, roughness_par_titles = get_dataset_info(Method_number=1, dataset_number=3)
# print(dataset_name)
dataset_names.append(dataset_name)
txtfilepaths.append(txtfilepathways)
roots = [root_2]

txtfilepathways, graphsave_root, dataset_name, roughness_par_names, roughness_par_labels, roughness_par_titles = get_dataset_info(Method_number=1, dataset_number=4)
# print(dataset_name)
dataset_names.append(dataset_name)
txtfilepaths.append(txtfilepathways)
roots = [root_2, root_2]

graphsave_root='F:/3.SiC Blinds/1.ILC/Graphs/'
sampletitles=findtitles(roots[0])

##Read data into dictionaries. Outlier data is not exported and as such an underscore is used to ignore that output.
all_data, all_clean_data, _ = read_all_data(roots, txtfilepaths, dataset_names, 'fixed_zero_stats_04.txt')
bendstrain_data=read_bendstrain()
bendstrain_stats=dictionarystats(bendstrain_data)

del dataset_name, txtfilepathways
graphsave_loc = graphsave_root+'Single Histograms/'
# os.chdir(graphsave_loc)    

# def bendstrain_graphs():
#     sampletitles=findtitles(roots[0])
#     tagline='A03 Anomaly Removed'
#     os.chdir('F:/3.SiC Blinds/')
#     mean_and_uncertainty_plot(bendstrain_data, 0, 0, sampletitles, 'Bendstrain Results', tagline, savefile=True)

def IQRoutliermethod(sample_data, column_number, method_number, anomaliestxtfilename, write_to_file=False):
    """Basic level function that identifies the frames that would be selected 
    by an IQR Filter of the parameter associated with the column number in 
    the .txt file.

    e.g. column number 2=Sa, 3=Sq, 4=Ssk, 5=Sku and 6=Sz"""
    
    # global outlier, Sa_outliers, Distance_outliers, chosen_outliers, outliernumber
    outliernumber=[]        
    Outliervalue=[]         #Creating the list of the frame's roughness values that would be identified as anomalous
    Distance_outliers=[]    #Creating the list of frames that will be identified as anomalous
    roughness_data = sample_data[roughness_par_names[column_number]]    #Finds the roughness parameter data
    
    if column_number==1:
        b='Average Height'
        xscale=1e9
        title='Average Height(nm)'
    
    elif column_number==2 and method_number==1:
        b='Sa'
        xscale=1e9
        title='Sa(nm)'
    elif column_number==2 and method_number==2:
        b='Ra'
        xscale=1e9
        
    elif column_number==3 and method_number==1:
        b='Sq'
        xscale=1e9
        title='Sq(nm)'    
    elif column_number==3 and method_number==2:
        b='Rq'
        xscale=1e9
        title='Rq(nm)'
   
    elif column_number==4:
        b='Skew'
        xscale=1
        title='Skew'
  
    elif column_number==5:
        b='Kurtosis'
        xscale=1
        title='Kurtosis'
        # High_cutoff=1
    
    elif column_number==6 and method_number==1:
        b='Sz'
        xscale=1e9
        title='Sz(nm)'
    
    elif column_number==6 and method_number==2:
        b='Rz'
        xscale=1e9
        title='Rz(nm)'
    else:
        pass
    
    roughness_data=[i * xscale for i in roughness_data]
    
    Q1=numpy.quantile(roughness_data, 0.25)
    Q3=numpy.quantile(roughness_data, 0.75)
    ##Define the Interquartile Range
    IQR = Q3-Q1
    
    ##Define the lower and upper cut-offs for the IQR filter, if undefined previously, (defined for skew and kurtosis)
    try:
        High_cutoff #Variable is defined already
    except NameError:
        High_cutoff=Q3+1.5*IQR      #Variable is undefined, and so we define it. 
    try:
        Low_cutoff #Variable is defined already
    except NameError:
        Low_cutoff=Q1-1.5*IQR      #Variable is undefined, and so we define it. 
    if write_to_file==True:
        with open(anomaliestxtfilename, "a") as file_object:
            file_object.write('\n  Lower 1.5IQR <%.3g'%(Low_cutoff)+' and Upper 1.5IQR >%.3g'%(High_cutoff)+
                              '\n     {:^10s}{:^3s}{:>10s}'.format('No.','Dist(um)',title))
            file_object.close()
    # print('  Lower 1.5IQR <%.3g'%(Low_cutoff)+' and Upper 1.5IQR >%.3g'%(High_cutoff))
    # print('     {:^10s}{:^3s}{:>10s}'.format('No.','Dist(um)',title))
    for i, row in enumerate(roughness_data):
        if row>High_cutoff or row<Low_cutoff:
            n=i+1
            outliernumber.append(n)                               #Row of the outlier
            Distance_outliers.append(sample_data["Distance (um)"][i])
            Outliervalue.append(int(row))
            if write_to_file==True:
                with open(anomaliestxtfilename, "a") as file_object:
                    file_object.write('\n     {:^10d}{:>5d}{:>12.4g}'.format(n, sample_data["Distance (um)"][i], row))
                    file_object.close()
            # print('     {:^10d}{:>5d}{:>12.4g}'.format(n, sample_data["Distance (um)"][i], row))
        else:
            pass
    with open(anomaliestxtfilename, "a") as file_object:
        file_object.write("\n       Number of identified outlier frames = %i" %len(outliernumber)+"\n")
        file_object.close()
    # print ("       Number of identified outlier frames = %i" %len(outliernumber)+"\n")
    return Distance_outliers

def all_samples_IQRMethod(dataset_number):
    """Higher Level function of the "IQRoutliermethod" function that filters" 
    dataset_number=1 corresponds to the "section 1" dataset where the data was collected with 1 tip throughout the whole scan
    dataset_number=2 corresponds to the "2nd Good Dataset" where the data was sometimes repeated and sometimes done on a different section, but with tip replacement
    dataset_number=3 corresponds to the "3rd Good Dataset" 
    dataset_number=4 corresponds to the "4th Dataset" - identical to 2nd Good Dataset, but A08 is taken from section 4.""" # and plots all single graphs using IQR Method
    
    dataset_name=dataset_names[dataset_number-1]
    # for x, dataset in enumerate(dataset_names):   #Selecting which dataset 
    dash = '-' * 30
    # print('\n'+str(dataset_name))
    sampletitles=findtitles(roots[dataset_number-1])
    
    anomaliestxtfilename="Anomalous Frames - "+dataset_name+".txt"
    f = open(anomaliestxtfilename, "w+") 
    f.write(str(dataset_name) +"\n")
    f.close()
    
    for countr, samplename in enumerate(sampletitles):
        # print('  '+ str(samplename))
        with open(anomaliestxtfilename, "a") as file_object:
            file_object.write('\n'+str(samplename)+' ('+ txtfilepaths[dataset_number-1][countr]+')\n')
            file_object.close()
        sample_data= all_data[dataset_name][samplename]
        Distance_outliers=[]
        for column_number, rough_par in enumerate(roughness_par_names):
                if column_number==0:
                    pass
                else:
                    parameter_outliers=IQRoutliermethod(sample_data, column_number, 1, anomaliestxtfilename, write_to_file=True)
                    for outlier in parameter_outliers:
                        Distance_outliers.append(outlier)
        ##Remove Duplicates ???
        Distance_outliers=list(set(Distance_outliers))
        ##Order the outliers
        Distance_outliers.sort()
        
        with open(anomaliestxtfilename, "a") as file_object:
            file_object.write('\n'+str(Distance_outliers)+'\n'+dash+'\n')
            file_object.close()
        
        # print("\n"+str(Distance_outliers))
        
        # print('\n'+dash+"\n\n")
        # time.sleep(5)
        
        #Save dimensions and profile positions to a .txt file
        
        
    # remove_chosen_outliers(0)
    # IQRoutliermethod(6)           #If a combination of IQR filters is to be examined...
    # remove_chosen_outliers()      #As above.
    # for j in [Avg, Ra, Rms, Skew, Kurtosis, Zmax]:
    #     n_sites=int(len(j))       #int(round(len(j)/2))
    #     histogram(j, sampletitles[x])
        # boxplot(j, sampletitles[x])
   


def mean_convergence(rootpath_of_dataset, include_offset=False, do_iterations=False, correct_order=False, savefile=True):
    """Old function used for early development. mean_convergence_2 is the new function (see later)
    Plots the cumulative mean of data in accordance to the current order showing it converge onto the sample mean.
        'include_offset' - if TRUE...includes the 2nd offset dataset to increase the number of measurements 
                           if FALSE...uses just the first dataset in the graph plotting
        'do_iterations' - if TRUE...scrambles the data order 20,000 times and outputs the average maximum deviation of the 
                                    cumulative mean away from the sample mean after 200 or 400 measurements for the 1 or 
                                    combined datasets respectively and prints out the results
                          if FALSE... ignores this
        'correct_order' - if TRUE... only plots the convergence graph for the parameters in the correct order"""
    
    sampletitles=findtitles(rootpath_of_dataset)
    # global xdata, ydata, cum_sum
    ###For each sample
    for sampledataset_number, samplename in enumerate(sampletitles):
        ###Iterate through roughness parameters for plotting graphs/estimating errors
        for column_number, roughness_parameter in enumerate(roughness_par_names):
            cum_sum=0
            cum_stdev=0
            txt_x_loc=200       #The maximum deviation calculated 
            # os.chdir(core_path+'Cumulative Mean Plots/Fibre 1 Dataset')
            roughness_data={}
            dataset="Section 1"
            if column_number==0:    #Skips the Distance column
                pass
            else:                            
                if column_number==1:
                    b='Average Height'
                    xscale=1e9
                    axislabel='Average Height (nm)'
                    units='nm'
                    countlim=5          #The number of graph iterations to plot for each sample's roughness parameter
                    allowed_err=10
                    valuescale=15
                elif column_number==2:
                    b='Sa'
                    xscale=1e9
                    axislabel='Sa Mean (nm)'
                    units='nm'
                    countlim=5
                    allowed_err=5
                    valuescale=5
                elif column_number==3:
                    b='Sq'
                    xscale=1e9
                    axislabel='Sq Mean (nm)'
                    units='nm'
                    countlim=5
                    allowed_err=5
                    valuescale=5 
                elif column_number==4:
                    b='Skew'
                    xscale=1
                    axislabel='Skew Mean'
                    countlim=0
                    units=''
                elif column_number==5:
                    b='Kurtosis'
                    xscale=1
                    axislabel='Kurtosis Mean'
                    countlim=0
                    units=''
                elif column_number==6:
                    b='Sz'
                    xscale=1e9
                    axislabel='Sz Mean (nm)'
                    units='nm'
                    countlim=5
                    allowed_err=10
                    valuescale=10
                else:
                    pass
                data=all_clean_data[dataset][samplename][roughness_parameter]
                n_sites_original=int(len(data))                   
                
                data=numpy.array(data)*xscale
                n_sites=int(len(data))
                mu = numpy.mean(data)
                sigma = numpy.std(data)
                counter=0
                max_err=[]
                max_err_pos=[]
                
                if correct_order==True: 
                ###Plots the mean convergence graph in the order that the data was collected without any random shuffling.
                    xdata=[]
                    yreaddata=[]
                    ydata=[]
                    ystdevdata=[]
                    ystderr=[]
                    exp_unc_upper=[]
                    exp_unc_lower=[]
                    unc_upper=[]
                    unc_lower=[]
                    ultraExp_unc_upper=[]
                    ultraExp_unc_lower=[]
                    n=[]
                    selecteddata=data
                    for n, row in enumerate(selecteddata,1):
                        if n>1:
                            cum_sum=cum_sum+row
                            cum_mean=cum_sum/n
                            ydata.append(cum_mean)
                        elif n==1:
                            cum_sum=row
                        xdata.append(n)
                        yreaddata.append(row)
                        cum_stdev=numpy.std(yreaddata)
                        sqrtn=math.sqrt(n)
                        cum_sterr=cum_stdev/sqrtn
                        ystdevdata.append(cum_stdev)
                        ystderr.append(cum_sterr)
                    cum_sum=0
                    suptitle=u'n= %.3g, \u03bc= %.3g, \u03C3= %.3g'%(n_sites, mu, sigma)
                    
                    ax1=plt.subplot(1,1,1)
                    plt.plot(xdata[1:],ydata[:], label='Cumulative Mean')
                    
                    
                    plt.xlim(0,(int(n_sites/50))*50+50)
                    
                    ##Set y scale limits
                    # ymin, ymax = plt.ylim()
                    # yrange=ymax-ymin
                    # # if ((mu-ymin)/mu)>((ymax-mu)/mu):     ##Sets Y limit to accommdate data
                    # #     plt.ylim(ymin,(mu+mu-ymin))
                    # # elif mu-(ymax-mu)<0:
                    # #     plt.ylim(0,ymax)
                    # # else:
                    # #     plt.ylim(mu-(ymax-mu),ymax)
                    
                    # ## OR
                    # plt.ylim(mu-allowed_err,mu+allowed_err)   ##Sets y limit to just the mean+error
                    # ##OR
                    # # plt.ylim(mu-allowed_err*valuescale,mu+allowed_err*valuescale) ##Sets the ylimit to mean+error + extra room for data.
                    # ymin, ymax = plt.ylim()
                    
                    plt.axhline(y=mu,c='green', alpha=0.4, label='Sample Mean')
                    
                    plt.xlabel('Number of Measurements, n')
                    plt.ylabel(axislabel)
                    plt.suptitle(samplename+', '+b+' - correct order', fontsize=14)
                    plt.title(suptitle)
                    
                    for n, row in enumerate(ydata, 1):
                        unc_upper.append(row+ystderr[n])    #upper band of 1u uncertainty band representing ~68% data
                        unc_lower.append(row-ystderr[n])    #lower band of 1u uncertainty band representing ~68% data
                        exp_unc_upper.append(row+(2*ystderr[n]))    #upper band of 2u uncertainty band representing 95% data
                        exp_unc_lower.append(row-(2*ystderr[n]))    #lower band of 2u uncertainty band representing 95% data
                        ultraExp_unc_upper.append(row+(3*ystderr[n]))   #upper band of 3u uncertainty band representing 99.7% data
                        ultraExp_unc_lower.append(row-(3*ystderr[n]))   #lower band of 3u uncertainty band representing 99.7% data
                    alpha_1=0.08
                    alpha_2=0.07
                    alpha_3=0.06
                    plt.fill_between(xdata[1:], unc_upper, unc_lower, alpha=alpha_1, color='b')
                    plt.fill_between(xdata[1:], exp_unc_upper, exp_unc_lower, alpha=alpha_2, color='b')
                    plt.fill_between(xdata[1:], ultraExp_unc_upper, ultraExp_unc_lower, alpha=alpha_3, color='b')
                    
                    ##Creating the legend
                    # where some data has already been plotted to ax
                    handles, labels = ax1.get_legend_handles_labels()
                    # manually define new patches for the different bands 
                    u1_band = mpatches.Patch(color='b', alpha=alpha_1+alpha_2+alpha_3, label='68% Coverage Probability')
                    u2_band = mpatches.Patch(color='b', alpha=alpha_2+alpha_3, label='95% Coverage Probability')
                    u3_band = mpatches.Patch(color='b', alpha=alpha_3, label='99.7% Coverage Probability')
                    handles.append(u1_band) 
                    handles.append(u2_band) 
                    handles.append(u3_band) 
                    plt.legend(handles=handles, prop={"size":8}, loc=1)
                    
                    plt.show()
                    # plt.pause(1)
                    if savefile==True:
                        plt.savefig('Correct Order Convergence Graph - '+samplename +' - '+ b +'.png') 
                    plt.close()
                
                else:
                ##Continues with convergence plot scrambling
                    iterations=20000
                    while counter<countlim:    #countlim determines number of scrambled graphs to plot, before everything done in background. 
                        xdata=[]
                        yreaddata=[]
                        ydata=[]
                        ystdevdata=[]
                        ystderr=[]
                        exp_unc_upper=[]
                        exp_unc_lower=[]
                        unc_upper=[]
                        unc_lower=[]
                        ultraExp_unc_upper=[]
                        ultraExp_unc_lower=[]
                        n=[]
                        selecteddata=data
                        random.shuffle(selecteddata)
                        for n, row in enumerate(selecteddata,1):
                            if n>1:
                                cum_sum=cum_sum+row
                                cum_mean=cum_sum/n
                                ydata.append(cum_mean)
                            elif n==1:
                                cum_sum=row
                            xdata.append(n)
                            yreaddata.append(row)
                            cum_stdev=numpy.std(yreaddata)
                            sqrtn=math.sqrt(n)
                            cum_sterr=cum_stdev/sqrtn
                            ystdevdata.append(cum_stdev)
                            ystderr.append(cum_sterr)
                        cum_sum=0
                        suptitle=u'n= %.3g, \u03bc= %.3g, \u03C3= %.3g'%(n_sites, mu, sigma)
                        
                        ax1=plt.subplot(1,1,1)
                        plt.plot(xdata[1:],ydata[:], label='Cumulative Mean')
                        
                        #plt.plot(xdata[1:],selecteddata[1:], 'r', alpha=0.1, label='Data')
                        
                        plt.xlim(0,(int(n_sites/50))*50+50)
                        
                        ##Set y scale limits
                        # ymin, ymax = plt.ylim()
                        # yrange=ymax-ymin
                        # # if ((mu-ymin)/mu)>((ymax-mu)/mu):     ##Sets Y limit to accommdate data
                        # #     plt.ylim(ymin,(mu+mu-ymin))
                        # # elif mu-(ymax-mu)<0:
                        # #     plt.ylim(0,ymax)
                        # # else:
                        # #     plt.ylim(mu-(ymax-mu),ymax)
                        
                        # ## OR
                        # plt.ylim(mu-allowed_err,mu+allowed_err)   ##Sets y limit to just the mean+error
                        # ##OR
                        # # plt.ylim(mu-allowed_err*valuescale,mu+allowed_err*valuescale) ##Sets the ylimit to mean+error + extra room for data.
                        # ymin, ymax = plt.ylim()
                        
                        ###Calculates the maximum deviation of the cumulative mean line from the sample mean... 
                        if include_offset==False:           ###...after 200 measurements for 1 dataset
                            n_err=200   #The number of measurements before the max deviation calculation is done
                            diff_200=abs(ydata[200]-mu)
                            curr_max=diff_200
                            for x, i in enumerate(ydata[200:],200):
                                if abs(i-mu)>curr_max:
                                    curr_max= abs(i-mu)
                                    pos=x
                                else:
                                    pass
                            max_err.append(curr_max)
                            # max_err_pos.append(pos)

                        # plt.text(txt_x_loc, 0.05*yrange+ymin, u'Max \u0394=%.3g'%(curr_max)+units)
                        plt.axhline(y=mu,c='green', alpha=0.4, label='Sample Mean')
                        # plt.axhline(y=ydata[250],c='orange', alpha=0.3)
                        
                        plt.xlabel('Number of Measurements, n')
                        plt.ylabel(axislabel)
                        # plt.suptitle('Convergence of the Mean for '+sampletitles[sampledataset_number-1]+'\n'+b, fontsize=14)
                        plt.suptitle(samplename+', '+b+' - '+str(counter+1), fontsize=14)
                        plt.title(suptitle)
                        
                        for n, row in enumerate(ydata, 1):
                            unc_upper.append(row+ystderr[n])    #upper band of 1u uncertainty band representing ~68% data
                            unc_lower.append(row-ystderr[n])    #lower band of 1u uncertainty band representing ~68% data
                            exp_unc_upper.append(row+(2*ystderr[n]))    #upper band of 2u uncertainty band representing 95% data
                            exp_unc_lower.append(row-(2*ystderr[n]))    #lower band of 2u uncertainty band representing 95% data
                            ultraExp_unc_upper.append(row+(3*ystderr[n]))   #upper band of 3u uncertainty band representing 99.7% data
                            ultraExp_unc_lower.append(row-(3*ystderr[n]))   #lower band of 3u uncertainty band representing 99.7% data
                        alpha_1=0.08
                        alpha_2=0.07
                        alpha_3=0.06
                        plt.fill_between(xdata[1:], unc_upper, unc_lower, alpha=alpha_1, color='b')
                        plt.fill_between(xdata[1:], exp_unc_upper, exp_unc_lower, alpha=alpha_2, color='b')
                        plt.fill_between(xdata[1:], ultraExp_unc_upper, ultraExp_unc_lower, alpha=alpha_3, color='b')
                        
                        ##Creating the legend
                        # where some data has already been plotted to ax
                        handles, labels = ax1.get_legend_handles_labels()
                        # manually define new patches for the different bands 
                        u1_band = mpatches.Patch(color='b', alpha=alpha_1+alpha_2+alpha_3, label='68% Coverage Probability')
                        u2_band = mpatches.Patch(color='b', alpha=alpha_2+alpha_3, label='95% Coverage Probability')
                        u3_band = mpatches.Patch(color='b', alpha=alpha_3, label='99.7% Coverage Probability')
                        handles.append(u1_band) 
                        handles.append(u2_band) 
                        handles.append(u3_band) 
                        plt.legend(handles=handles, prop={"size":9}, loc=1)
                        
                        plt.show()
                        plt.pause(1)
                        # plt.savefig(sampletitles[sampledataset_number-1] +' - '+ b +' - '+str(counter+1)+'.png') 
                        # plt.savefig('Stdev - '+sampletitles[sampledataset_number-1] +' - '+ b +' - '+str(counter+1)+'.png')
                        if savefile==True:
                            plt.savefig('Mean+err-labelled - '+sampletitles[sampledataset_number-1] +' - '+ b +' - '+str(counter+1)+'.png') 
                        plt.close()
                        
                        # print(sampletitles[sampledataset_number-1]+' - '+ b)
                        # print('Mean error (above n=%i) after %i iterations = %.3f'%(n_err, iterations, mean_err)+units+' (%.3f%%)'%percent_mean_err)
                        # print('Max error (above n=%i) after %i iterations = %.3f'%(n_err, iterations, max_iter_err)+units+' (%.3f%%)'%percent_max_iter_err +'\n')
                        
                        # suptitle=u'n= %.3g'%(n_sites)
                        # plt.plot(xdata[1:],ystderr[1:], label='Standard Error')
                        # plt.title(suptitle)
                        # plt.ylim(ymin=0)
                        # plt.xlim(xmin=0)
                        # plt.suptitle('Standard Error - '+sampletitles[sampledataset_number-1]+', '+b+' - '+str(counter+1), fontsize=14)
                        # plt.xlabel('Number of Measurements, n')
                        # plt.grid(b=True, which='major', axis='y', linestyle='-', alpha=0.5)
                        # # plt.grid(b=True, which='minor', axis='y', linestyle='-', color='gray')
                        # plt.ylabel(axislabel)
                        
                        # plt.legend(prop={"size":9}, loc=1)
                    
                        # plt.show()
                        # plt.pause(1)
                        # plt.savefig('Sterr - '+sampletitles[sampledataset_number-1] +' - '+ b +' - '+str(counter+1)+'.png')
                        # plt.close()
                        # print(sampletitles[sampledataset_number-1]+' - '+ b+u'   Max \u0394=%.3g'%(curr_max)+units)
                        counter+=1
                    """Do iterations of 20,000 reorders and plot mean/max difference"""
                    if do_iterations==True:     #if continuing with iterations above countlim, do them. 
                        while counter>=countlim and counter>0 and counter<iterations:
                            n=[]
                            xdata=[]
                            ydata=[]
                            selecteddata=data
                            random.shuffle(selecteddata)
                            for n, row in enumerate(selecteddata,1):
                                if n>1:
                                    xdata.append(n)
                                    cum_sum=cum_sum+row
                                    cum_mean=cum_sum/n
                                    ydata.append(cum_mean) 
                                elif n==1:
                                    cum_sum=row
                            cum_sum=0
                            if include_offset==False:
                                diff_200=abs(ydata[200]-mu)
                                curr_max=diff_200
                                for x, i in enumerate(ydata[200:],200):
                                    if abs(i-mu)>curr_max:
                                        curr_max= abs(i-mu)
                                        pos=x
                                    else:
                                        pass
                                max_err.append(curr_max)
                                max_err_pos.append(pos)
                            if include_offset==True:
                                diff_400=abs(ydata[400]-mu)
                                curr_max=diff_400
                                for x, i in enumerate(ydata[400:],400):
                                    if abs(i-mu)>curr_max:
                                        curr_max= abs(i-mu)
                                        pos=x
                                    else:
                                        pass
                                max_err.append(curr_max)
                                max_err_pos.append(pos)
                            counter+=1
                        if counter==iterations:
                            mean_err=numpy.mean(max_err)
                            percent_mean_err=(mean_err/mu)*100
                            max_iter_err=max(max_err)
                            percent_max_iter_err=(max_iter_err/mu)*100
                            print(sampletitles[sampledataset_number-1]+' - '+ b)
                            print('Mean error (above n=%i) after %i iterations = %.3f'%(n_err, iterations, mean_err)+units+' (%.3f%%)'%percent_mean_err)
                            print('Max error (above n=%i) after %i iterations = %.3f'%(n_err, iterations, max_iter_err)+units+' (%.3f%%)'%percent_max_iter_err +'\n')
                            counter=0
                        else:
                            pass
    print ('Graphs Completed')  
 
def convergence_point_finder(scrambled_data, mu, final_lower, final_upper, convergence_pointlimit, last_convergence_points, counter, first_or_last_convergence='first'):
    n=[]
    ydata=[]
    old_n=-2
    yreaddata=[]
    convergence_counter=[]
    convergence_pointlist=[]
    cum_sum=0
    for n, row in enumerate(scrambled_data,1):
        if n>1:
            cum_sum=cum_sum+row
            cum_mean=cum_sum/n
            ydata.append(cum_mean)
        elif n==1:
            cum_sum=row
            cum_mean=cum_sum
        yreaddata.append(row)
        cum_stdev=numpy.std(yreaddata)
        sqrtn=math.sqrt(n)
        cum_sterr=cum_stdev/sqrtn
        
        #Define the multiplying factor for the error bands, according to standard deviations. i.e. 1=67% confidence, etc. The smaller the better. 
        mult_fac=1
        
        #Define the current upper and lower error band
        upper_errorband=cum_mean+mult_fac*cum_sterr
        lower_errorband=cum_mean-mult_fac*cum_sterr
        
        # #If the "true" roughness value is inbetween the error bands, save the value
        # if lower_errorband<mu<upper_errorband:
        #     if convergence_counter==[]:
        #         convergence_counter=[n]
        #         # print('first time', n)
        #     elif old_n==n-1 and len(convergence_counter)<=convergence_pointlimit:
        #         convergence_counter.append(n)
        #         # print('registered ',n)
        #     elif old_n==n-1 and len(convergence_counter)>convergence_pointlimit:
        #         pass
        #     else:    
        #         # print('new ',n)
        #         convergence_counter=[n]
        #     old_n=n    
        #     if len(convergence_counter)==convergence_pointlimit:
        #         # print('saved ',n,convergence_counter[0])
        #         convergence_pointlist.append(convergence_counter[0])
        
        #If the "true" roughness value's error bands are inbetween the current error bands, save the value
        if lower_errorband<final_lower<upper_errorband or lower_errorband<final_upper<upper_errorband:
            if convergence_counter==[]:
                convergence_counter=[n]
                # print('first time', n)
            elif old_n==n-1 and len(convergence_counter)<=convergence_pointlimit:
                convergence_counter.append(n)
                # print('registered ',n)
            elif old_n==n-1 and len(convergence_counter)>convergence_pointlimit:
                pass
            else:   
                # print('new ',n)
                convergence_counter=[n]
            old_n=n    
            if len(convergence_counter)==convergence_pointlimit:
                convergence_pointlist.append(convergence_counter[0])
                
    cum_sum=0
    if convergence_pointlist==[]:
        pass
    else:
        if first_or_last_convergence=='first':
            last_convergence_points.append(convergence_pointlist[0])
        elif first_or_last_convergence=='last':
            last_convergence_points.append(convergence_pointlist[-1])
        else:
            print('please select either "first" or "last" for first_or_last_convergence parameter')
            exit()
        convergence_pointlist=[]
        counter+=1
    # print('complete',counter)
    return last_convergence_points, counter
        
def confidence_of_overlap_finder(scrambled_data, best_mu, best_std_unc, convergence_point_to_test, list_of_confidence_overlaps, reverse_counter):
    """Takes the mean and standard uncertainty of the final measurement and compares it to the cum. mean and std unc. at the 'convergence point'.
    """
    
    #Calculate the mean and std uncertainty at the convergence point for this scrambled iteration of data.
    #Goes through the data until the 90th percentile has been reached.
    measured_mean=numpy.mean(scrambled_data[0:convergence_point_to_test])
    measured_stdev=numpy.std(scrambled_data[0:convergence_point_to_test])
    sqrtn=math.sqrt(convergence_point_to_test)
    measured_sterr=measured_stdev/sqrtn
    
    #Calculates the overlap region between large sample distribution and the current iteration distribution at the 90th percentile convergence point.
    current_overlap_amount=NormalDist(mu=best_mu, sigma=best_std_unc).overlap(NormalDist(mu=measured_mean, sigma=measured_sterr))
    
    list_of_confidence_overlaps.append(current_overlap_amount)
    reverse_counter+=1
    
    # print('complete ', reverse_counter)
    return list_of_confidence_overlaps, reverse_counter


def mean_convergence_2(rootpath_of_dataset, first_or_last_convergence='first', no_of_converged_points=20, plot_graph=False, do_iterations=False, iteration_no=1000, reverse_calc_accuracy=False, savefile=True):
    """Plots the cumulative mean of data in accordance to the current order showing it converge onto the sample mean.
        'do_iterations' - if TRUE...scrambles the data order 20,000 times and outputs the average maximum deviation of the 
                                    cumulative mean away from the sample mean after 200 or 400 measurements for the 1 or 
                                    combined datasets respectively and prints out the results
                          if FALSE... ignores this"""
    # global data
    sampletitles=findtitles(rootpath_of_dataset)
    iterations=iteration_no                             #Define the total number of iterations
    convergence_pointlimit=no_of_converged_points       #Define the number of consecutive cumulative mean roughnesses points need to have error bands 
                                                        #that envelope the "true" roughness before being identified as a convergence point. 
                                                        # e.g. 1=just the first point, 10 = first 10 points
    
    graphsave_loc= graphsave_root+'Test_2/'
    os.chdir(graphsave_loc)
    ###For each sample
    for sampledataset_number, samplename in enumerate(sampletitles):
        if do_iterations==False:
            print(samplename)
        counter=0
        reverse_counter=0
        ###Iterate through roughness parameters for plotting graphs/estimating errors
        for column_number, roughness_parameter in enumerate(roughness_par_names):
            cum_sum=0
            cum_stdev=0
            countlim=4              #The number of graph iterations to plot for each sample's roughness parameter
            if do_iterations==False:
                print('\n',roughness_parameter)
                
            if column_number==0 or column_number==4 or column_number==5:    #Skips the Distance, Kurtosis and Skew columns
                pass
            else:                            
                if column_number==1:
                    b='Average Height'
                    xscale=1e9
                    axislabel='Average Height (nm)'
                    units='nm'
                    allowed_err=10
                    valuescale=15
                elif column_number==2:
                    b='Sa'
                    xscale=1e9
                    axislabel='Sa Mean (nm)'
                    units='nm'
                    allowed_err=5
                    valuescale=5
                elif column_number==3:
                    b='Sq'
                    xscale=1e9
                    axislabel='Sq Mean (nm)'
                    units='nm'
                    allowed_err=5
                    valuescale=5 
                elif column_number==4:
                    b='Skew'
                    xscale=1
                    axislabel='Skew Mean'
                    countlim=0              #Ignore skew
                    units=''
                elif column_number==5:
                    b='Kurtosis'
                    xscale=1
                    axislabel='Kurtosis Mean'
                    countlim=0              #Ignore Kurtosis
                    units=''
                elif column_number==6:
                    b='Sz'
                    xscale=1e9
                    axislabel='Sz Mean (nm)'
                    units='nm'
                    allowed_err=10
                    valuescale=10
                else:
                    pass
                #Combine the data from the two datasets
                data=all_clean_data[dataset_names[2]][samplename][roughness_parameter] + all_clean_data[dataset_names[3]][samplename][roughness_parameter]
                
                
                """=========================Changed datasets to 2 and 3 in the line above ^^^ ===================================================="""
                
                n_sites_original=int(len(data)) 
                counter=0
                
                data=numpy.array(data)*xscale
                convergence_pointlist=[]     #List of all the first convergence points for each iteration
                last_convergence_points=[]   #List of last convergence points
                
                ##Get the mean, st dev and number of measurements for the dataset as well as the sample mean uncertainties.
                n_sites=int(len(data))
                mu = numpy.mean(data)
                sigma = numpy.std(data)
                std_unc=sigma/(math.sqrt(n_sites))
                final_error=sigma/(math.sqrt(n_sites))
                # t_score=t.ppf(1-0.025, n_sites-1)
                # print(t_score)
                # final_upper=mu+(t_score*final_error)   
                # final_lower=mu-(t_score*final_error)
                final_upper=mu+final_error
                final_lower=mu-final_error
                
                ##Continues with convergence plot scrambling
                while counter<countlim:    
                    xdata=[]
                    yreaddata=[]
                    ydata=[]
                    ystdevdata=[]
                    ystderr=[]
                    exp_unc_upper=[]
                    exp_unc_lower=[]
                    unc_upper=[]
                    unc_lower=[]
                    ultraExp_unc_upper=[]
                    ultraExp_unc_lower=[]
                    n=[]
                    old_n=-2
                    convergence_counter=[]  #List of all consecutive points that have converged.
                    selecteddata=data
                    #Scramble the order of the data
                    # st_time = time()
                    random.shuffle(selecteddata)
                    for n, row in enumerate(selecteddata,1):
                        if n>1:
                            cum_sum=cum_sum+row
                            cum_mean=cum_sum/n
                            ydata.append(cum_mean)
                        elif n==1:
                            cum_sum=row
                            cum_mean=cum_sum
                        xdata.append(n)
                        yreaddata.append(row)
                        cum_stdev=numpy.std(yreaddata)
                        sqrtn=math.sqrt(n)
                        cum_sterr=cum_stdev/sqrtn
                        ystdevdata.append(cum_stdev)
                        ystderr.append(cum_sterr)
                        
                        #Define the multiplying factor for the error bands, according to standard deviations. i.e. 1=67% confidence, etc. The smaller the better. 
                        mult_fac=1
                        
                        #Define the current upper and lower error band
                        upper_errorband=cum_mean+mult_fac*cum_sterr
                        lower_errorband=cum_mean-mult_fac*cum_sterr
                        
                        # #If the "true" roughness value is inbetween the error bands, save the value
                        # if lower_errorband<mu<upper_errorband:
                        #     if convergence_counter==[]:
                        #         convergence_counter=[n]
                        #     elif old_n==n-1 and len(convergence_counter)<=convergence_pointlimit:
                        #         convergence_counter.append(n)
                        #     elif old_n==n-1 and len(convergence_counter)>convergence_pointlimit:
                        #         pass
                        #     else:    
                        #         convergence_counter=[n]
                        #     old_n=n    
                        #     if len(convergence_counter)==convergence_pointlimit:
                        #         convergence_pointlist.append(convergence_counter[0])
                        #         if do_iterations==False:
                        #             print('\n %i Convergence Point Recorded. n ='%(convergence_pointlimit),convergence_counter[0])
                        #             print(' Iteration number =',counter+1)
                                    
                        #If the "true" roughness value's error bands are inbetween the current error bands, save the value
                        if lower_errorband<final_lower<upper_errorband or lower_errorband<final_upper<upper_errorband:
                            if convergence_counter==[]:
                                convergence_counter=[n]
                            elif old_n==n-1 and len(convergence_counter)<=convergence_pointlimit:
                                convergence_counter.append(n)
                            elif old_n==n-1 and len(convergence_counter)>convergence_pointlimit:
                                pass
                            else:    
                                convergence_counter=[n]
                            old_n=n    
                            if len(convergence_counter)==convergence_pointlimit:
                                convergence_pointlist.append(convergence_counter[0])
                                if do_iterations==False:
                                    print('\n %i Convergence Point Recorded. n ='%(convergence_pointlimit),convergence_counter[0])
                                    print(' Iteration number =',counter+1)

                    cum_sum=0
                    if convergence_pointlist==[]:
                        pass
                    else:
                        if first_or_last_convergence=='first':
                            last_convergence_points.append(convergence_pointlist[0])
                        elif first_or_last_convergence=='last':
                            last_convergence_points.append(convergence_pointlist[-1])
                        else:
                            print('please select either "first" or "last" for first_or_last_convergence parameter')
                            exit()
                        convergence_pointlist=[]
                        

                        if plot_graph==True:
                            suptitle=u'n= %.3g, \u03bc= %.3g, \u03C3= %.3g'%(n_sites, mu, sigma)
                            
                            ax1=plt.subplot(1,1,1)
                            plt.plot(xdata[1:],ydata[:], label='Cumulative Mean')
                            
                            #plt.plot(xdata[1:],selecteddata[1:], 'r', alpha=0.1, label='Data')
                            
                            plt.xlim(0,(int(n_sites/50))*50+50)
                            
                            ##Set y scale limits
                            # ymin, ymax = plt.ylim()
                            # yrange=ymax-ymin
                            # # if ((mu-ymin)/mu)>((ymax-mu)/mu):     ##Sets Y limit to accommdate data
                            # #     plt.ylim(ymin,(mu+mu-ymin))
                            # # elif mu-(ymax-mu)<0:
                            # #     plt.ylim(0,ymax)
                            # # else:
                            # #     plt.ylim(mu-(ymax-mu),ymax)
                            
                            # ## OR
                            # plt.ylim(mu-allowed_err,mu+allowed_err)   ##Sets y limit to just the mean+error
                            # ##OR
                            plt.ylim(mu-allowed_err*valuescale,mu+allowed_err*valuescale) ##Sets the ylimit to mean+error + extra room for data.
                            # ymin, ymax = plt.ylim()
                            
                            plt.axhline(y=mu,c='green', alpha=0.4, label='Sample Mean')
                            # upper_mu_err=[final_upper]*n_sites
                            # lower_mu_err=[final_lower]*n_sites
                            # plt.fill_between(xdata[:], upper_mu_err, lower_mu_err, alpha=0.2, color='g')
                            plt.axhline(y=final_upper,c='green', linestyle=':', alpha=0.3)
                            plt.axhline(y=final_lower,c='green', linestyle=':', alpha=0.3)
                            # plt.axhline(y=ydata[250],c='orange', alpha=0.3)
                            
                            plt.xlabel('Number of Measurements, n')
                            plt.ylabel(axislabel)
                            # plt.suptitle('Convergence of the Mean for '+sampletitles[sampledataset_number-1]+'\n'+b, fontsize=14)
                            plt.suptitle(samplename+', '+b+' - '+str(counter+1), fontsize=14)
                            plt.title(suptitle)
                            
                            for n, row in enumerate(ydata, 1):
                                unc_upper.append(row+ystderr[n])    #upper band of 1u uncertainty band representing ~68% data
                                unc_lower.append(row-ystderr[n])    #lower band of 1u uncertainty band representing ~68% data
                                
                                #Uses the t-score to expand the error bands  
                                t_score_95=(t.ppf(1-0.025, n-1))
                                t_score_99_7=(t.ppf(1-0.0015, n-1))
                                exp_unc_upper.append(row+(t_score_95*ystderr[n]))    #upper band of 2u uncertainty band representing 95% data
                                exp_unc_lower.append(row-(t_score_95*ystderr[n]))    #lower band of 2u uncertainty band representing 95% data
                                ultraExp_unc_upper.append(row+(t_score_99_7*ystderr[n]))   #upper band of 3u uncertainty band representing 99.7% data
                                ultraExp_unc_lower.append(row-(t_score_99_7*ystderr[n]))   #lower band of 3u uncertainty band representing 99.7% data
                            alpha_1=0.08
                            alpha_2=0.08
                            alpha_3=0.08
                            plt.fill_between(xdata[1:], unc_upper, unc_lower, alpha=alpha_1, color='b')
                            plt.fill_between(xdata[1:], exp_unc_upper, exp_unc_lower, alpha=alpha_2, color='b')
                            plt.fill_between(xdata[1:], ultraExp_unc_upper, ultraExp_unc_lower, alpha=alpha_3, color='b')
                            
                            ##Creating the legend
                            # where some data has already been plotted to ax
                            handles, labels = ax1.get_legend_handles_labels()
                            # manually define new patches for the different bands 
                            u1_band = mpatches.Patch(color='b', alpha=alpha_1+alpha_2+alpha_3, label='68% Coverage Probability')
                            u2_band = mpatches.Patch(color='b', alpha=alpha_2+alpha_3, label='95% Coverage Probability')
                            u3_band = mpatches.Patch(color='b', alpha=alpha_3, label='99.7% Coverage Probability')
                            handles.append(u1_band) 
                            handles.append(u2_band) 
                            handles.append(u3_band) 
                            plt.legend(handles=handles, prop={"size":9}, loc=1)
                            
                            plt.show()
                            plt.pause(5)
                            if savefile==True:
                                plt.savefig('Mean+err - '+samplename +' - '+ b +' - '+str(counter+1)+'.png') 
                            plt.close()
                        counter+=1
                    
                if do_iterations==False:
                    print('\nLast Convergence points: ', last_convergence_points) 
                
                
                """Do iterations of 20,000 reorders and plot mean/max difference"""
                if do_iterations==True:
                    while counter>=countlim and counter>0 and counter<iterations:
                        n=[]
                        ydata=[]
                        old_n=-2
                        yreaddata=[]
                        convergence_counter=[]
                        selecteddata=data
                        random.shuffle(selecteddata)
                        
                        convergence_points, counter = convergence_point_finder(selecteddata, mu, final_lower, final_upper, convergence_pointlimit, last_convergence_points, counter, first_or_last_convergence=first_or_last_convergence)
   
                    if counter==iterations:
                        mean_convergence=numpy.mean(convergence_points)
                        max_convergence=max(convergence_points)
                        # percentile_80=numpy.percentile(convergence_points, 80)
                        # percentile_85=numpy.percentile(convergence_points, 85)
                        percentile_90=numpy.percentile(convergence_points, 90)
                        percentile_95=numpy.percentile(convergence_points, 95)
                        percentile_97=numpy.percentile(convergence_points, 97)
                        # percentile_997=numpy.percentile(convergence_points, 99.7)
                        print(samplename+' - '+ b)
                        print('Mean convergence after %i iterations = %i'%(iterations, mean_convergence))
                        # print('80th Percentile after %i iterations = %i'%(iterations, percentile_80))
                        # print('85th Percentile after %i iterations = %i'%(iterations, percentile_85))
                        print('90th Percentile after %i iterations = %i'%(iterations, percentile_90))
                        print('95th Percentile after %i iterations = %i'%(iterations, percentile_95))
                        # print('97th Percentile after %i iterations = %i'%(iterations, percentile_97))
                        # print('99.7th Percentile after %i iterations = %i'%(iterations, percentile_997))
                        print('Max convergence after %i iterations = %i'%(iterations, max_convergence)+'\n')
                        # print (time() - st_time)
                        if first_or_last_convergence=='first':
                            tagline = '- %i iterations, %i points - first convergence point'%(iterations, convergence_pointlimit)
                        elif first_or_last_convergence=='last':
                            tagline = '- %i iterations, %i points - last convergence point'%(iterations, convergence_pointlimit)
                        histogram(last_convergence_points, column_number, 1, samplename, tagline, n_sites_original, other_data=[True,percentile_95,percentile_97], first_or_last=first_or_last_convergence, savefile=True)
                        counter=0
                    if reverse_calc_accuracy==True:
                        global overlap_amounts
                        reverse_counter_lim=iteration_no
                        
                        print(percentile_90)
                        convergence_point_to_test=int(percentile_90)
                        print(convergence_point_to_test)
                        overlap_amounts=[]
                        while reverse_counter < reverse_counter_lim:
                            selecteddata=data
                            random.shuffle(selecteddata)
                            
                            ##Calculate the overlap here.
                            overlap_amounts, reverse_counter = confidence_of_overlap_finder(selecteddata, mu, std_unc, convergence_point_to_test, overlap_amounts, reverse_counter)
                            
                        if reverse_counter==iterations:
                            mean_overlap=numpy.mean(overlap_amounts)
                            max_overlap=max(overlap_amounts)
                            # ov_percentile_80=numpy.percentile(overlap_amounts, 80)
                            # ov_percentile_85=numpy.percentile(overlap_amounts, 85)
                            ov_percentile_90=numpy.percentile(overlap_amounts, 90)
                            ov_percentile_95=numpy.percentile(overlap_amounts, 95)
                            ov_percentile_97=numpy.percentile(overlap_amounts, 97)
                            # ov_percentile_997=numpy.percentile(overlap_amounts, 99.7)
                            median_overlap=numpy.percentile(overlap_amounts, 50)
                            print(samplename+' - '+ b)
                            print('Mean overlap after %i iterations = %f'%(iterations, mean_overlap))
                            print('Median overlap after %i iterations = %f'%(iterations, median_overlap))
                            # print('80th Percentile after %i iterations = %i'%(iterations, ov_percentile_80))
                            # print('85th Percentile after %i iterations = %i'%(iterations, ov_percentile_85))
                            print('90th Percentile after %i iterations = %f'%(iterations, ov_percentile_90))
                            # print('95th Percentile after %i iterations = %f'%(iterations, ov_percentile_95))
                            # print('97th Percentile after %i iterations = %f'%(iterations, ov_percentile_97))
                            # print('99.7th Percentile after %i iterations = %i'%(iterations, ov_percentile_997))
                            
                            print('Max overlap after %i iterations = %f'%(iterations, max_overlap)+'\n')
                            # print (time() - st_time)
                            
                            tagline = '- %i iterations - confidence at %i measurements'%(iterations, convergence_point_to_test)
                            histogram(overlap_amounts, column_number, 1, samplename, tagline, convergence_point_to_test, other_data=[True, mean_overlap, median_overlap], first_or_last='Overlap', savefile=True)
                            reverse_counter=0
    print ('Graphs Completed') 
    
def normality_test():
    """Runs through various statistical normality tests for all samples to see how close they fit to a Gaussian distribution"""
    data_dict=all_clean_data
    for x, dataset in enumerate(dataset_names):   #Selecting which dataset    
        sampletitles=findtitles(roots[x])
        print(str(dataset))
        for samplename in sampletitles:
            print("   Sample "+str(samplename))
            for column_number, rough_par in enumerate(roughness_par_names):
                if column_number==0:
                    pass
                else:
                    roughness_data=data_dict[dataset][samplename][rough_par]
                    
                    ###Choose which of the 3 normality tests to run: (https://machinelearningmastery.com/a-gentle-introduction-to-normality-tests-in-python/)
                    ##Shapiro normality test
                    # stat, p = shapiro(roughness_data)
                    ##D’Agostino’s K^2 Test normality test
                    # stat, p = normaltest(roughness_data)
                    
                    
                    # print('      '+str(rough_par))#+' - Statistics=%.3f, p=%.3f' % (stat, p))
                    # ## Interpret
                    # alpha = 0.05
                    # if p > alpha:
                    #  	print('         Sample looks Gaussian (fail to reject H0)')
                    # else:
                    #  	print('         Sample does not look Gaussian (reject H0)')
                        
                    ##Anderson normality test
                    result = anderson(roughness_data)
                    # print('Statistic: %.3f' % result.statistic)
                    p = 0
                    print(str(rough_par))
                    for i in range(len(result.critical_values)):
                     	sl, cv = result.significance_level[i], result.critical_values[i]
                     	if result.statistic < result.critical_values[i]:
                             print('%.3f: %.3f, data looks normal (fail to reject H0)' % (sl, cv))
                     	else:
                    		# print('%.3f: %.3f, data does not look normal (reject H0)' % (sl, cv))
                            pass
                    # time.sleep(0.1)
    print('Complete')

                  
##Choose which function to run:
    
# all_samples_IQRMethod(4)
# mean_convergence(root_1, include_offset=False, do_iterations=False, correct_order=False, savefile=False)
# normality_test()
mean_convergence_2(root_1, first_or_last_convergence='last', no_of_converged_points=20, plot_graph=False, do_iterations=True, iteration_no=20000, reverse_calc_accuracy=True, savefile=False)

                        