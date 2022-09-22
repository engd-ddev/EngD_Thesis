import os, numpy
from scipy.stats import pearsonr, spearmanr
# from timeit import default_timer as timer

from Blinds_Core_Functions_0030_test import findtitles, boxplot, three_combined_uncertainty_plotter, combined_uncertainty_plotter, histogram, read_all_data, append_value, samplecomparison_histogram, bootstrap_confidence_intervals, random_sampling_within_confidence_intervals, dictionarystats, uncertainty_plotter, mean_and_uncertainty_plot, read_bendstrain, comparison_plot #readstats_readlastframe, remove_chosen_outliers
# import Core_Functions_0001 #Imports functions and the modules. 
from Dataset_file_locations import get_dataset_info 
from time import time
# first_start=timer()
# start = timer()

tablefilename = 'fixed_zero_stats_04.txt'   #The .txt file name to read data from
txtfilepaths=[]
dataset_names=[]
root_1='E:/2.SiC 140um Blinds/1.ILC/1.Section 1/'
root_2='E:/2.SiC 140um Blinds/3.ILC Repeat/1.Section 1/'

##Get txtfilepaths for the different samples, the graph saving file root location and the dataset name.
"""dataset_number=3 corresponds to the "1st Dataset" 
    dataset_number=4 corresponds to the "2nd Dataset"
    For purposes of naming, the 1st and 2nd datasets are 
    called 3rd and 4th respectively in the code due to obsolete 
    datasets being previously used. 
    """
txtfilepathways, graphsave_root, dataset_name, roughness_par_names, roughness_par_labels, roughness_par_titles = get_dataset_info(Method_number=1, dataset_number=4)

print(dataset_name)

dataset_names.append(dataset_name)
txtfilepaths.append(txtfilepathways)
roots = [root_2]

txtfilepathways, graphsave_root, dataset_name, roughness_par_names, roughness_par_labels, roughness_par_titles = get_dataset_info(Method_number=1, dataset_number=3)
print(dataset_name)
dataset_names.append(dataset_name)
txtfilepaths.append(txtfilepathways)
roots = [root_2, root_2]


##Read data into dictionaries. Outlier data is not exported and as such an underscore is used to ignore that output.
all_data, all_clean_data, _ = read_all_data(roots, txtfilepaths, dataset_names, 'fixed_zero_stats_04.txt')
bendstrain_data=read_bendstrain()
bendstrain_stats=dictionarystats(bendstrain_data)

del dataset_name, txtfilepathways
graphsave_loc = graphsave_root+'Single Histograms/'
os.chdir(graphsave_loc)  

def bendstrain_graphs(bootstrap_method=False):
    """A function that plots the bend strains graphs and their uncertainties according to either assuming 
    a t-distribution or using the bootstrap method. The t-distribution is advised/preferred."""
    sampletitles=findtitles(roots[0])
    tagline='A03 Anomaly Removed - t-score'
    os.chdir('F:/3.SiC Blinds/')
    if bootstrap_method==True:
        tagline='A03 Anomaly Removed - Bootstrap'
        uncertainty_plotter(bendstrain_data, 0, 0, sampletitles, 'Bend Strain Results', tagline, savefile=True)
    else:
        mean_and_uncertainty_plot(bendstrain_data, 0, 0, sampletitles, 'Bend Strain Results', tagline, savefile=True)


def plot_all_roughness_histograms(remove_anomalies=True):
    """A function that plots histograms of each area roughness parameter to show the 
       distribution of roughnesses of the monofilaments"""
    # first_start=timer()
    if remove_anomalies==True:
        data_dict=all_clean_data
        tagline=' - no anomalies'
    elif remove_anomalies==False:
        data_dict=all_data
        tagline=' - unfiltered' 
    # start = timer()
    for x, dataset in enumerate(dataset_names):   #Selecting which dataset    
        sampletitles=findtitles(roots[x])
        graphsave_loc = 'F:/3.SiC Blinds/1.ILC/Graphs/4.1 - 4th Dataset - Method 1 Graphs/'+'Single Histograms/'
        os.chdir(graphsave_loc)
        for samplename in sampletitles:
            title=samplename+' - '+dataset        #the title for the histogram
            for column_number, rough_par in enumerate(roughness_par_names):
                if column_number==0:
                    pass
                else:
                    roughness_data=data_dict[dataset][samplename][rough_par]
                    n_sites=len(roughness_data)
                    histogram(roughness_data, column_number, 1, title, tagline, n_sites, savefile=True)
                    boxplot(roughness_data, column_number, 1, title, tagline, n_sites, savefile=True)
    # print("Line completed in:", timer()-start)
    # print("Total time to run: ",timer()-first_start)

def all_samplecomparison_histograms(remove_anomalies=True):
    """A function that plots all the sample-comparison histogram graphs to visually examine which samples are rougher/smoother"""
    global roughness_data
    if remove_anomalies==True:
        data_dict=all_clean_data
        tagline=' - no anomalies'
    elif remove_anomalies==False:
        data_dict=all_data
        tagline=' - unfiltered' 
    graphsave_loc= graphsave_root+'Sample Comparison Histograms/'
    for x, dataset in enumerate(dataset_names):   #Selecting which dataset    
        sampletitles=findtitles(roots[x])
        os.chdir(graphsave_loc)
        for column_number, rough_par in enumerate(roughness_par_names):
            if column_number==0:    #Skips the Distance column
                pass
            else:
                ##Creating the dictionary of roughness data for each roughness parameter for each dataset.
                roughness_data={}
                for samplename in sampletitles:
                    rough_para_data=data_dict[dataset][samplename][rough_par]
                    append_value(roughness_data, samplename, rough_para_data)    
                
                #Plot joint histogram once the data has been collected into a single dictionary
                samplecomparison_histogram(roughness_data, column_number, 1, sampletitles, dataset, tagline, savefile=True)    

def all_mean_and_uncertainty_plots(remove_anomalies=True, bootstrap_method=False, table=False):
    """A function that plots the means and uncertainties for the roughnesses of the samples either using the bootstrap 
    method to estimate uncertainties or not."""
    if remove_anomalies==True:
        data_dict=all_clean_data
        if bootstrap_method==True:
            tagline='- bootstrap no anomalies'
        else:
            tagline=' - no anomalies'
    elif remove_anomalies==False:
        data_dict=all_data
        if bootstrap_method==True:
            tagline='- bootstrap unfiltered'
        else:
            tagline=' - unfiltered' 
    graphsave_loc= graphsave_root+'Sample Comparison Mean and Uncertainties/'
    for x, dataset in enumerate(dataset_names):   #Selecting which dataset    
        sampletitles=findtitles(roots[x])
        os.chdir(graphsave_loc)
        ##For each roughness parameter, plot the graph.
        for column_number, rough_par in enumerate(roughness_par_names):
            if column_number==0:
                pass
            else:
                roughness_data={}
                for samplename in sampletitles:
                    rough_para_data=data_dict[dataset][samplename][rough_par]
                    append_value(roughness_data, samplename, rough_para_data)
                if bootstrap_method==True:
                    uncertainty_plotter(roughness_data, column_number, 1, sampletitles, dataset, tagline, savefile=True)
                else:
                    if table==False:
                        mean_and_uncertainty_plot(roughness_data, column_number, 1, sampletitles, dataset, tagline, savefile=True, only_errors=False)
                    else:
                        means, err, st_unc, st_dev, numbers = mean_and_uncertainty_plot(roughness_data, column_number, 1, sampletitles, dataset, tagline, savefile=False, only_errors=True)
                        print(rough_par)
                        for i,j in enumerate(means):
                            print(str(j) + ' '+str(st_dev[i])+' '+ str(numbers[i]))

def dataset_comparison_mean_and_uncertainty_plots(remove_anomalies=True, bootstrap_method=True, table=False):
    """A function that plots the means and uncertainties for the roughnesses of the samples for two datasets either using the bootstrap 
    method to estimate uncertainties or not. Please ensure only two datasets are read in at the top of the python file."""
    if not len(dataset_names)==2:
        print('Please read in only two datasets for comparisons. (See top of the python file).')
        return
    if remove_anomalies==True:
        data_dict=all_clean_data
        if bootstrap_method==True:
            tagline='- bootstrap no anomalies, 60k runs'
        else:
            tagline=' - no anomalies'
    elif remove_anomalies==False:
        data_dict=all_data
        if bootstrap_method==True:
            tagline='- bootstrap unfiltered'
        else:
            tagline=' - unfiltered' 
    graphsave_loc= 'F:/3.SiC Blinds/1.ILC/Graphs/Dataset Comparison Graphs/'
    
    sampletitles=findtitles(roots[0])
    os.chdir(graphsave_loc)
    
    dataset_1=dataset_names[0]
    print('Dataset 1 is: ', dataset_1)
    dataset_2=dataset_names[1]
    print('Dataset 2 is: ', dataset_2)
    
    ##For each roughness parameter, plot the graph.
    for column_number, rough_par in enumerate(roughness_par_names):
        if column_number==0 or column_number==4 or column_number==5:
            print(roughness_par_names[column_number]+' skipped.')
            pass
        else:
            dataset_1_roughness_data={}
            dataset_2_roughness_data={}
            st_time = time()
            for samplename in sampletitles:
                
                #Create roughness data dict for dataset 1
                rough_para_data_1=data_dict[dataset_1][samplename][rough_par]
                append_value(dataset_1_roughness_data, samplename, rough_para_data_1)
                #Create roughness data dict for dataset 2
                rough_para_data_2=data_dict[dataset_2][samplename][rough_par]
                append_value(dataset_2_roughness_data, samplename, rough_para_data_2)
            
            dataset_comparison_title=dataset_1+' and '+dataset_2
            two_dataset_labels=[dataset_1, dataset_2]
            two_dataset_labels=['1st Dataset', '2nd Dataset']
            combined_uncertainty_plotter(dataset_1_roughness_data, dataset_2_roughness_data, column_number, 1, sampletitles, two_dataset_labels, dataset_comparison_title, tagline, savefile=True)
            print(roughness_par_names[column_number]+' completed in: '+ "%.2f s" %(time() - st_time))
            # if bootstrap_method==True:
            #     # uncertainty_plotter(roughness_data, column_number, 1, sampletitles, dataset, tagline, savefile=True)
            #     pass
            # else:
            #     if table==False:
            #         pass
            #         # mean_and_uncertainty_plot(roughness_data, column_number, 1, sampletitles, dataset, tagline, savefile=True, only_errors=False)
            #     else:
            #         means, err, st_unc, st_dev, numbers = mean_and_uncertainty_plot(roughness_data, column_number, 1, sampletitles, dataset, tagline, savefile=False, only_errors=True)
            #         print(rough_par)
            #         for i,j in enumerate(means):
            #             print(str(j) + ' '+str(st_dev[i]), ' '+str(numbers[i]))

def dataset_comparison_mean_and_uncertainty_plots_repeated_A09(remove_anomalies=True, bootstrap_method=True, table=False):
    """Same function as above, but plots the A09 repeat too in another colour to show repeatability"""
    if not len(dataset_names)==3:
        print('Please read in only three datasets for comparisons. (See top of the python file).')
        return
    if remove_anomalies==True:
        data_dict=all_clean_data
        if bootstrap_method==True:
            tagline='- bootstrap no anomalies, 60k runs-A09 Repeat'
        else:
            tagline=' - no anomalies'
    elif remove_anomalies==False:
        data_dict=all_data
        if bootstrap_method==True:
            tagline='- bootstrap unfiltered'
        else:
            tagline=' - unfiltered' 
    graphsave_loc= 'F:/3.SiC Blinds/1.ILC/Graphs/Dataset Comparison Graphs/'
    
    sampletitles=findtitles(roots[0])
    os.chdir(graphsave_loc)
    
    dataset_1=dataset_names[0]
    dataset_2=dataset_names[1]
    dataset_3=dataset_names[2]
    
    ##For each roughness parameter, plot the graph.
    for column_number, rough_par in enumerate(roughness_par_names):
        if column_number==0 or column_number==4 or column_number==5:
            print(roughness_par_names[column_number]+' skipped.')
            pass
        else:
            dataset_1_roughness_data={}
            dataset_2_roughness_data={}
            dataset_3_roughness_data={}
            st_time = time()
            for samplename in sampletitles:
                #Create roughness data dict for dataset 1
                rough_para_data_1=data_dict[dataset_1][samplename][rough_par]
                append_value(dataset_1_roughness_data, samplename, rough_para_data_1)
                #Create roughness data dict for dataset 2
                rough_para_data_2=data_dict[dataset_2][samplename][rough_par]
                append_value(dataset_2_roughness_data, samplename, rough_para_data_2)
            
            #Appends only the A09 data from Dataset 3 as the "repeat" sample.
            rough_para_data_3=data_dict[dataset_3]['A09'][rough_par]
            append_value(dataset_3_roughness_data, 'A09', rough_para_data_3)
            
            dataset_comparison_title= dataset_1+' and '+ dataset_2
            three_dataset_labels=['1st Dataset', '2nd Dataset', 'Repeated A09']
            
            three_combined_uncertainty_plotter(dataset_1_roughness_data, dataset_2_roughness_data, dataset_3_roughness_data, column_number, 1, sampletitles, three_dataset_labels, dataset_comparison_title, tagline, savefile=True)
            print(roughness_par_names[column_number]+' completed in: '+ "%.2f s" %(time() - st_time))


def roughness_bendstrain_comparison(remove_anomalies=True, bootstrap_method=True, savefile=False):
    """ A function that plots the roughnesses against the bend strain results to see if there is a correlation. 
    It also calculates the PMCC and Spearman's rank Correlation coefficient to quantify the relationships. 
    
    If bootstrap_method is True:
        The function will plot [and save] the graphs according to the bootstrap method for uncertainty measurement for roughnesses and
        print out the correlation results according to the bootstrapped means.
    If false,:
        The function will plot [and save] the graphs using the NPL standard uncertainty measurement, but ALSO randomly sample from this 
        uncertainty to transfer the uncertainties to the correlation results, which will print the means with their 95% CI lower and upper bands. 
    """
    global x_data, x_err, x_st_unc, y_data, y_err, y_st_unc
    sampletitles=findtitles(roots[0])
    if remove_anomalies==True:
        data_dict=all_clean_data
        if bootstrap_method==True:
            tagline='- bootstrap no anomalies'
        else:
            tagline=' - no anomalies - NPL Method'
    elif remove_anomalies==False:
        data_dict=all_data
        if bootstrap_method==True:
            tagline='- bootstrap unfiltered'
        else:
            tagline=' - unfiltered' 
    graphsave_loc= graphsave_root+'Bendstrain Comparison/'
    
    # x_data, x_err = uncertainty_plotter(bendstrain_data, 0, 0, sampletitles, 'Bendstrain Results', tagline, savefile=False, only_errors=True)
    
    x_data, x_err, x_st_unc = mean_and_uncertainty_plot(bendstrain_data, 0, 0, sampletitles, 'Bendstrain Results', tagline, savefile=False, only_errors=True)
    # x_data.pop(7)                   #Remove sample A08
    # x_err=numpy.delete(x_err, 7, 1) #Remove sample A08
    # print(x_err)

    
    xlabel = 'Bend Strain (%)'
    for x, dataset in enumerate(dataset_names):   #Selecting which dataset    
        sampletitles=findtitles(roots[x])
        os.chdir(graphsave_loc)
        # sampletitles.pop(7)             #Remove sample A08
        ##For each roughness parameter, plot the graph.
        for column_number, rough_par in enumerate(roughness_par_names):
            if column_number==0:
                pass
            else:
                roughness_data={}
                for samplename in sampletitles:
                    rough_para_data=data_dict[dataset][samplename][rough_par]
                    append_value(roughness_data, samplename, rough_para_data)
                if bootstrap_method==True:
                    y_data, y_err = uncertainty_plotter(roughness_data, column_number, 1, sampletitles, dataset, tagline, savefile=False, only_errors=True)
                    ylabel=roughness_par_labels[column_number]
                    figure_title = 'Correlation Plot of Bend Strain vs ' + roughness_par_titles[column_number]
                    figure_name = figure_title+' '+tagline
                    if savefile==True:
                        comparison_plot(x_data, x_err, y_data, y_err, xlabel, ylabel, figure_title, figure_name, savefile=True)
                    else:
                        comparison_plot(x_data, x_err, y_data, y_err, xlabel, ylabel, figure_title, figure_name, savefile=False)
                    
                    #Calculate Pearson's Correlation Coefficient
                    corr, _ = pearsonr(x_data, y_data)
                    print(rough_par)
                    print('Pearsons correlation: %.3f' %corr)
                    
                    #Calculate Spearman's Rank correlation Coefficient
                    corr1, _ = spearmanr(x_data, y_data)
                    print('Spearmans correlation: %.3f' % corr1)
                else:
                    y_data, y_err, y_st_unc = mean_and_uncertainty_plot(roughness_data, column_number, 1, sampletitles, dataset, tagline, savefile=False, only_errors=True)
                    ylabel=roughness_par_labels[column_number]
                    figure_title = 'Correlation Plot of Bend Strain vs ' + roughness_par_titles[column_number]
                    figure_name = figure_title+' '+tagline
                    if savefile==True:
                        comparison_plot(x_data, x_err, y_data, y_err, xlabel, ylabel, figure_title, figure_name, savefile=True)
                    est, ci = random_sampling_within_confidence_intervals(x_data, x_st_unc, y_data, y_st_unc, pearsonr, [2.5, 97.5], number_of_samples=10, runs=10000)
                    print(rough_par)
                    print("     Pearson mean and CI:         {:.3g} [{:.3g}, {:.3g}]".format(est, ci[0], ci[1]))
                    
                    est, ci = random_sampling_within_confidence_intervals(x_data, x_st_unc, y_data, y_st_unc, spearmanr, [2.5, 97.5], number_of_samples=10, runs=10000)
                    print("     Spearmans Rank mean and CI:   {:.3g} [{:.3g}, {:.3g}]".format(est, ci[0], ci[1]))
                    
def conf_int_bootstrapmethod(plot_graphs=False):
    """A function that estimates the 95 % confidence interval for the roughness data by sampling by replacement the original data 
    and iterating it 1000 times.
    """
    data_dict=all_data  
    graphsave_loc= graphsave_root+'Test/'
    tagline='- test' 
    if plot_graphs==False:
        for x, dataset in enumerate(dataset_names):   #Selecting which dataset    
            sampletitles=findtitles(roots[x])
            os.chdir(graphsave_loc)
            print('\n-=-=-=-=-=-=-=-=-=-=-=-=-'+dataset+'=-=-=-=-=-=-=-=-=-=-=-=-=--=-')
            
            for samplename in sampletitles:
                print(' '+samplename)
                for column_number, rough_par in enumerate(roughness_par_names):
                    if column_number==0:
                        pass
                    else:
                        print('   '+rough_par)
                        roughness_data=data_dict[dataset][samplename][rough_par]
                        print("     Mean: {:.3g}".format(numpy.mean(roughness_data)))
                        est, ci = bootstrap_confidence_intervals(roughness_data, column_number, numpy.mean, [2.5, 97.5])
                        print("     Bootstrapped mean and CI: {:.3g} [{:.3g}, {:.3g}]".format(est, ci[0], ci[1]))
    if plot_graphs==True:       
        for x, dataset in enumerate(dataset_names):   #Selecting which dataset    
            sampletitles=findtitles(roots[x])
            print(sampletitles)
            os.chdir(graphsave_loc)
            ##For each roughness parameter, plot the graph.
            for column_number, rough_par in enumerate(roughness_par_names):
                if column_number==0:
                    pass
                else:
                    roughness_data={}
                    for samplename in sampletitles:
                        rough_para_data=data_dict[dataset][samplename][rough_par]
                        append_value(roughness_data, samplename, rough_para_data)
                    uncertainty_plotter(roughness_data, column_number, 1, sampletitles, dataset, tagline, savefile=False)

def print_number_of_frames():
    """Prints the number of total frames - anomalous frames - filtered frames for all datasets"""                    
    for dataset in dataset_names:
        print(dataset)
        for sample_no in ['A01', 'A02','A03','A04','A05','A06','A07','A08','A09','A10']:
            print(sample_no+': '+str(len(all_data[dataset][sample_no]['Distance (um)']))+ ' '+
                  str(len(all_data[dataset][sample_no]['Distance (um)'])-len(all_clean_data[dataset][sample_no]['Distance (um)']))+
                  ' '+str(len(all_clean_data[dataset][sample_no]['Distance (um)'])))

##Choose which function to run: 
# bendstrain_graphs(bootstrap_method=False)                   
# plot_all_roughness_histograms(remove_anomalies=False)
# all_samplecomparison_histograms(remove_anomalies=False)
# conf_int_bootstrapmethod(plot_graphs=True)

all_mean_and_uncertainty_plots(remove_anomalies=False, bootstrap_method=False, table=True)
# dataset_comparison_mean_and_uncertainty_plots(remove_anomalies=True, bootstrap_method=True, table=True)
# dataset_comparison_mean_and_uncertainty_plots_repeated_A09(remove_anomalies=True, bootstrap_method=True, table=False)
# print_number_of_frames()
# roughness_bendstrain_comparison(remove_anomalies=True, bootstrap_method=False, savefile=False)

