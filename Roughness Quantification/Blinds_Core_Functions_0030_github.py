import os, fnmatch, csv, numpy, math, seaborn as sns, matplotlib.pyplot as plt
from scipy.stats import norm, t
from scipy.stats import pearsonr, spearmanr
import matplotlib.patches as mpatches



"""==========BASIC LEVEL FUNCTIONS: Simple plotting, simple functions used numerous times=========="""
def findtitles(root_path):           #Gets titles of the samples, removing "?." from the beginning
    """
    Finds the titles of the samples in the 'root' folder path to allow for easy labelling of graphs. 
    Be aware: Assumes all the folders within the root folder are the names of the samples.
    """
    sampletitles=[]
    os.chdir(root_path)
    foldertitles=next(os.walk('.'))[1]
    for title in foldertitles:
        if fnmatch.fnmatch(title,'A??'):
            sampletitles.append(title)
        else:
            pass
    return sampletitles

def append_value(dict_obj, key, value):
    """
    A function that adds a value to a dictionary's key, checking first whether that key exists first or not, 
    and then adds the value accordingly as a new key:value pair, or as a new value
    """
    # Check if key exist in dict or not
    if key in dict_obj:
        # Key exist in dict.
        # Add the value to the key
        dict_obj[key].update(value)
    else:
        # As key is not in dict, add key-value pair
        dict_obj[key] = value

def dictionarystats(datadict):
    """
    A function that calculated the mean, median, st dev and st unc for a given dictionary of form:
        datadict["Samplename"][data]
    """
    dictionarystats={}
    for key in datadict:
        sampledata=datadict[key]
        #Count the number of datapoints
        n_measurements=int(len(sampledata))
        #Calculate its mean
        mu = numpy.mean(sampledata)    
        #Calculate the standard deviation
        sigma = numpy.std(sampledata)  
        #Calculate the median
        med = numpy.median(sampledata) 
        #Calculate u by first doing the sqrt of n
        sqrtn=math.sqrt(n_measurements)
        u = sigma/sqrtn
        append_value(dictionarystats, key, {
                            "Mean": mu,
                            "Standard Deviation": sigma,
                            "Median": med,
                            "Number of Measurements": n_measurements,
                            "Standard Uncertainty": u
                            })
    return dictionarystats

def readstats_readlastframe(tablefilename):
    """
    Reads the 'tablefilename' .txt file and imports the data from the last frame at each site, 
    to minimise chance of Loss of Contact frames being used. 
    Returns the data as a dictionary called 'data'.
    """
    global readstats
    if readstats==False:
        Distance=[]
        Avg=[]
        Ra=[]
        Rms=[]
        Skew=[]
        Kurtosis=[]
        Zmax=[]
        with open(tablefilename, 'r') as csvfile:
            plots= csv.reader(csvfile, delimiter=' ')
            line_count=0
            old_a=-1
            counter=0
            for row in plots:
                if line_count==0:
                    line_count += 1           #skips header line
                else:
                    a=int(row[0])
                    if a==old_a and counter==0:      #Reading second frame 
                        Distance.append(int(row[0]))
                        Avg.append(float(row[1]))
                        Ra.append(float(row[2]))
                        Rms.append(float(row[3]))
                        Skew.append(float(row[4]))
                        Kurtosis.append(float(row[5]))
                        Zmax.append(float(row[6]))
                        counter=1
                        
                    elif a==old_a and counter==1:  #If 3 (or more) frames exist, deletes 2nd frame's data from list. Currently on 3rd Frame.
                        #Delete last value of variables
                        del Distance[-1]
                        del Avg[-1]
                        del Ra[-1]
                        del Rms[-1]
                        del Skew[-1]
                        del Kurtosis[-1]
                        del Zmax[-1]
                        counter=2
                        
                    elif a==old_a and counter==2:   #On the 4th frame, save value. 
                        Distance.append(int(row[0]))
                        Avg.append(float(row[1]))
                        Ra.append(float(row[2]))
                        Rms.append(float(row[3]))
                        Skew.append(float(row[4]))
                        Kurtosis.append(float(row[5]))
                        Zmax.append(float(row[6]))
                        
                    else:           #Reading the first frame of a new site
                        counter=0
                        pass
                    old_a=a     #Saving the current location to compare with the next location       
            readstats=True
        data = {
            "Distance (um)": Distance,    #Distance Along Monofilament (estimate based on step size of 25um)
            "Average Height (nm)": Avg,    #Average Height of frame
            "AM, Ra/Sa (nm)": Ra,         #Arithmetic Mean Roughness
            "RMS, Rq/Sq (nm)": Rms,       #Root Mean Square Roughness
            "Rsk/Ssk": Skew,              #Skewness
            "Rku/Sku": Kurtosis,          #Kurtosis
            "Rz/Sz (nm)": Zmax            #Height Range
            }
    else:
        print("Data not read as readstats==False")
        pass    
    return data

def remove_chosen_outliers(data_dict_of_roughness_parameters, dataname, tablefilename):
    """
    Reads the .txt file and removes the manually chosen outliers (the variable) from the dataset. 
    This can be done following a filtering method, such as IQRoutliermethod to get the 'clean' dataset.
    """
    
    global chosen_outliers, path, readstats
    Distance_outliers=[]
    Avg_height_outliers=[]
    Ra_outliers=[]
    Rms_outliers=[]
    Skew_outliers=[]
    Kurtosis_outliers=[]
    Zmax_outliers=[]
    bad_row=[]
    chosen_outliers=[]
  
    if dataname=='4th Dataset - A01':
        chosen_outliers=[1200, 1400, 1450, 1775, 1800, 2100, 2175, 2450, 2475, 
                         2500, 2525, 2550, 2800, 3225, 3475, 3600, 3875, 4100, 
                         4150, 4300, 4450, 4600, 4675, 4725, 4850, 4900, 4975, 
                         5175, 5375, 5425, 5650, 5950, 5975, 6100, 6150, 6175, 
                         6300, 6325, 6350, 6400, 6525, 6575, 6625, 6775, 6800, 
                         6825, 6975, 7050, 7325, 7425, 7450, 7500, 7525, 7500, 
                         7875, 8025, 8150]                                           #Manually Selected bad frames for A1. 
        
    elif dataname=='4th Dataset - A02':
        chosen_outliers=[300, 325, 550, 1000, 1050, 1200, 1525, 1725, 1775, 
                         1850, 1875, 1950, 1975, 2100, 2300, 2425, 2450, 2475, 
                         2525, 2550, 2625, 2775, 2925, 2950, 3375, 3450, 3500, 
                         3800, 3850, 3975, 4250, 4325, 4825, 4850, 5200, 5450, 
                         6150, 7500, 7525, 7625, 7650, 7675, 7700, 7725, 7750, 
                         7775, 7800, 7825]                                            #Manually Selected bad frames for A2. 
        
    elif dataname=='4th Dataset - A03':  
        chosen_outliers=[600, 650, 750, 1075, 1600, 2125, 2275, 2725, 2750, 
                         2825, 3075, 3225, 3275, 3300, 3475, 3550, 4200, 4225, 
                         4250, 4700, 4725, 4750, 4825, 4850, 5150, 5200, 5325, 
                         5425, 5750, 5775, 5825, 7275, 7725]                                    #Manually Selected bad frames for A3. 
    
    elif dataname=='4th Dataset - A04':  
        chosen_outliers=[275, 350, 400, 575, 600, 1050, 1100, 1150, 1775, 1800, 
                         1825, 2300, 2325, 2450, 2550, 2600, 2625, 2725, 2750, 
                         2775, 2975, 3000, 3025, 3050, 3075, 3100, 3225, 3300, 
                         3450, 3525, 3775, 3800, 4450, 4725, 4925, 4950, 5000, 
                         5025, 5125, 5475, 5500, 5825, 6000, 6075, 6100, 6275, 
                         6300, 6325, 6700, 6800, 6850, 6925, 6975, 7050, 7125, 
                         7325]                                    #Manually Selected bad frames for A4.
    
    elif dataname=='4th Dataset - A05':  
        chosen_outliers=[950, 1175, 1600, 2075, 2150, 2450, 2525, 2550, 2575, 
                         2600, 2625, 2650, 2675, 2725, 2825, 2875, 2925, 2950, 
                         3075, 3100, 3225, 3250, 3325, 3425, 3450, 3525, 3575,
                         3600, 3650, 3725, 3875, 4250, 4275, 4325, 4450, 4575, 
                         4700, 4725, 5125, 5200, 5450, 6200, 6350, 6550, 6875, 
                         6950, 7350, 8075, 8200]                                    #Manually Selected bad frames for A5.
    
    elif dataname=='4th Dataset - A06':  
        chosen_outliers=[125, 150, 1750, 1775, 2025, 3775, 4150, 5700, 6100, 
                         6675, 6925]                                    #Manually Selected bad frames for A06. 
    
    elif dataname=='4th Dataset - A07':  
        chosen_outliers=[200, 525, 750, 1375, 1925, 1975, 2075, 2125, 2150, 2175, 
                         2200, 2225, 2500, 2775, 3650, 3750, 3975, 4125, 4625, 
                         4750, 4825, 5225, 5275, 5525, 5750, 5900, 6425, 7350, 
                         7375, 7750]                                    #Manually Selected bad frames for A07
    
    elif dataname=='4th Dataset - A08':  
        chosen_outliers=[75, 100, 150, 175, 225, 275, 325, 350, 450, 550, 650,
                         700, 725, 750, 775, 800, 825, 925, 950, 975, 1100, 1125,
                         1150, 1200, 1225, 1250, 1275, 1300, 1325, 1350, 1375, 
                         1400, 1450, 1475, 1500, 1525, 1575, 1600, 1675, 1700, 
                         1875, 1975, 2025, 2050, 2400, 2450, 2475, 2625, 2700,
                         2925, 3450, 3475, 3550, 3575, 3675, 3700, 3750, 3775, 
                         3825, 3925, 4175, 4275, 4600, 4625, 4650, 4675, 4800, 
                         5075, 5100, 5125, 5150, 5200, 5225, 5475]      #Manually Selected bad frames for A08. 
    
    elif dataname=='4th Dataset - A09':  
        chosen_outliers=[1450, 1475, 2975, 3075, 3125, 3225, 3500, 3550, 3575, 
                         3700, 4075, 4100, 4125, 4150, 4175, 4200, 4225, 4250, 
                         4275, 4300, 4325, 4350, 4375, 4400, 4425, 4450, 4475, 
                         4500, 4525, 4550, 4575, 4600, 4625, 4650, 4675, 4700,
                         4725, 4750, 4775, 4800, 4825, 4850, 4875, 4900, 4925,
                         4950, 4975, 5000, 5025, 5050, 5075, 5100, 5125, 5200, 
                         5225, 5250, 5275, 5300, 5325, 5350, 5375, 5400, 5425, 
                         5450, 5500, 5575, 5600, 6450, 6475, 6500, 6525, 6550,
                         6575, 6600, 6625, 6650, 6675, 6700, 6725, 6750, 6775,
                         6800, 6825, 6850, 6875, 6900, 6925, 6950, 6975, 7000,
                         7025, 7050, 7075, 7100, 7125, 7150, 7175, 7200]  #Manually Selected bad frames for A09. 
    
    
    elif dataname=='4th Dataset - A10':  
        chosen_outliers=[50, 250, 375, 500, 1750, 1775, 1950, 2050, 2125, 2650, 
                         2975, 3575, 3725, 3900, 4025, 4050, 4450, 4525, 4700, 
                         4750, 4975, 5175, 5225, 5250, 5650, 5950, 5975, 6025, 
                         6025, 6775, 6875]                                    #Manually Selected bad frames for A10. 
        
    elif dataname=='3rd Good Dataset - A01':
        chosen_outliers=[250, 300, 375, 400, 500, 625, 675, 700, 775, 925, 1000, 
                         1125, 1200, 1225, 1300, 1425, 1450, 1525, 1600, 1675, 
                         1725, 1850, 1950, 1975, 2000, 2050, 2150, 2175, 2225, 
                         2300, 2550, 2625, 2650, 2700, 2950, 3025, 3125, 3300, 
                         3325, 3350, 3575, 3675, 3875, 3975, 4200, 4225, 4275, 
                         4325, 4350, 4375, 4400, 4425, 4450, 4475, 4525, 4625, 
                         4675, 4750, 4800, 4825, 4850, 4900, 5000, 5075, 5250, 
                         5275, 5375, 5800, 5925, 5975, 6000, 6025, 6050, 6125, 
                         6150, 6225, 6275, 6350, 6400, 6550, 6575, 6600, 6650, 
                         6700, 6825, 6850, 6875, 6900, 6975, 7025, 7075, 7100, 
                         7150, 7200, 7325, 7425, 7575, 7600, 7625, 7725, 7800, 
                         7825, 7850, 8525, 8550]                                    #Manually Selected bad frames for A01. 
        
    elif dataname=='3rd Good Dataset - A02':
        chosen_outliers=[25, 175, 375, 450, 550, 625, 675, 700, 850, 875, 950, 
                         975, 1000, 1050, 1075, 1175, 1275, 1500, 1525, 1700, 
                         1725, 1750, 1825, 1900, 1950, 1975, 2025, 2075, 2175, 
                         2225, 2275, 2300, 2525, 2550, 2600, 2725, 2775, 2825, 
                         2875, 2900, 2925, 2950, 2975, 3000, 3025, 3125, 3250, 
                         3425, 3475, 3500, 3550, 3575, 3600, 3650, 3725, 3800, 
                         3825, 3850, 3975, 4225, 4625, 5075, 5400, 5550, 5875, 
                         6000, 6325, 6350, 6725, 6800, 6825, 6875, 7300, 7325, 
                         7525, 7800, 7900, 7950]                                    #Manually Selected bad frames for A02. 
        
    elif dataname=='3rd Good Dataset - A03':  
        chosen_outliers=[25, 50, 200, 225, 250, 275, 300, 325, 375, 400, 425, 
                         650, 750, 775, 800, 825, 875, 900, 975, 1025, 1050, 
                         1075, 1175, 1200, 1225, 1375, 1475, 1550, 1675, 1725, 
                         1775, 1825, 1975, 2000, 2025, 2350, 2400, 2525, 2625, 
                         2725, 2750, 2875, 2900, 2975, 3100, 3225, 3250, 3275, 
                         3300, 3425, 3475, 3550, 3600, 3700, 3725, 3750, 3850, 
                         3900, 3925, 4000, 4025, 4050, 4075, 4125, 4250, 4375, 
                         4550, 4575, 4700, 4725, 4775, 4850, 5100, 5125, 5150, 
                         5200, 5325, 5350, 5375, 5400, 5425, 5525, 5700, 5750, 
                         5800, 5825, 6000, 6100, 6300, 6450, 6475, 6500, 6825, 
                         6850, 6875, 7050, 7075, 7250, 7275, 7300, 7325, 7350, 
                         7375, 7400, 7425, 7450, 7475, 7500, 7525, 7550, 7575, 
                         7600, 7625, 7650, 7675, 7700, 7725, 7750, 7775, 7800]      #Manually Selected bad frames for A03.
    
    elif dataname=='3rd Good Dataset - A04':  
        chosen_outliers=[275, 350, 425, 500, 550, 575, 600, 1000, 1025, 1075, 
                         1550, 1725, 1975, 2000, 2225, 2250, 2275, 2300, 2325, 
                         2350, 2375, 2425, 2450, 2475, 2525, 2700, 2725, 2750, 
                         2775, 2800, 2875, 2900, 2925, 2950, 2975, 3000, 3025, 
                         3050, 3075, 3100, 3150, 3175, 3200, 3225, 3250, 3275, 
                         3350, 3375, 3400, 3450, 3500, 3700, 3750, 3875, 4425, 
                         4450, 4525, 4625, 4750, 4850, 4875, 4925, 4975, 5000, 
                         5050, 5150, 5200, 5225, 5325, 5500, 5550, 5650, 5700, 
                         6100, 6325, 6350, 6450, 6575, 6750, 6800, 6975, 7375, 
                         7600]                                                      #Manually Selected bad frames for A04.  
    
    elif dataname=='3rd Good Dataset - A05':  
        chosen_outliers=[50, 150, 200, 275, 325, 550, 675, 700, 750, 775, 875, 
                         975, 1000, 1050, 1075, 1225, 1300, 1350, 1600, 1650, 
                         1675, 1725, 1800, 1825, 2050, 2075, 2200, 2275, 2325, 
                         2375, 2575, 2600, 2725, 2750, 2775, 2875, 2900, 2950, 
                         3000, 3050, 3125, 3150, 3250, 3275, 3400, 3450, 3475, 
                         3500, 3525, 3600, 3625, 3675, 3700, 3825, 3850, 3950, 
                         4000, 4025, 4325, 4375, 4450, 4475, 4500, 4975, 5050, 
                         5200, 5225, 5250, 5375, 5575, 5650, 5850, 5875, 5900, 
                         5950, 6025, 6150, 6350, 6575, 4625, 6750, 6775, 6850, 
                         7025, 7050, 7125, 7250, 7375, 7450, 7475, 7700, 7750, 
                         8075, 8200]                                                #Manually Selected bad frames for A05. 
        
    elif dataname=='3rd Good Dataset - A06':  
        chosen_outliers=[0, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 
                         300, 325, 600, 675, 700, 725, 750, 775, 800, 825, 
                         850, 875, 900, 1075, 1100, 1275, 1325, 1775, 2000, 
                         2250, 2825, 3050, 3100, 3300, 3350, 3550, 4125, 5225, 
                         5675, 5700, 5925, 6725, 6775, 6850, 7675, 7725]            #Manually Selected bad frames for A06. 
    
    elif dataname=='3rd Good Dataset - A07':  
        chosen_outliers=[425, 450, 650, 825, 1150, 1250, 1325, 1425, 1450, 1525, 
                         1675, 1825, 2575, 2625, 2650, 2975, 3925, 4150, 4375, 
                         4425, 4475, 4700, 4800, 5000, 5050, 5250, 5275, 5400, 
                         5450, 5600, 5625, 5650, 5750, 5800, 6025, 6075, 6475, 
                         6550, 7025, 7300, 7600, 7750, 7850]                         #Manually Selected bad frames for A07
    
    elif dataname=='3rd Good Dataset - A08':  
        chosen_outliers=[50, 75, 100, 225, 275, 300, 325, 350, 375, 450, 525, 550, 
                         600, 625, 675, 700, 725, 775, 800, 825, 850, 875, 900, 
                         925, 1000, 1025, 1050, 1100, 1125, 1175, 1200, 1225,
                         1250, 1275, 1300, 1350, 1400, 1450, 1500, 1550, 1575, 
                         1600, 1650, 1675, 1700, 1725, 1750, 1800, 1825, 1850, 
                         1875, 1900, 1975, 2425, 2450, 2475, 2600, 2925, 2975, 
                         3800, 4125, 4250, 4650, 4675, 4700, 5050, 5075, 5125, 
                         5425]                                    #Manually Selected bad frames for A08. 
    
    elif dataname=='3rd Good Dataset - A09':  
        chosen_outliers=[50, 1875, 1900, 2450, 2450, 2475, 2500, 2525, 2550,
                         2575, 2600, 2625, 2650, 2675, 2700, 2725, 2750, 2775,
                         2800, 2825, 2850, 2875, 2900, 2925, 2950, 2975, 3000,
                         3025, 3050, 3075, 3100, 3125, 3225, 3250, 3475, 3500, 
                         3525, 3725, 4300, 5775, 5800, 5825, 5850, 5875, 5900, 
                         5925, 5950, 5975, 6200, 6225, 6250, 6275, 6300, 6325,
                         6350, 6375, 6400, 6425, 6450, 6475, 6500, 6525, 6550,
                         6575, 6600, 6625, 6650, 6675, 6700, 6725, 6750, 6775,
                         6800, 6825, 6850, 6875, 6900, 6925, 6950, 6975, 7000, 
                         7025, 7050]                                    #Manually Selected bad frames for A09.
    
    elif dataname=='3rd Good Dataset - A10':  
        chosen_outliers=[0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275,
                         300, 325, 350, 375, 400, 425, 450, 475, 500, 525, 550,
                         575, 600, 625, 650, 675, 700, 725, 750, 775, 800, 825,
                         850, 875, 900, 925, 950, 975, 1000, 1025, 1050, 1075,
                         1100, 1125, 1150, 1175, 1200, 1225, 1250, 1275, 1300,
                         1525, 2150, 2200, 2225, 2250, 2275, 2425, 2475, 2500, 
                         2575, 2600, 2625, 2650, 2675, 2700, 2725, 2750, 2775,
                         2875, 2975, 3075, 3200, 3500, 3650, 3700, 4100, 4400, 
                         4425, 4475, 4500, 4575, 4600, 4650, 4675, 4750, 4800,
                         4975, 5100, 5150, 5225, 5550, 5775]            #Manually Selected bad frames for A10.
        

    
    elif dataname==0:      #Used if removing frames from a pre-defined 'chosen_outliers' variable (i.e. frames already identified as anomalous from a different method
                    #and saved under the 'chosen outliers' variable name)
        pass
    else:
        print('Incorrect value put into "remove_chosen_outliers(dataname)". Dataname = '+ dataname)
     
    old_b=-1
    for k in chosen_outliers: 
        bad_row.append(data_dict_of_roughness_parameters["Distance (um)"].index(k))     #finds the row of the bad data
    for a, b in enumerate(bad_row):
        if a==0:                                    #Removes those bad data rows from all of the data
            nth_term=b
        elif b==old_b:
            pass
        else: 
            nth_term=b-a
        
        Distance_outliers.append(data_dict_of_roughness_parameters["Distance (um)"].pop(nth_term))
        Avg_height_outliers.append(data_dict_of_roughness_parameters["Average Height (nm)"].pop(nth_term))
        Ra_outliers.append(data_dict_of_roughness_parameters["AM, Ra/Sa (nm)"].pop(nth_term))
        Rms_outliers.append(data_dict_of_roughness_parameters["RMS, Rq/Sq (nm)"].pop(nth_term))
        Skew_outliers.append(data_dict_of_roughness_parameters["Rsk/Ssk"].pop(nth_term))
        Kurtosis_outliers.append(data_dict_of_roughness_parameters["Rku/Sku"].pop(nth_term))
        Zmax_outliers.append(data_dict_of_roughness_parameters["Rz/Sz (nm)"].pop(nth_term))
        
        old_b=b  
    
    clean_data=data_dict_of_roughness_parameters
    
    outlier_data={
            "Distance outliers(um)": Distance_outliers,             
            "Average Height outliers(nm)": Avg_height_outliers,     
            "AM, Ra/Sa outliers(nm)": Ra_outliers,                  
            "RMS, Rq/Sq outliers(nm)": Rms_outliers,                
            "Rsk/Ssk outliers": Skew_outliers,                      
            "Rku/Sku outliers": Kurtosis_outliers,                  
            "Rz/Sz outliers(nm)": Zmax_outliers                     
            }
    
    return clean_data, outlier_data

def read_all_data(All_roots, All_txtfilepaths, dataset_names, tablefilename):
    """"A function that reads all the data from all samples and all datasets, returning:
         == a dictionary containing all the data for the various roughness parameters under their dataset and sample names in the format: 
            ' all_data["Fibre X"]["Sample name"]["Roughness Parameter Name"] '
         == a dictionary that has all manually selected anomalous frames removed as per the frames listed in 'remove_chosen_outliers()'
         == and a dictionary of all the frames that have been removed.
    """
    global readstats, all_data, all_clean_data, all_outlier_data
    all_data={}
    all_clean_data={}
    all_outlier_data={}
    
    """Example inputs are:
    roots = [root_1, root_2, root_3, root_4]
    txtfilepaths=[txtfilepaths_1, txtfilepaths_2, txtfilepaths_3, txtfilepaths_4]
    dataset_names = ["Fibre 1", "Fibre 1 Repeat", "Fibre 2", "Fibre 2 - Used Tips"]
    tablefilename='Line_data_2.txt'
    """
    for a, txtfilepath in enumerate(All_txtfilepaths):
        sampletitles=findtitles(All_roots[a])   #Reading the name of the samples
        # print(sampletitles)
        for b, path in enumerate(txtfilepath):  #Iterating through the txtfilepaths for each sample
            ##Reading all of the data and saving it to the all_data dictionary
            readstats=False
            os.chdir(path)                      #Changes to the directory of each sample
            data = readstats_readlastframe(tablefilename)   #Reads all the roughness parameters for that sample
            append_value(all_data, dataset_names[a], {sampletitles[b]:data})        #Adds that sample roughness data to the 'all_data' dictionary with the associated sample name

            ##Reading all of the data, removing anomalous frames and saving it to the clean_data and outlier_data dictionary
            readstats=False
            currdata = readstats_readlastframe(tablefilename)   #Reads all the roughness parameters for that sample
            dataname = dataset_names[a] + " - " + sampletitles[b]
            # print(dataname+'  :  '+path)
            clean_data, outlier_data = remove_chosen_outliers(currdata, dataname, tablefilename) #Removes anomalous data from the data just read
            append_value(all_clean_data, dataset_names[a], {sampletitles[b]:clean_data})
            append_value(all_outlier_data, dataset_names[a], {sampletitles[b]:outlier_data})
            
    return all_data, all_clean_data, all_outlier_data


def histogram(roughness_data, column_number, method_number, title_of_sample, tagline, no_of_measurements, other_data=[False,False], first_or_last='first', savefile=False):
    """
    Plots the histogram for roughness data and its corresponding .txtfile column number (1 - Average Height, 2 - Ra/Sa and sample title etc.), 
    the method number (1=Area roughness, 2=Line roughness), the title of the sample, the tagline to add onto the end of the graph's file name
    and the number of n measurements made.
    """
    
    if method_number!=1 and method_number!=2:
        print("Please input either 1 or 2 for the method_number input in the histogram plotting function")
        return
    elif type(title_of_sample)!=str or type(tagline)!=str:
        print("Please input a string as the 4th and 5th inputs to histogram(). Sample title and file tagline respectively.")
        return
    elif type(no_of_measurements)!=int:
        print("Please input an integer as the number of measurements (6th input) to histogram()")
        return
    else:
        pass
    
    if column_number==1:
        b='Average Height'
        xscale=1e9
        xlabel='Average Height(nm)'
        units='nm'
    
    elif column_number==2 and method_number==1:
        b='Sa'
        xscale=1e9
        xlabel='Sa(nm)'
        units='nm'  
    elif column_number==2 and method_number==2:
        b='Ra'
        xscale=1e9
        xlabel='Ra(nm)'
        units='nm'
        
    elif column_number==3 and method_number==1:
        b='Sq'
        xscale=1e9
        xlabel='Sq(nm)'
        units='nm'      
    elif column_number==3 and method_number==2:
        b='Rq'
        xscale=1e9
        xlabel='Rq(nm)'
        units='nm'
   
    elif column_number==4:
        b='Skew'
        xscale=1
        xlabel='Skew'
        units=''
  
    elif column_number==5:
        b='Kurtosis'
        xscale=1
        xlabel='Kurtosis'
        units=''
    
    elif column_number==6 and method_number==1:
        b='Sz'
        xscale=1e9
        xlabel='Sz(nm)'
        units='nm'
    
    elif column_number==6 and method_number==2:
        b='Rz'
        xscale=1e9
        xlabel='Rz(nm)'
        units='nm'
    else:
        pass

    if other_data[0]==False:
        roughness_data=numpy.array(roughness_data)*xscale
    mu = numpy.mean(roughness_data)
    sigma = numpy.std(roughness_data)
    med = numpy.median(roughness_data)
    
    n, bins, patches = plt.hist(roughness_data, bins='auto', range=None, density=True,\
                             facecolor='blue', alpha=0.6, rwidth=0.8)
    # pdf=norm.pdf(bins, mu, sigma)
    # plt.plot(bins,pdf,'r--', label='PDF')

    sns.distplot(roughness_data, hist=False, kde=True, 
             bins=bins, color = 'purple',
             kde_kws={'linewidth': 2, 'alpha':0.8},label='KDE')

    if other_data[0]==False:
        plt.xlabel(xlabel)
        plt.axvline(x=med,c='red', alpha=0.7, label='Median: %.3g'%(med)+units)
        plt.axvline(x=mu,c='orange', alpha=0.8, label='Mean: %.3g'%(mu)+units) 
    else:
        if first_or_last=='Overlap':
            plt.xlabel('Overlap Amount')
            plt.xlim(left=0)
            plt.axvline(x=other_data[1],c='red', alpha=0.7, label='Mean: %.3g'%(other_data[1])+units)
            plt.axvline(x=other_data[2],c='orange', alpha=0.8, label='Median: %.3g'%(other_data[2])+units)
        else:
            if first_or_last=='last':
                plt.xlabel('Last Convergence Points')
            elif first_or_last=='first':
                plt.xlabel('First Convergence Points')

            plt.xlim(left=0)
            plt.axvline(x=other_data[1],c='red', alpha=0.7, label='95th Percentile')
            plt.axvline(x=other_data[2],c='orange', alpha=0.8, label='97th Percentile')
        
    plt.ylabel('Probability Density')
    plt.suptitle(title_of_sample +' - '+b, fontsize=14)
    plt.title(u'\u03C3 = %.3g' %(sigma)+units+', n = %i'%(no_of_measurements))
    plt.legend()
    plt.pause(0.5)
    if savefile==True:
        plt.savefig(title_of_sample +' - Histogram of '+ b + tagline +'.png')
    # print(title_of_sample +' - Histogram of '+ b + tagline +'.png')
    plt.close()

def boxplot(roughness_data, column_number, method_number, title_of_sample, tagline, no_of_measurements, savefile=False):
    """Plots a boxplot for 'g' data and 'w' sample title"""
    if method_number!=1 and method_number!=2:
        print("Please input either 1 or 2 for the method_number input in the histogram plotting function")
        return
    elif type(title_of_sample)!=str or type(tagline)!=str:
        print("Please input a string as the 4th and 5th inputs to histogram(). Sample title and file tagline respectively.")
        return
    elif type(no_of_measurements)!=int:
        print("Please input an integer as the number of measurements (6th input) to histogram()")
        return
    else:
        pass
    
    if column_number==1:
        b='Average Height'
        xscale=1e9
        xlabel='Average Height(nm)'
        units='nm'
    
    elif column_number==2 and method_number==1:
        b='Sa'
        xscale=1e9
        xlabel='Sa(nm)'
        units='nm'  
    elif column_number==2 and method_number==2:
        b='Ra'
        xscale=1e9
        xlabel='Ra(nm)'
        units='nm'
        
    elif column_number==3 and method_number==1:
        b='Sq'
        xscale=1e9
        xlabel='Sq(nm)'
        units='nm'      
    elif column_number==3 and method_number==2:
        b='Rq'
        xscale=1e9
        xlabel='Rq(nm)'
        units='nm'
   
    elif column_number==4:
        b='Skew'
        xscale=1
        xlabel='Skew'
        units=''
  
    elif column_number==5:
        b='Kurtosis'
        xscale=1
        xlabel='Kurtosis'
        units=''
    
    elif column_number==6 and method_number==1:
        b='Sz'
        xscale=1e9
        xlabel='Sz(nm)'
        units='nm'
    
    elif column_number==6 and method_number==2:
        b='Rz'
        xscale=1e9
        xlabel='Rz(nm)'
        units='nm'
    else:
        pass
    
    roughness_data=numpy.array(roughness_data)*xscale
    mu = numpy.mean(roughness_data)
    sigma = numpy.std(roughness_data)

    
    Q1=numpy.quantile(roughness_data, 0.25)
    Q3=numpy.quantile(roughness_data, 0.75)
    IQR = Q3-Q1
    Low_cutoff=Q1-1.5*IQR
    High_cutoff=Q3+1.5*IQR
    low_stdev=mu-(2*sigma)
    upper_stdev=mu+(2*sigma)
    low_stdev3=mu-(3*sigma)
    upper_stdev3=mu+(3*sigma)
    
    sns.boxplot(x=roughness_data)

    plt.axvline(x=Low_cutoff,c='red', alpha=0.3, label='Lower 1.5 IQR < %.3g'%(Low_cutoff)+units)
    plt.axvline(x=High_cutoff,c='red', alpha=0.3, label='Upper 1.5 IQR > %.3g'%(High_cutoff)+units) 
    plt.axvline(x=mu,c='blue', alpha=0.6, label=u'\u03bc = %.3g'%(mu) + units) 
    plt.axvline(x=low_stdev,c='blue', alpha=0.3, label=u'Lower 2\u03C3 < %.3g'%(low_stdev)+units)
    plt.axvline(x=upper_stdev,c='blue', alpha=0.3, label=u'Upper 2\u03C3 > %.3g'%(upper_stdev)+units)
    plt.axvline(x=low_stdev3,c='purple', alpha=0.3, label=u'Lower 3\u03C3 < %.3g'%(low_stdev3)+units)
    plt.axvline(x=upper_stdev3,c='purple', alpha=0.3, label=u'Upper 3\u03C3 > %.3g'%(upper_stdev3)+units)

    plt.xlabel(xlabel)
    plt.suptitle(title_of_sample+' - '+b, fontsize=14)
    plt.title('IQR =%.3g' %(IQR)+units)
    plt.legend()
    plt.pause(0.5)
    if savefile==True:
        plt.savefig(title_of_sample +' - Boxplot of '+ b +'.png')     #Change title accordingly
    plt.close()
    
def samplecomparison_histogram(roughness_data, column_number, method_number, sample_names, fibre_section_title, tagline, savefile=False):
    """Plots a graph for each roughness parameter, comparing the distributions of the various samples.
       Plots the roughness data and its corresponding .txtfile column number (1 - Average Height, 2 - Ra/Sa and sample title etc.), 
       the method number (1=Area roughness, 2=Line roughness), a lsit of sample names,
       fibre_section_title to add to the title, the tagline to add onto the end of the graph's file name.
       
       The roughness data needs to be a dictionary containing the data for the parameter to be plotted for the various samples (in order) in the format of:
           'roughness_data["Sample name"]["Roughness Parameter Name"]'
       
        """
    ##Check if the some of the inputs are in the correct format/order
    if method_number!=1 and method_number!=2:
        print("Please input either 1 or 2 for the method_number input in the histogram plotting function")
        return
    elif type(sample_names)!=list or type(tagline)!=str:
        print("Please input a list as the 4th and string as the 5th inputs to samplecomparison_histogram(). Sample titles and file tagline respectively.")
        return
    else:
        pass
    
    colour_list=['blue','orange','green','red','purple','brown','pink','gray','olive','cyan']
    
    if column_number==1:
        b='Average Height'
        xscale=1e9
        xlabel='Average Height(nm)'
        units='nm'
        xmin=0
        xlimit=500
    
    elif column_number==2 and method_number==1:
        b='Sa'
        xscale=1e9
        xlabel='Sa(nm)'
        units='nm'   
        xmin=0
        xlimit=200
    elif column_number==2 and method_number==2:
        b='Ra'
        xscale=1e9
        xlabel='Ra(nm)'
        units='nm'
        
    elif column_number==3 and method_number==1:
        b='Sq'
        xscale=1e9
        xlabel='Sq(nm)'
        units='nm' 
        xmin=0
        xlimit=200
    elif column_number==3 and method_number==2:
        b='Rq'
        xscale=1e9
        xlabel='Rq(nm)'
        units='nm'
   
    elif column_number==4:
        b='Skew'
        xscale=1
        xlabel='Skew'
        units=''
        xmin=-3
        xlimit=3
  
    elif column_number==5:
        b='Kurtosis'
        xscale=1
        xlabel='Kurtosis'
        units=''
        xmin=-4
        xlimit=5
    
    elif column_number==6 and method_number==1:
        b='Sz'
        xscale=1e9
        xlabel='Sz(nm)'
        units='nm'
        xmin=0
        xlimit=700
    elif column_number==6 and method_number==2:
        b='Rz'
        xscale=1e9
        xlabel='Rz(nm)'
        units='nm'
    else:
        pass
    
    for i, samplename in enumerate(sample_names):     
        sampledata=roughness_data[samplename]
        
        #Correct the scale of the data
        x= [j * xscale for j in sampledata]
        #Calculate its mean
        mu = numpy.mean(x)    
        #Calculate the standard deviation
        sigma = numpy.std(x)  
        #Calculate the median
        med = numpy.median(x) 
        #Create the  label for the histogram legend
        label= samplename   #Creating the curve labels
    
        histo=False
        rugo=False
        #Add the histogram to the plot with a different colour set by colours_list
        sns.distplot([x], hist=histo, kde=True, rug=rugo, bins='auto', color = colour_list[i], 
                 hist_kws={'edgecolor':'black'}, kde_kws={'linewidth': 1.5, "alpha": 0.5}, label=label)
        
    plt.xlabel(xlabel)
    plt.ylabel('Probability Density')
    plt.suptitle(b, fontsize=16)
    plt.xlim(left=xmin, right=xlimit)
    plt.title('KDE plots for All Fibre Samples - '+fibre_section_title)
    plt.legend(prop={"size":9})
    # plt.pause(1)
    if savefile==True:
        plt.savefig(fibre_section_title+' - Histogram of '+ b + tagline+'.png')          #Change the title accordingly
    plt.close()
    
def mean_and_uncertainty_plot(roughness_data, column_number, method_number, sample_names, fibre_section_title, tagline, savefile=False, only_errors=False):
    """Plots the mean roughnesses for each of the fibres in one graph with uncertainty bars to check if there is overlap"""
    # error_range = input("Enter what confidence for interval (68%, 95% or 99.7%): \n")   
    # if error_range!=68 and error_range!=95 and error_range!=99.7:
    #     print('Please enter a confidence interval of "68", "95" or "99.7"')
    
    if method_number!=1 and method_number!=2 and method_number!=0:
        print("Please input either 1 or 2 for the method_number input in the histogram plotting function. or 0 if using the function generally")
        return
    elif type(sample_names)!=list or type(tagline)!=str:
        print("Please input a list as the 4th and string as the 5th inputs to samplecomparison_histogram(). Sample titles and file tagline respectively.")
        return
    else:
        pass
    roughness_means=[]
    standard_uncertainties=[]
    n_sites_list=[]
    xscale=1
    
    if column_number==1:
        b='Average Height'
        xscale=1e9
        ylabel='Average Height(nm)'
        units='nm'
    
    elif column_number==2 and method_number==1:
        b='Sa'
        xscale=1e9
        ylabel='Sa(nm)'
        units='nm'   

    elif column_number==2 and method_number==2:
        b='Ra'
        xscale=1e9
        ylabel='Ra(nm)'
        units='nm'
        
    elif column_number==3 and method_number==1:
        b='Sq'
        xscale=1e9
        ylabel='Sq(nm)'
        units='nm' 

    elif column_number==3 and method_number==2:
        b='Rq'
        xscale=1e9
        ylabel='Rq(nm)'
        units='nm'
   
    elif column_number==4:
        b='Skew'
        ylabel='Skew'
        units=''
  
    elif column_number==5:
        b='Kurtosis'
        ylabel='Kurtosis'
        units=''
    
    elif column_number==6 and method_number==1:
        b='Sz'
        xscale=1e9
        ylabel='Sz(nm)'
        units='nm'

    elif column_number==6 and method_number==2:
        b='Rz'
        xscale=1e9
        ylabel='Rz(nm)'
        units='nm'
    
    elif column_number==0 and method_number==0:         #The inputs that choose the general labels
        ylabel='Bend Strain (%)'
    else:
        pass
    
    for i, samplename in enumerate(sample_names):     
        sampledata=roughness_data[samplename]
        
        #Correct the scale of the data
        x= [j * xscale for j in sampledata]
        #Count the number of datapoints
        n_sites=int(len(x))
        n_sites_list.append(n_sites)
        #Calculate its mean
        mu = numpy.mean(x)    
        #Calculate the standard deviation
        sigma = numpy.std(x)  
        #Calculate the median
        med = numpy.median(x) 
        #Calculate u by first doing the sqrt of n
        sqrtn=math.sqrt(n_sites)
        u = sigma/sqrtn
        roughness_means.append(mu)
        standard_uncertainties.append(u)


    x_pos = numpy.arange(len(sample_names))
    
    if n_sites<30:
        #Uses the t-score if the number of measurements is too small
        t_score_95=[]
        t_score_99_7=[]
        for i, sample_n_site in enumerate(n_sites_list):     
            t_score_95.append(t.ppf(1-0.025, sample_n_site-1))
            t_score_99_7.append(t.ppf(1-0.0015, sample_n_site-1))
        expanded_errors_95percent=[a*b for a,b in zip(t_score_95,standard_uncertainties)]
        expanded_errors_99_7percent=[a*b for a,b in zip(t_score_99_7,standard_uncertainties)]
        print('t-score used')
    else:
        expanded_errors_95percent=[j * 1.96 for j in standard_uncertainties]
        expanded_errors_99_7percent=[j * 3 for j in standard_uncertainties]

    error=expanded_errors_95percent
    # print(error)
    if only_errors==True:
        return(roughness_means, error, standard_uncertainties)
    fig, ax = plt.subplots()
    ax.errorbar(x_pos, roughness_means, yerr=error, fmt='x', alpha=1, ecolor='black', capsize=10)
    ax.set_ylabel(ylabel)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(sample_names)
    ax.set_xlabel('Samples')
    plt.xlim(-1,len(sample_names))
    if method_number==0:
        ax.set_title(fibre_section_title)
    else:
        ax.set_title('Comparison of Samples - ' + fibre_section_title +' - '+b)
    ax.yaxis.grid(True)
    
    # Save the figure and show
    plt.tight_layout()
    plt.show()
    # plt.pause(5)
    if savefile==True and method_number==0:
        plt.savefig(fibre_section_title +' - '+tagline+'.png')
    elif savefile==True:
        plt.savefig('Comparison of Sample Means - '+ fibre_section_title +' - '+b+' ' +tagline+'.png')
    plt.close()  
    
def bootstrap_confidence_intervals(data, column_number, estimator, percentiles, runs=10000):
    """    
    1 Take a large number of samples of the same size as our original dataset, by sampling with replacement
    2 For each sample, calculate the parameter of interest (eg. the mean)
    3 Calculate the relevant percentile from the distribution of the parameter of interest"""
    if column_number==1:
        xscale=1e9
    elif column_number==2:
        xscale=1e9    
    elif column_number==3:
        xscale=1e9
    elif column_number==4:
        xscale=1
    elif column_number==5:
        xscale=1
    elif column_number==6:
        xscale=1e9
    else:
        pass
    # data= [j * xscale for j in data]
    replicates = numpy.empty(runs)
    for i in range(runs):
        replicates[i] = estimator(numpy.random.choice(data, len(data), replace=True))
    est = numpy.mean(replicates)
    ci = numpy.percentile(numpy.sort(replicates), percentiles)
    return (est, ci)

def random_sampling_within_confidence_intervals(x_data_means, x_uncertainties, y_data_means, y_uncertainties, correlationfunc, percentiles, number_of_samples=10, runs=10000):
    """A function used to carry across the uncertainty in means to the correlation function of interest through random selection according to the uncertainty distributions.    
    1 Randomly select values of the means according to a normal distribution where the st unc of mean is the st dev.
    2 For each set of randomly selected values, calculate the correlation parameter of interest (eg. the Pearsons Correlation Coefficient)
    3 Calculate the relevant percentile from the distribution of the parameter of interest"""
    replicates = numpy.empty(runs)
    x_means = numpy.empty((runs, number_of_samples))
    y_means = numpy.empty((runs, number_of_samples))
    
    #For each x sample mean...
    for num1, xmu in enumerate(x_data_means):
        #Find its associated uncertainty
        sigma = x_uncertainties[num1]
        #And then
        for i in range(runs):
            #Randomly select a value in the normal distribution for each x_data_mean
            x_means[i][num1]=numpy.random.normal(xmu, sigma)
    
    #For each y sample mean...
    for num2, ymu in enumerate(y_data_means):
        #Find its associated uncertainty
        sigma = y_uncertainties[num2]
        #And then
        for i in range(runs):
            #Randomly select a value in the normal distribution for each y_data_mean
            y_means[i][num2]=numpy.random.normal(ymu, sigma)
    
    #Now read each iteration (i.e row in the arrays) 
    for j, xrow in enumerate(x_means):
        x_data=xrow
        y_data=y_means[j]
        #...and see if they have a correlation (correlationfunc)
        replicates[j], _ =correlationfunc(x_data, y_data)
    
    #'est' is estimated mean of correlation function
    est = numpy.mean(replicates)
    #'ci' is the confidence intervals for the correlation function according to the input values in 'percentiles'
    ci = numpy.percentile(numpy.sort(replicates), percentiles)
    return (est, ci)


def uncertainty_plotter(roughness_data, column_number, method_number, sample_names, fibre_section_title, tagline, savefile=False, only_errors=False):
    """A mean and uncertainty plot that compares between the samples for each roughness parameter. 
    It uses the bootstrap confidence interval function to identify the CI to plot."""
    if method_number!=1 and method_number!=2 and method_number!=0:
        print("Please input either 1 or 2 for the method_number input in the histogram plotting function. or 0 if using the function generally, e.g. bend strain results")
        return
    elif type(sample_names)!=list or type(tagline)!=str:
        print("Please input a list as the 4th and string as the 5th inputs to samplecomparison_histogram(). Sample titles and file tagline respectively.")
        return
    else:
        pass
    
    roughness_means=[]
    CI_95=[]
    CI_99_7=[]
    xscale=1
    
    if column_number==1:
        b='Average Height'
        xscale=1e9
        ylabel='Average Height(nm)'
        units='nm'
    
    elif column_number==2 and method_number==1:
        b='Sa'
        xscale=1e9
        ylabel='Sa(nm)'
        units='nm'   

    elif column_number==2 and method_number==2:
        b='Ra'
        xscale=1e9
        ylabel='Ra(nm)'
        units='nm'
        
    elif column_number==3 and method_number==1:
        b='Sq'
        xscale=1e9
        ylabel='Sq(nm)'
        units='nm' 

    elif column_number==3 and method_number==2:
        b='Rq'
        xscale=1e9
        ylabel='Rq(nm)'
        units='nm'
   
    elif column_number==4:
        b='Skew'
        ylabel='Skew'
        units=''
  
    elif column_number==5:
        b='Kurtosis'
        ylabel='Kurtosis'
        units=''
    
    elif column_number==6 and method_number==1:
        b='Sz'
        xscale=1e9
        ylabel='Sz(nm)'
        units='nm'

    elif column_number==6 and method_number==2:
        b='Rz'
        xscale=1e9
        ylabel='Rz(nm)'
        units='nm'
    
    elif column_number==0 and method_number==0:         #The inputs that choose the label for the bendstrain graph
        ylabel='Bend Strain (%)'
    else:
        pass
    
    for i, samplename in enumerate(sample_names):     
        sampledata=roughness_data[samplename]
        
        #Correct the scale of the data
        x= [j * xscale for j in sampledata]

        mu, conf_int_95percent = bootstrap_confidence_intervals(x, column_number, numpy.mean, [2.5, 97.5])
        _, conf_int_99_7percent = bootstrap_confidence_intervals(x, column_number, numpy.mean, [0.15, 99.85])
        roughness_means.append(mu)
        CI_95.append([mu- conf_int_95percent[0], conf_int_95percent[1]-mu])
        CI_99_7.append([mu- conf_int_99_7percent[0], conf_int_99_7percent[1]-mu])

    # CI_95_array=numpy.array(CI_95)
    
    x_pos = numpy.arange(len(sample_names))
    error=numpy.transpose(CI_95)
    # print(error)
    
    if only_errors==True:
        return(roughness_means, error)
    
    fig, ax = plt.subplots()
    ax.errorbar(x_pos, roughness_means, yerr=error, fmt='x', alpha=1, ecolor='black', capsize=10)
    ax.set_ylabel(ylabel, weight='bold')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(sample_names)
    ax.set_xlabel('Samples', weight='bold')
    plt.xlim(-1,len(sample_names))
    if method_number==0:
        ax.set_title(fibre_section_title)
    else:
        ax.set_title('Comparison of Samples - ' + fibre_section_title +' - '+b)
    ax.yaxis.grid(True)
    
    # Save the figure and show
    plt.tight_layout()
    plt.show()
    # plt.pause(5)
    if savefile==True and method_number==0:
        plt.savefig(fibre_section_title +' - '+tagline+'.png')
    elif savefile==True:
        plt.savefig('Comparison of Sample Means - '+ fibre_section_title +' - '+b+' ' +tagline+'.png')
    plt.close()

def combined_uncertainty_plotter(roughness_data_1, roughness_data_2, column_number, method_number, sample_names, two_dataset_labels, dataset_comparison_title, tagline, savefile=False):
    """A mean and uncertainty plot that compares between the datasets for each sample. 
    It uses the bootstrap confidence interval function to identify the CI to plot."""
    if method_number!=1 and method_number!=2 and method_number!=0:
        print("Please input either 1 or 2 for the method_number input in the histogram plotting function. or 0 if comparing between the two")
        return
    elif type(sample_names)!=list or type(tagline)!=str:
        print("Please input a list as the 4th and string as the 5th inputs to samplecomparison_histogram(). Sample titles and file tagline respectively.")
        return
    else:
        pass
    
    roughness_means_1=[]
    roughness_means_2=[]
    CI_95_dataset1=[]
    CI_99_7_dataset1=[]
    CI_95_dataset2=[]
    CI_99_7_dataset2=[]
    xscale=1
    
    if column_number==1:
        b='Average Height'
        xscale=1e9
        ylabel='Average Height (nm)'
        units='nm'
    
    elif column_number==2 and method_number==1:
        b='Sa'
        xscale=1e9
        ylabel='Sa (nm)'
        units='nm'   

    elif column_number==2 and method_number==2:
        b='Ra'
        xscale=1e9
        ylabel='Ra (nm)'
        units='nm'
        
    elif column_number==3 and method_number==1:
        b='Sq'
        xscale=1e9
        ylabel='Sq (nm)'
        units='nm' 

    elif column_number==3 and method_number==2:
        b='Rq'
        xscale=1e9
        ylabel='Rq (nm)'
        units='nm'
   
    elif column_number==4:
        b='Skew'
        ylabel='Skew'
        units=''
  
    elif column_number==5:
        b='Kurtosis'
        ylabel='Kurtosis'
        units=''
    
    elif column_number==6 and method_number==1:
        b='Sz'
        xscale=1e9
        ylabel='Sz (nm)'
        units='nm'

    elif column_number==6 and method_number==2:
        b='Rz'
        xscale=1e9
        ylabel='Rz (nm)'
        units='nm'
    
    elif column_number==0 and method_number==0:         #The inputs that choose the label for the bendstrain graph
        ylabel='Bend Strain (%)'
    else:
        pass
    
    for i, samplename in enumerate(sample_names):     
        sampledata_1=roughness_data_1[samplename]
        sampledata_2=roughness_data_2[samplename]
        
        #Correct the scale of the data
        x1= [j * xscale for j in sampledata_1]
        x2= [j * xscale for j in sampledata_2]

        mu1, conf_int_95percent1 = bootstrap_confidence_intervals(x1, column_number, numpy.mean, [2.5, 97.5], runs=60000)
        # _, conf_int_99_7percent1 = bootstrap_confidence_intervals(x1, column_number, numpy.mean, [0.15, 99.85])
        roughness_means_1.append(mu1)
        CI_95_dataset1.append([mu1- conf_int_95percent1[0], conf_int_95percent1[1]-mu1])
        # CI_99_7_dataset1.append([mu1- conf_int_99_7percent1[0], conf_int_99_7percent1[1]-mu1])
        
        mu2, conf_int_95percent2 = bootstrap_confidence_intervals(x2, column_number, numpy.mean, [2.5, 97.5], runs=60000)
        # _, conf_int_99_7percent2 = bootstrap_confidence_intervals(x2, column_number, numpy.mean, [0.15, 99.85])
        roughness_means_2.append(mu2)
        CI_95_dataset2.append([mu2- conf_int_95percent2[0], conf_int_95percent2[1]-mu2])
        # CI_99_7_dataset2.append([mu2- conf_int_99_7percent2[0], conf_int_99_7percent2[1]-mu2])

    # CI_95_array1=numpy.array(CI_95_dataset1)
    # CI_95_array2=numpy.array(CI_95_dataset2)
    
    x_pos = numpy.arange(len(sample_names))
    error1=numpy.transpose(CI_95_dataset1)
    error2=numpy.transpose(CI_95_dataset2)
    
    #=============================Plotting the graph ===============================================
    fig, ax = plt.subplots()
    ax.errorbar(x_pos, roughness_means_1, yerr=error1, fmt='x', alpha=0.6, color='blue', ecolor='blue', capsize=10)
    ax.errorbar(x_pos, roughness_means_2, yerr=error2, fmt='x', alpha=0.6, color='red', ecolor='red', capsize=10)
    ax.set_ylabel(ylabel, weight='bold')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(sample_names)
    ax.set_xlabel('Samples', weight='bold')
    plt.xlim(-1,len(sample_names))
    if method_number==0:
        ax.set_title(dataset_comparison_title)
    else:
        # ax.set_title('Comparison of Sample Means Between Two Repeat Datasets - '+b)
        pass
    ax.yaxis.grid(True)
    
    # Save the figure and show
    ##Creating the legend
    handles, labels = ax.get_legend_handles_labels()
    # manually define new patches for the different bands 
    u1_band = mpatches.Patch(color='b', alpha=0.4, label=two_dataset_labels[0])
    u2_band = mpatches.Patch(color='r', alpha=0.4, label=two_dataset_labels[1])
    handles.append(u1_band) 
    handles.append(u2_band) 
    plt.legend(handles=handles, prop={"size":8})
    plt.tight_layout()
    plt.show()
    plt.pause(0.1)
    if savefile==True and method_number==0:
        plt.savefig(dataset_comparison_title +' - '+tagline+'.png')
    elif savefile==True:
        plt.savefig('Comparison of Sample Means - '+ dataset_comparison_title +' - '+b+' ' +tagline+'.png')
    plt.close()

def three_combined_uncertainty_plotter(roughness_data_1, roughness_data_2, roughness_data_3, column_number, method_number, sample_names, three_dataset_labels, dataset_comparison_title, tagline, savefile=False):
    """A mean and uncertainty plot that compares between the datasets for each sample. 
    It uses the bootstrap confidence interval function to identify the CI to plot."""
    if method_number!=1 and method_number!=2 and method_number!=0:
        print("Please input either 1 or 2 for the method_number input in the histogram plotting function. or 0 if comparing between the two")
        return
    elif type(sample_names)!=list or type(tagline)!=str:
        print("Please input a list as the 4th and string as the 5th inputs to samplecomparison_histogram(). Sample titles and file tagline respectively.")
        return
    else:
        pass
    
    roughness_means_1=[]
    roughness_means_2=[]
    roughness_means_3=[]
    CI_95_dataset1=[]
    CI_99_7_dataset1=[]
    CI_95_dataset2=[]
    CI_99_7_dataset2=[]
    CI_95_dataset3=[]
    CI_99_7_dataset3=[]
    xscale=1
    
    if column_number==1:
        b='Average Height'
        xscale=1e9
        ylabel='Average Height (nm)'
        units='nm'
    
    elif column_number==2 and method_number==1:
        b='Sa'
        xscale=1e9
        ylabel='Sa (nm)'
        units='nm'   

    elif column_number==2 and method_number==2:
        b='Ra'
        xscale=1e9
        ylabel='Ra (nm)'
        units='nm'
        
    elif column_number==3 and method_number==1:
        b='Sq'
        xscale=1e9
        ylabel='Sq (nm)'
        units='nm' 

    elif column_number==3 and method_number==2:
        b='Rq'
        xscale=1e9
        ylabel='Rq (nm)'
        units='nm'
   
    elif column_number==4:
        b='Skew'
        ylabel='Skew'
        units=''
  
    elif column_number==5:
        b='Kurtosis'
        ylabel='Kurtosis'
        units=''
    
    elif column_number==6 and method_number==1:
        b='Sz'
        xscale=1e9
        ylabel='Sz (nm)'
        units='nm'

    elif column_number==6 and method_number==2:
        b='Rz'
        xscale=1e9
        ylabel='Rz (nm)'
        units='nm'
    
    elif column_number==0 and method_number==0:         #The inputs that choose the label for the bendstrain graph
        ylabel='Bend Strain (%)'
    else:
        pass
    
    for i, samplename in enumerate(sample_names):     
        sampledata_1=roughness_data_1[samplename]
        sampledata_2=roughness_data_2[samplename]
        
        #Correct the scale of the data
        x1= [j * xscale for j in sampledata_1]
        x2= [j * xscale for j in sampledata_2]
        


        mu1, conf_int_95percent1 = bootstrap_confidence_intervals(x1, column_number, numpy.mean, [2.5, 97.5], runs=60000)
        # _, conf_int_99_7percent1 = bootstrap_confidence_intervals(x1, column_number, numpy.mean, [0.15, 99.85])
        roughness_means_1.append(mu1)
        CI_95_dataset1.append([mu1- conf_int_95percent1[0], conf_int_95percent1[1]-mu1])
        # CI_99_7_dataset1.append([mu1- conf_int_99_7percent1[0], conf_int_99_7percent1[1]-mu1])
        
        mu2, conf_int_95percent2 = bootstrap_confidence_intervals(x2, column_number, numpy.mean, [2.5, 97.5], runs=60000)
        # _, conf_int_99_7percent2 = bootstrap_confidence_intervals(x2, column_number, numpy.mean, [0.15, 99.85])
        roughness_means_2.append(mu2)
        CI_95_dataset2.append([mu2- conf_int_95percent2[0], conf_int_95percent2[1]-mu2])
        # CI_99_7_dataset2.append([mu2- conf_int_99_7percent2[0], conf_int_99_7percent2[1]-mu2])
        
        if samplename=='A09':
            sampledata_3=roughness_data_3[samplename]
            x3= [j * xscale for j in sampledata_3]
            mu3, conf_int_95percent3 = bootstrap_confidence_intervals(x3, column_number, numpy.mean, [2.5, 97.5], runs=60000)
            roughness_means_3.append(mu3)
            CI_95_dataset3.append([mu3- conf_int_95percent3[0], conf_int_95percent3[1]-mu3])
            
    # CI_95_array1=numpy.array(CI_95_dataset1)
    # CI_95_array2=numpy.array(CI_95_dataset2)
    
    x_pos = numpy.arange(len(sample_names))
    error1=numpy.transpose(CI_95_dataset1)
    error2=numpy.transpose(CI_95_dataset2)
    error3=numpy.transpose(CI_95_dataset3)
    
    #=============================Plotting the graph ===============================================
    fig, ax = plt.subplots()
    ax.errorbar(x_pos, roughness_means_1, yerr=error1, fmt='x', alpha=0.6, color='blue', ecolor='blue', capsize=10)
    ax.errorbar(x_pos, roughness_means_2, yerr=error2, fmt='x', alpha=0.6, color='red', ecolor='red', capsize=10)
    ax.errorbar(    8, roughness_means_3, yerr=error3, fmt='x', alpha=0.6, color='black', ecolor='black', capsize=10) #Plot values for repeat A09.
    ax.set_ylabel(ylabel, weight='bold')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(sample_names)
    ax.set_xlabel('Samples', weight='bold')
    plt.xlim(-1,len(sample_names))
    if method_number==0:
        ax.set_title(dataset_comparison_title)
    else:
        # ax.set_title('Comparison of Sample Means Between Two Repeat Datasets - '+b)
        pass
    ax.yaxis.grid(True)
    
    # Save the figure and show
    ##Creating the legend
    handles, labels = ax.get_legend_handles_labels()
    # manually define new patches for the different bands 
    u1_band = mpatches.Patch(color='b', alpha=0.6, label=three_dataset_labels[0])
    u2_band = mpatches.Patch(color='r', alpha=0.6, label=three_dataset_labels[1])
    u3_band = mpatches.Patch(color='black', alpha=0.6, label=three_dataset_labels[2])
    handles.append(u1_band) 
    handles.append(u2_band) 
    handles.append(u3_band) 
    plt.legend(handles=handles, prop={"size":8})
    plt.tight_layout()
    plt.show()
    plt.pause(7)
    if savefile==True and method_number==0:
        plt.savefig(dataset_comparison_title +' - '+tagline+'.png')
    elif savefile==True:
        plt.savefig('Comparison of Sample Means - '+ dataset_comparison_title +' - '+b+' ' +tagline+'.png')
    plt.close()

def give_me_a_straight_line(x,y):
    m, c  = numpy.polyfit(x,y,deg=1)
    # line  = m * x + c
    return (m,c)

def comparison_plot(x_data, x_err, y_data, y_err, xlabel, ylabel, figure_title, figure_name, savefile=False):
    fig, ax = plt.subplots()
    ax.errorbar(x_data, y_data, yerr=y_err, xerr=x_err, fmt='x', alpha=1, ecolor='black', capsize=10)
    ax.set_ylabel(ylabel, weight='bold')
    m, c = give_me_a_straight_line(x_data, y_data)
    plt.plot(x_data,m*numpy.asarray(x_data)+c,'r')
    # ax.set_xticks(x_pos)
    # ax.set_xticklabels(sample_names)
    ax.set_xlabel(xlabel, weight='bold')
    # plt.xlim(-1,len(sample_names))
    ax.set_title(figure_title)
    # ax.yaxis.grid(True)
    
    # Save the figure and show
    plt.tight_layout()
    plt.show()
    plt.pause(6)
    if savefile==True:
        plt.savefig(figure_name+'.png')
    plt.close()
    
