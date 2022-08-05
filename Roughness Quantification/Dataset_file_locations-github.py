import os

def dataset_folder_check(path):
    if not os.path.exists(path):
        os.makedirs(path)
        os.makedirs(path+'/Sample Comparison Mean and Uncertainties')
        os.makedirs(path+'/Single Histograms')
        os.makedirs(path+'/Convergence Plots')
        os.makedirs(path+'/Test')
        print('\nCreated folder and subdirectories for: '+path+'\n')
    else:
        # print('\nPath exists\n')
        pass

def get_dataset_info(Method_number=1, dataset_number=3):
    """dataset_number=3 corresponds to the "1st Dataset" 
       dataset_number=4 corresponds to the "2nd Dataset"
       For purposes of naming, the 1st and 2nd datasets are 
       called 3rd and 4th respectively in the code due to obsolete 
       datasets being previously used. 
       """
    
    root_1='F:/3.SiC Blinds/1.ILC/1.Section 1/'
    root_2='F:/3.SiC Blinds/3.ILC Repeat/1.Section 1/'
    graph_root='F:/3.SiC Blinds/1.ILC/Graphs/'
    roughness_par_names=["Distance (um)","Average Height (nm)", "AM, Ra/Sa (nm)",
                       "RMS, Rq/Sq (nm)","Rsk/Ssk","Rku/Sku","Rz/Sz (nm)"]
            
    if dataset_number==3:
        dataset_name='3rd Good Dataset'
        txtfilepathways=['F:/3.SiC Blinds/4.ILC (original, but changing tips)/3.Section 3/A01',\
                'F:/3.SiC Blinds/4.ILC (original, but changing tips)/1.Section 1/A02',\
                'F:/3.SiC Blinds/3.ILC Repeat/3.Section 3/A03',\
                'F:/3.SiC Blinds/4.ILC (original, but changing tips)/1.Section 1/A04',\
                'F:/3.SiC Blinds/4.ILC (original, but changing tips)/1.Section 1/A05',\
                'F:/3.SiC Blinds/4.ILC (original, but changing tips)/1.Section 1/A06',\
                'F:/3.SiC Blinds/4.ILC (original, but changing tips)/1.Section 1/A07',\
                'F:/3.SiC Blinds/4.ILC (original, but changing tips)/4.Section 4/A08',\
                'F:/3.SiC Blinds/3.ILC Repeat/2.Section 2/A09',\
                'F:/3.SiC Blinds/3.ILC Repeat/2.Section 2/A10']
            
    elif dataset_number==4:
        dataset_name='4th Dataset'
        txtfilepathways=['F:/3.SiC Blinds/3.ILC Repeat/3.Section 3/A01',\
                root_2+'A02',\
                'F:/3.SiC Blinds/1.ILC/3.Section 3/A03',\
                root_2+'A04',\
                root_2+'A05',\
                root_2+'A06',\
                root_2+'A07',\
                'F:/3.SiC Blinds/3.ILC Repeat/4.Section 4/A08',\
                'F:/3.SiC Blinds/1.ILC/2.Section 2/A09',\
                'F:/3.SiC Blinds/1.ILC/2.Section 2/A10']
            
    #Corrects the txtfilepaths for the correct method number and defines the correct roughness graph labels    
    if Method_number==1:
        #Add 'SS Flattened' to the end of path list.
        txtfilepathways=[filepath+'/SS Flattened' for filepath in txtfilepathways]
        roughness_par_labels=["Distance (\u03BCm)", "Average Height (nm)", "Sa (nm)", "Sq (nm)", "Ssk", "Sku", "Sz (nm)"]
        roughness_par_titles=["Distance", "Average Height", "Sa", "Sq", "Ssk", "Sku", "Sz"]  
    elif Method_number==2:
        #Add 'Unflattened' to the end of the path list.
        txtfilepathways=[filepath+'/Unflattened' for filepath in txtfilepathways]
        roughness_par_labels=["Distance (\u03BCm)", "Average Height (nm)", "Ra (nm)", "Rq (nm)", "Rsk", "Rku", "Rz (nm)"]
        roughness_par_titles=["Distance", "Average Height", "Ra", "Rq", "Rsk", "Rku", "Rz"]
     
    #Define the file location for the graph saving
    graphsave_root=graph_root+str(dataset_number)+'.'+str(Method_number)+' - '+dataset_name +' - Method '+str(Method_number) + ' Graphs/'
    
    #Checks if the file location exists, if not, creates it. 
    dataset_folder_check(graphsave_root)  

    return (txtfilepathways, graphsave_root, dataset_name, roughness_par_names, roughness_par_labels, roughness_par_titles)
