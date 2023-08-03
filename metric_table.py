
# import required libraries
import os
import pandas as pd
from pathlib import Path
import glob
import numpy as np
import re

# Set the directory path that contains each samples folder, comparison name, and threshold
directory = '/home/mbxha18/allcombinations/New_illumina/exeter/exeter_ill_vs_exeter_ont_NT002/'
comparison = "exeter_ill_vs_exeter_ont_NT002.mincov1.nostrandbiasfilter.mpileup2snp.varscan"
threshold = "optimal"

# Read data from the control file
All_control = pd.read_csv(f'/home/mbxha18/OvsI_prom/ONT_match_prom/All_control_with_indels.tsv',sep="\t")
# Remove indel rows and convert the 'POS' column to int data type
Control_remove_indel_positions = All_control[~All_control['POS'].str.startswith(('del', 'ins'))].astype({'POS': int})
# Define the list of 'misalignment_errors_positions'
misalignment_errors_positions = [21766, 22034, 6515, 21987, 21988, 22196, 28363]
# Filter out rows with positions listed in 'misalignment_errors_positions'
Control = Control_remove_indel_positions[~Control_remove_indel_positions['POS'].isin(misalignment_errors_positions)]

# Read data from the TWIST_synthetic_mixes_sample_sheet that contains the sample information
a_conc = pd.read_csv("/home/mbxha18/allcombinations/controls_old/TWIST_synthetic_mixes_sample_sheet_150222_edited.txt", sep="\t",index_col=0)

# Initialize an empty list to store data for each sample
unfiltered = []

# Iterate through each subdirectory in the 'directory'
all_subdirs = glob.glob(f"{directory}/*/")
for subdir in all_subdirs:
    dir = subdir.split('/')[7]
    print(f"{dir}")
    files = os.listdir(subdir)
   
    # Check if both illumina and ont required varscan files are present in the subdirectory
    ill_file_present = False
    ont_file_present = False
    for file in files:
        if file.endswith(".ill.mpileup2snp.mincov4.minvar0.06.nostrandbiasfilter.varscan.tsv"):
            ill_file_present = True
        elif file.endswith(".ont.mpileup.mincov4.minvar0.06.varscan.tsv"):
            ont_file_present = True
    if ill_file_present and ont_file_present:
        print("Run")

        # Set the working directory to the current subdirectory
        os.chdir(f'{directory}/{dir}')
        # Create sample names for "ill" and "ont"
        sample_ill = f"{dir}_ill"
        sample_ont = f"{dir}_ont"

        # Get a list of needed ont and illumina files and read into a dataframe 
        ont_files = glob.glob(f"{directory}/{dir}/*.ont.mpileup.mincov4.minvar0.06.varscan.tsv") #varscan
        ont = pd.read_csv(ont_files[0], sep="\t")
        ont_pileup = glob.glob(f"{directory}/{dir}/*.ont.pileup") #pileup
        av_cov_ont_files = glob.glob(f"{directory}/{dir}/*.ont.mosdepth.summary.txt") #mosdepth
        av_cov_ont = pd.read_csv(av_cov_ont_files[0], sep="\t")
        ont_depth_files = glob.glob(f"{directory}/{dir}/*.ont.samtools.cov") # samtools
        ont_depth = pd.read_csv(ont_depth_files[0], sep="\t", header = None, names = ["REGION","POS","DEPTH_ont"], usecols = ["POS","DEPTH_ont"])

        ill_files = glob.glob(f"{directory}/{dir}/*.ill.mpileup2snp.mincov4.minvar0.06.nostrandbiasfilter.varscan.tsv") 
        ill = pd.read_csv(ill_files[0], sep="\t")
        ill_pileup = glob.glob(f"{directory}/{dir}/*.ill.pileup")
        av_cov_ill_files = glob.glob(f"{directory}/{dir}/*.ill.mosdepth.summary.txt")
        av_cov_ill = pd.read_csv(av_cov_ill_files[0], sep="\t")
        ill_depth_files = glob.glob(f"{directory}/{dir}/*.ill.samtools.cov")
        ill_depth = pd.read_csv(ill_depth_files[0], sep="\t", header = None, names = ["REGION","POS","DEPTH_ill"], usecols = ["POS","DEPTH_ill"])
        # Check if file is empty
        if ill.empty or ont.empty:
            print("Don't Run")
            continue  # Skip the rest of the code and move to the next iteration if the file is empty
        else:   
            #Convert the varscan output to needed format
            ont[['Cons', 'Cov', 'Reads1', 'Reads2', 'Freq', 'P-value']] = ont['Cons:Cov:Reads1:Reads2:Freq:P-value'].str.split(':', expand=True)
            ont = ont.rename(columns={"Var": "ALT", "Position": "POS", "Ref": "REF", "Reads1": "REF_DP", "Reads2": "ALT_DP", "Freq" : "ALT_F", "Chrom" : "REGION"}) # change column names
            ont['ALT'] = ont['ALT'].map(lambda x: x.lstrip('+-'))  # remove the indels in the varscan file
            ont["ALT_F"] = ont["ALT_F"].str.rstrip('%').astype('float')
            ont = ont.sort_values('ALT_F').drop_duplicates(subset='POS', keep='last', ignore_index= True).sort_values('POS').reset_index(drop=True)
            ont = ont[["POS", "REF", "ALT", "REF_DP", "ALT_DP", "ALT_F"]] # extract needed columns 
            ont = ont.fillna(value={"REF_DP": 0, "ALT_DP": 0, "ALT_F": 0}) # file NAs with 0
    
            # match the sample to its lineage in the TWIST_synthetic_mixes_sample_sheet by matching the directoryname to its respective plate well
            try: 
                for index, row in a_conc.iterrows():
                    if (row["plate_well"]  == dir.split('-')[0]):
                        print(row["lineage"])
                        ont["lineage"] = row["lineage"]
                        sample_lineage = ont.lineage.values.tolist()[0]
            except IndexError:
                pass
            control_for_sample = Control[Control["lineage"].eq(sample_lineage)].drop_duplicates(ignore_index= True) #extract the control for that lineage

            # varscan df modification for the positions with multiple snps [23604 and 28881] using extracted data of the second snp from the pileup file
            pos_values_in_ont = ont['POS'].values  # extract all positions called in the varscan df
            # Define the list of lineages with positions with multiple snps
            msnp_list = ["AlphaDeltaC23", "AlphaDeltaAY.2", "AlphaBetaDeltaC23", "AlphaBetaDeltaC23DeltaAY.2", "DeltaC23Omicron"]
            # check if the sample is one of the msnps lineage and if any of those positions have been called in the varscan df
            if sample_lineage in msnp_list and (23604 in pos_values_in_ont.tolist() or 28881 in pos_values_in_ont.tolist()):
                # Initialize an empty list to store the extracted data of the snps
                data = []
                # Iterate over each pileup file
                for file_name in ont_pileup:
                    # Read the file line by line
                    with open(file_name, 'r') as file:
                        for line in file:
                            # Check if positions contains "23604" or "28881"
                            if re.search(r'23604|28881', line):
                            # Split the line into columns
                                columns = line.split('\t')       
                                # Extract the desired values of snps at "23604" or "28881"
                                pos = pd.to_numeric(columns[1], errors='coerce') 
                                ref = columns[2]
                                depth_ont = pd.to_numeric(columns[3], errors='coerce') 
                                alt_bases = re.findall(r'[A-Za-z]', columns[4].upper()) 
                                alt = set(re.findall(r'[A-Za-z]', columns[4].upper())) # Convert the string to a set of individual characters
                                # Append the extracted values to the data list
                                data.append([pos, ref, depth_ont, alt_bases, alt])
                
                # Create a DataFrame from the extracted data
                snp_ont = pd.DataFrame(data, columns=['POS', 'REF', 'DEPTH_ont','ALT_BASES', 'ALT'])
                # Count occurrences of each character in ALT and store in a new column
                snp_ont['ALT_DP'] = snp_ont['ALT_BASES'].apply(lambda x: {i: x.count(i) for i in set(x)})
                
                #Create individual snp rows and modify to fit existing varscan df
                #Create an empty DataFrame to store the new rows
                new_rows = []
                # Iterate over the rows of msnp_ont DataFrame
                for index, row in snp_ont.iterrows():
                    pos = row["POS"]
                    ref = row["REF"]
                    depth_ont = row ["DEPTH_ont"]
                    alt_options = list(row["ALT"])  # Convert the set to a list to access elements
                    alt_DP = row["ALT_DP"]
                    # Create two new rows with the specified conditions
                    for alt in alt_options:
                        new_row = {
                            "POS": pos,
                            "REF": ref,
                            "ALT": alt,
                            "lineage": sample_lineage,
                            "REF_DP": 0,
                            "ALT_DP": alt_DP[alt],
                            "ALT_F": alt_DP[alt]/ depth_ont * 100.0,
                            "DEPTH_ont": depth_ont
                            }   
                        new_rows.append(new_row) # Append the new rows to the new_rows list
                df_snps = pd.DataFrame(data=new_rows)   # mke it into a DataFrame 
            
                ont = ont[~ont["POS"].isin([23604])] # remove the row with these positions in the varscan df
                ont = ont[~ont["POS"].isin([28881])]
                pos_23604_row = df_snps.loc[df_snps['POS'] == 23604] # Define each rown in the msnp df with the postion number
                pos_28881_row = df_snps.loc[df_snps['POS'] == 28881]
                
                # Check if the POS value in msnp is either 23604 or 28881 and concatenate these positions to varscan df
                if 23604 in pos_values_in_ont.tolist():
                    frames = [ont, pos_23604_row]
                    ont = pd.concat(frames).reset_index(drop=True)
                
                if 28881 in pos_values_in_ont.tolist():
                    frames = [ont, pos_28881_row]
                    ont = pd.concat(frames).reset_index(drop=True)
            else:
                print("No multiple snp")
                pass

            # merge the varscan df with the control for this sample
            ont_true = pd.merge(control_for_sample, ont, how= "left", on=["POS","REF","ALT","lineage"])
            ont_true_depth = pd.merge(ont_true, ont_depth, how= "left", on=["POS"]) # merge the resulting df with the depth file
            if "DEPTH_ont_x" in ont_true_depth.columns: 
                ont_true_depth["DEPTH_ont"] = ont_true_depth["DEPTH_ont_x"].combine_first(ont_true_depth["DEPTH_ont_y"])
                ont_true_depth.drop("DEPTH_ont_x",axis=1, inplace=True)
                ont_true_depth.drop("DEPTH_ont_y",axis=1, inplace=True)
            ont_table = ont_true_depth.fillna(0)
            A = ont_table
            #for optimal threshold, remove snp with cov less than 4 from the missing snp positions added [23604 and 28881], because this will originally not get called in the varscan at that threshold.
            try: 
                if threshold == "optimal":
                    A['ALT_DP'] = pd.to_numeric(A['ALT_DP'], errors='coerce')
                    A = A[(A['POS'] == 23604) | (A['POS'] == 28881) & (A["ALT_DP"] >= 4)]
            except IndexError:
                pass
            print()
         
 


            
            ## same codes for illumina files analysis
            ill[['Cons', 'Cov', 'Reads1', 'Reads2', 'Freq', 'P-value']] = ill['Cons:Cov:Reads1:Reads2:Freq:P-value'].str.split(':', expand=True)
            ill = ill.rename(columns={"Var": "ALT", "Position": "POS", "Ref": "REF", "Reads1": "REF_DP", "Reads2": "ALT_DP", "Freq" : "ALT_F", "Chrom" : "REGION"})
            ill['ALT'] = ill['ALT'].map(lambda x: x.lstrip('+-'))
            ill["ALT_F"] = ill["ALT_F"].str.rstrip('%').astype('float')
            ill = ill.sort_values('ALT_F').drop_duplicates(subset='POS', keep='last', ignore_index= True).sort_values('POS').reset_index(drop=True)
            ill = ill[["POS", "REF", "ALT", "REF_DP", "ALT_DP", "ALT_F"]]
            ill = ill.fillna(value={"REF_DP": 0, "ALT_DP": 0, "ALT_F": 0})
    
            try: 
                for index, row in a_conc.iterrows():
                    if (row["plate_well"]  == dir.split('-')[0]):
                        ill["lineage"] = row["lineage"]
                        sample_lineage = ill.lineage.values.tolist()[0]
            except IndexError:
                pass
            control_for_sample = Control[Control["lineage"].eq(sample_lineage)].drop_duplicates(ignore_index= True)
            
            pos_values_in_ill = ill['POS'].values
            msnp_list = ["AlphaDeltaC23", "AlphaDeltaAY.2", "AlphaBetaDeltaC23", "AlphaBetaDeltaC23DeltaAY.2", "DeltaC23Omicron"]
            
            if sample_lineage in msnp_list and (23604 in pos_values_in_ill.tolist() or 28881 in pos_values_in_ill.tolist()):
                data = []
                for file_name in ill_pileup:
                    with open(file_name, 'r') as file:
                        for line in file:
                            if re.search(r'23604|28881', line):
                                columns = line.split('\t')
                                pos = pd.to_numeric(columns[1], errors='coerce') 
                                ref = columns[2]
                                depth_ill = pd.to_numeric(columns[3], errors='coerce') 
                                alt_bases = re.findall(r'[A-Za-z]', columns[4].upper()) 
                                alt = set(re.findall(r'[A-Za-z]', columns[4].upper())) #
                                
                                data.append([pos, ref, depth_ill, alt_bases, alt])
                
                snp_ill = pd.DataFrame(data, columns=['POS', 'REF', 'DEPTH_ill','ALT_BASES', 'ALT'])
                snp_ill['ALT_DP'] = snp_ill['ALT_BASES'].apply(lambda x: {i: x.count(i) for i in set(x)})
          
                for index, row in snp_ill.iterrows():
                    pos = row["POS"]
                    ref = row["REF"]
                    depth_ill = row ["DEPTH_ill"]
                    alt_options = list(row["ALT"]) 
                    alt_DP = row["ALT_DP"]
                    for alt in alt_options:
                        new_row = {
                            "POS": pos,
                            "REF": ref,
                            "ALT": alt,
                            "lineage": sample_lineage,
                            "REF_DP": 0,
                            "ALT_DP": alt_DP[alt],
                            "ALT_F": alt_DP[alt]/ depth_ill * 100.0,
                            "DEPTH_ill": depth_ill
                            }   
                        new_rows.append(new_row)
                df_snps = pd.DataFrame(data=new_rows)
            
                ill = ill[~ill["POS"].isin([23604])]
                ill = ill[~ill["POS"].isin([28881])]
                pos_23604_row = df_snps.loc[df_snps['POS'] == 23604]
                pos_28881_row = df_snps.loc[df_snps['POS'] == 28881]
                
                if 23604 in pos_values_in_ill.tolist():
                    frames = [ill, pos_23604_row]
                    ill = pd.concat(frames).reset_index(drop=True)
                
                if 28881 in pos_values_in_ill.tolist():
                    frames = [ill, pos_28881_row]
                    ill = pd.concat(frames).reset_index(drop=True)
            else:
                pass
                
            ill_true = pd.merge(control_for_sample, ill, how= "left", on=["POS","REF","ALT","lineage"])
            ill_true_depth = pd.merge(ill_true, ill_depth, how= "left", on=["POS"])
            if "DEPTH_ill_x" in ill_true_depth.columns: 
                ill_true_depth["DEPTH_ill"] = ill_true_depth["DEPTH_ill_x"].combine_first(ill_true_depth["DEPTH_ill_y"])
                ill_true_depth.drop("DEPTH_ill_x",axis=1, inplace=True)
                ill_true_depth.drop("DEPTH_ill_y",axis=1, inplace=True)
            ill_table = ill_true_depth.fillna(0)
            B = ill_table
            try: 
                if threshold == "optimal":
                    B['ALT_DP'] = pd.to_numeric(B['ALT_DP'], errors='coerce')
                    B = B[(B['POS'] == 23604) | (B['POS'] == 28881) & (B["ALT_DP"] >= 4)]
            except IndexError:
                pass
     

            
            # merge both df A(ont) and B(ill)
            All = pd.merge(A, B, how="outer", on=["POS","REF","ALT","lineage"], suffixes=("_ont", "_ill"))
            # Check if the DataFrame "All" is not empty
            if not All.empty:
                # Define a function to categorize rows in the DataFrame based on coverage  
                def categorise_ill(row): 
                    if row['lineage'] == "neg":
                        return 'neg'
                    if row['ALT_DP_ill'] == 0 and row["DEPTH_ill"] < 1: 
                        return 'no_coverage'
                    elif row['ALT_DP_ill']== 0 and row["DEPTH_ill"] >= 1 :
                        return 'not_called'
                    else: 
                        return 'called'
                # Apply the function to the "All" DataFrame, creating a new column "REASON"
                All["REASON_ill"] = All.apply(lambda row: categorise_ill(row), axis=1)
                        
                def categorise_ont(row): 
                    if row['lineage'] == "neg":
                        return 'neg'
                    if row['ALT_DP_ont'] == 0 and row["DEPTH_ont"] < 1: 
                        return 'no_coverage'
                    elif row['ALT_DP_ont']== 0 and row["DEPTH_ont"] >= 1 :
                        return 'not_called'
                    else: 
                        return 'called'
                All["REASON_ont"] = All.apply(lambda row: categorise_ont(row), axis=1)
                # Save the modified "All" DataFrame to a tsv file for each sample
                All.to_csv(f'{directory}/{dir}/{dir}.variant.{comparison}.varscan.tsv',sep="\t",index = False)
    
                # Calculate the proportion of the genome with minimum coverage 
                count = (ill_depth['DEPTH_ill'] > 3).sum()
                proportion_ill = (count/29903) * 100
                    
                count = (ont_depth['DEPTH_ont'] > 3).sum()
                proportion_ont = (count/29903) * 100
    
                # Define TP, FP, FN, TN, FNwc as metrics for classification performance
                # TP = candidates identified by truth and ONT as true
                # FP = candidates identified by ONT as true but not true
                # FN = candidates not identified by ONT as true
                # TN = candidates not identified by both truth and ONT
                # FNwc = candidates not identified by ONT and have coverage
        
                # Filter the "All" DataFrame to include only positions with minimum coverage for both ONT and ILL datasets
                All_Cov = All[~All['REASON_ont'].isin(["no_coverage"]) ]
                All_Cov2 = All_Cov[~All_Cov['REASON_ill'].isin(["no_coverage"]) ]
               
                # Loop through rows in DataFrame to calculate various metrics 
                for index, row in ont_table.iterrows():
                    if (row["lineage"]  == "neg"):
                        True_variant = 0
                        TP_ont = 0
                        FP_ont = 0
                        FN_ont = 0
                        TN_ont = 0
                        FNwc_ont = 0
                    else:
                        True_variant = len(control_for_sample["POS"]) 
                        TP_ont = len(control_for_sample["POS"]) - len(A[A["ALT_DP"].eq(0)])
                        FP_ont = len(ont["POS"].unique()) - len(A[A["ALT_DP"].ne(0)])
                        FN_ont = True_variant - len(A[A["ALT_DP"].ne(0)])
                        TN_ont = 29903 - len(ont["POS"].unique()) - FN_ont
                        FNwc_ont = len(All_Cov2[All_Cov2["REASON_ont"].eq("not_called")]) 
                # Calculate Sensitivity, Precision, and Jaccard similarity 
                try: 
                    Sensitivity_ONT = TP_ont/(TP_ont + FN_ont)
                    Precision_ONT = TP_ont/(TP_ont + FP_ont)
                    Jaccard_similarity_ONT = TP_ont/(TP_ont + FP_ont + FN_ont)
    
                    Sensitivity_ONT_WC = len(All_Cov2[All_Cov2["REASON_ont"].eq("called")])/len(All_Cov2)
                    Jaccard_similarity_ONT_WC = len(All_Cov2[All_Cov2["REASON_ont"].eq("called")])/(len(All_Cov2) + FNwc_ont)
    
                except ZeroDivisionError:
                    Precision_ONT = 0
                    Sensitivity_ONT = 0
                    Jaccard_similarity_ONT = 0
                    Sensitivity_ONT_WC = 0
                    Jaccard_similarity_ONT_WC = 0
    
                for index, row in ill_table.iterrows():
                    if (row["lineage"]  == "neg"):
                        True_variant = 0
                        TP_ill = 0
                        FP_ill = 0
                        FN_ill = 0
                        TN_ill = 0
                        FNwc_ill = 0
                    else:
                        True_variant = len(control_for_sample["POS"])
                        TP_ill = len(control_for_sample["POS"]) - len(B[B["ALT_DP"].eq(0)])
                        FP_ill = len(ill["POS"].unique()) - len(B[B["ALT_DP"].ne(0)])
                        FN_ill = True_variant - len(B[B["ALT_DP"].ne(0)])
                        TN_ill = 29903 - len(ill["POS"].unique()) - FN_ill
                        FNwc_ill = len(All_Cov2[All_Cov2["REASON_ill"].eq("not_called")]) 
                
    
                try:
                    Sensitivity_ILL = TP_ill/(TP_ill + FN_ill)
                    Precision_ILL = TP_ill/(TP_ill + FP_ill)
                    Jaccard_similarity_ILL = TP_ill/(TP_ill + FP_ill + FN_ill)
    
                    Sensitivity_ILL_WC = len(All_Cov2[All_Cov2["REASON_ill"].eq("called")])/len(All_Cov2)
                    Jaccard_similarity_ILL_WC = len(All_Cov2[All_Cov2["REASON_ill"].eq("called")])/(len(All_Cov2) + FNwc_ill)
    
                except ZeroDivisionError:
                    Precision_ILL = 0
                    Sensitivity_ILL = 0
                    Jaccard_similarity_ILL = 0
                    Sensitivity_ILL_WC = 0
                    Jaccard_similarity_ILL_WC = 0

            # Loop through rows in "a_conc" DataFrame to extract the lineage concentration in the sample
            for index, row in a_conc.iterrows():
                if (row["plate_well"]  == dir.split('-')[0]):
                    alpha_conc = row["Alpha-C15"]
                    alpha_exp_frequency = row ["alpha_exp_frequency"]
                    beta_conc = row["Beta-C16"]
                    beta_exp_frequency = row["beta_exp_frequency"]
                    delta_conc = row["Delta-C23"]
                    delta_exp_frequency = row["deltaC23_exp_frequency"]
                    deltaAY2_conc = row["DeltaAY.2-C29"]
                    deltaAY2_exp_frequency = row["deltaAY2_exp_frequency"]
                    omicron_conc = row["Omicron-C48"]
                    omicron_frequency = row["omicron_exp_frequency"]
                    lineage = row["lineage"]
            # Create a dictionary "d" to store all the calculated metrics and extracted values for each dataset and lineage
            d = {'sample': [sample_ont, sample_ill],'True_variant': [True_variant, True_variant], 'TP':[TP_ont,TP_ill], 'FP':[FP_ont, FP_ill], 'FN':[FN_ont,FN_ill], 'TN':[TN_ont,TN_ill], 'FN_wc':[FNwc_ont,FNwc_ill], 'sensitivity': [Sensitivity_ONT, Sensitivity_ILL], 'precision': [Precision_ONT, Precision_ILL], 'Jaccard_similarity': [Jaccard_similarity_ONT, Jaccard_similarity_ILL],'sensitivity_wc': [Sensitivity_ONT_WC, Sensitivity_ILL_WC], 'Jaccard_similarity_wc': [Jaccard_similarity_ONT_WC, Jaccard_similarity_ILL_WC],'lineage' : [lineage, lineage], 'alpha_conc': [alpha_conc, alpha_conc], 'beta_conc': [beta_conc,beta_conc], 'delta_conc': [delta_conc,delta_conc],'deltaAY.2_conc':[deltaAY2_conc,deltaAY2_conc], 'omicron_conc': [omicron_conc,omicron_conc], 'alpha_exp_frequency': [alpha_exp_frequency,alpha_exp_frequency],'beta_exp_frequency': [beta_exp_frequency,beta_exp_frequency], 'delta_exp_frequency':[delta_exp_frequency,delta_exp_frequency], 'deltaAY2_exp_frequency':[deltaAY2_exp_frequency,deltaAY2_exp_frequency], 'omicron_frequency' : [omicron_frequency,omicron_frequency], 'mean_cov': [av_cov_ont.loc[0]['mean'],av_cov_ill.loc[0]['mean']],'proportion': [proportion_ont,proportion_ill]}
            
            df = pd.DataFrame(data=d) # Create a DataFrame "df" using the dictionary "d"
            unfiltered.append(df)  # Append "df" to the "unfiltered" list
    # If the "All" DataFrame is empty, skip this block and do nothing
    else:
        pass    
# Concatenate all DataFrames in the "unfiltered" list into a single DataFrame 
unfiltered_table = pd.concat(unfiltered).reset_index(drop=True)
# Remove rows with 'lineage' value of 'neg' that is negative samples
unfiltered_table_no_neg = unfiltered_table[unfiltered_table['lineage'] != 'neg']
print (unfiltered_table_no_neg)
# Save the final DataFrame to a TSV file
unfiltered_table_no_neg.to_csv(f'/home/mbxha18/ww_paper/Baseline/bigtable.{comparison}.tsv',sep="\t")



