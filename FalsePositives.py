# import required libraries
import os
import pandas as pd
from pathlib import Path
import glob
import numpy as np
import re

# Set the directory path that contains each samples folder, comparison name, and threshold
directory = '/home/mbxha18/allcombinations/New_illumina/exeter/exeter_ill_vs_exeter_ont_NT002'
table = "FP_table_NT002"
threshold = "baseline"
# Read data from the control file
All_control = pd.read_csv(f'/home/mbxha18/OvsI_prom/ONT_match_prom/All_control_with_indels.tsv',sep="\t")
# remove indel positions
Control_remove_indel_positions = All_control[~All_control['POS'].str.startswith(('del', 'ins'))].astype({'POS': int})
# remove misalignment_errors_positions 
misalignment_errors_positions = [21766, 22034, 6515, 21987, 21988, 22196, 28363]
Control = Control_remove_indel_positions[~Control_remove_indel_positions['POS'].isin(misalignment_errors_positions)]
# Read data from the TWIST_synthetic_mixes_sample_sheet that contains the sample information
a_conc = pd.read_csv("/home/mbxha18/allcombinations/controls_old/TWIST_synthetic_mixes_sample_sheet_150222_edited.txt", sep="\t",index_col=0)

# Initialize an empty list to store data for each sample
FP_ont_table=[]
FP_ill_table=[]

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
        if file.endswith(".ill.mincov1.mpileup2snp.nostrandbiasfilter.varscan.tsv"):
            ill_file_present = True
        elif file.endswith(".ill.mincov1.mpileup2snp.varscan.tsv"):
            ont_file_present = True
    
    if ill_file_present and ont_file_present:
        print("Run")
        # Set the working directory to the current subdirectory
        os.chdir(f'{directory}/{dir}')
        sample_ill = f"{dir}_ill"
        sample_ont = f"{dir}_ont"

        # Get a list of needed files and read them in dfs
        ont_files = glob.glob(f"{directory}/{dir}/*.ont.mincov1.mpileup2snp.varscan.tsv") 
        print(ont_files)
        ont = pd.read_csv(ont_files[0], sep="\t")
        ont_pileup = glob.glob(f"{directory}/{dir}/*.ont.pileup")
        av_cov_ont_files = glob.glob(f"{directory}/{dir}/*.ont.mosdepth.summary.txt")
        av_cov_ont = pd.read_csv(av_cov_ont_files[0], sep="\t")
        ont_depth_files = glob.glob(f"{directory}/{dir}/*.ont.samtools.cov")
        ont_depth = pd.read_csv(ont_depth_files[0], sep="\t", header = None, names = ["REGION","POS","DEPTH_ont"], usecols = ["POS","DEPTH_ont"])

        ill_files = glob.glob(f"{directory}/{dir}/*.ill.mincov1.mpileup2snp.nostrandbiasfilter.varscan.tsv") 
        print(ont_files)
        ill = pd.read_csv(ill_files[0], sep="\t")
        ill_pileup = glob.glob(f"{directory}/{dir}/*.ill.pileup")
        av_cov_ill_files = glob.glob(f"{directory}/{dir}/*.ill.mosdepth.summary.txt")
        av_cov_ill = pd.read_csv(av_cov_ill_files[0], sep="\t")
        ill_depth_files = glob.glob(f"{directory}/{dir}/*.ill.samtools.cov")
        ill_depth = pd.read_csv(ill_depth_files[0], sep="\t", header = None, names = ["REGION","POS","DEPTH_ill"], usecols = ["POS","DEPTH_ill"])
        # Check if file is empty
        if ill.empty or ont.empty:
            print("Don't Run")
            continue  # Skip the rest of the code and move to the next iteration
        else:   
            #Convert te varscan output to needed format
            ont[['Cons', 'Cov', 'Reads1', 'Reads2', 'Freq', 'P-value']] = ont['Cons:Cov:Reads1:Reads2:Freq:P-value'].str.split(':', expand=True)
            ont = ont.rename(columns={"Var": "ALT", "Position": "POS", "Ref": "REF", "Reads1": "REF_DP", "Reads2": "ALT_DP", "Freq" : "ALT_F", "Chrom" : "REGION"})
            ont['ALT'] = ont['ALT'].map(lambda x: x.lstrip('+-'))
            ont["ALT_F"] = ont["ALT_F"].str.rstrip('%').astype('float')
            ont = ont.sort_values('ALT_F').drop_duplicates(subset='POS', keep='last', ignore_index= True).sort_values('POS').reset_index(drop=True)
            ont = ont[["POS", "REF", "ALT", "REF_DP", "ALT_DP", "ALT_F"]]
            ont = ont.fillna(value={"REF_DP": 0, "ALT_DP": 0, "ALT_F": 0})
    
            #match the sample to its lineage in the TWIST_synthetic_mixes_sample_sheet
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
                data = [] # Initialize an empty list to store the extracted data of the snps
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
                        new_rows.append(new_row)  # Append the new rows to the new_rows lis
                df_snps = pd.DataFrame(data=new_rows) # make it into a DataFrame 
            
                ont = ont[~ont["POS"].isin([23604])] # remove the rows with these positions in the varscan df
                ont = ont[~ont["POS"].isin([28881])]
                pos_23604_row = df_snps.loc[df_snps['POS'] == 23604]  # Define each rows in the msnp df with the postion number
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
            # merge the varscan df with depth file and the control for this sample
            ont_depth = pd.merge(ont, ont_depth, how= "left", on=["POS"])
            if "DEPTH_ont_x" in ont_depth.columns: 
                ont_depth["DEPTH_ont"] = ont_depth["DEPTH_ont_x"].combine_first(ont_depth["DEPTH_ont_y"])
                ont_depth.drop("DEPTH_ont_x",axis=1, inplace=True)
                ont_depth.drop("DEPTH_ont_y",axis=1, inplace=True)
            ont_depth_true = pd.merge(control_for_sample, ont_depth, on=["POS","REF","ALT","lineage"], how="outer", indicator=True)
            FP_O = ont_depth_true[(ont_depth_true["_merge"] != "both") & (ont_depth_true["_merge"] != "left_only")]
            FP_ont = FP_O.drop("_merge", axis=1)
            FP_ont["ALT_DP"] = pd.to_numeric(FP_ont["ALT_DP"], errors="coerce")
            #for optimal threshold, remove snp with cov less than 4 from the missing snp positions added [23604 and 28881], because this will originally not get called in the varscan at that threshold.
            try: 
                if threshold == "optimal":
                    FP_ont = FP_ont[(FP_ont['POS'] == 23604) | (FP_ont['POS'] == 28881) & (FP_ont["ALT_DP"] >= 4)]
            except IndexError:
                pass
     
 


            
            ## same codes for illumina files analysis
            ill[['Cons', 'Cov', 'Reads1', 'Reads2', 'Freq', 'P-value']] = ill['Cons:Cov:Reads1:Reads2:Freq:P-value'].str.split(':', expand=True)
            ill = ill.rename(columns={"Var": "ALT", "Position": "POS", "Ref": "REF", "Reads1": "REF_DP", "Reads2": "ALT_DP", "Freq" : "ALT_F", "Chrom" : "REGION"})
            ill['ALT'] = ill['ALT'].map(lambda x: x.lstrip('+-'))
            ill["ALT_F"] = ill["ALT_F"].str.rstrip('%').astype('float')
            ill = ill.sort_values('ALT_F').drop_duplicates(subset='POS', keep='last', ignore_index= True).sort_values('POS').reset_index(drop=True)
            ill = ill[["POS", "REF", "ALT", "REF_DP", "ALT_DP", "ALT_F"]]
            ill = ill.fillna(value={"REF_DP": 0, "ALT_DP": 0, "ALT_F": 0})
    
            #match the sample to its lineage in the TWIST_synthetic_mixes_sample_sheet
            try: 
                for index, row in a_conc.iterrows():
                    if (row["plate_well"]  == dir.split('-')[0]):
                        ill["lineage"] = row["lineage"]
                        sample_lineage = ill.lineage.values.tolist()[0]
            except IndexError:
                pass
            control_for_sample = Control[Control["lineage"].eq(sample_lineage)].drop_duplicates(ignore_index= True)
            
            pos_values_in_ill = ill['POS'].values
            # pileup modification
            # extract needed rows from pileup file of targed lineages
            # Initialize an empty list to store the extracted data
            msnp_list = ["AlphaDeltaC23", "AlphaDeltaAY.2", "AlphaBetaDeltaC23", "AlphaBetaDeltaC23DeltaAY.2", "DeltaC23Omicron"]
            
            if sample_lineage in msnp_list and (23604 in pos_values_in_ill.tolist() or 28881 in pos_values_in_ill.tolist()):
                data = []
                # Iterate over each file
                for file_name in ill_pileup:
                    # Read the file line by line
                    with open(file_name, 'r') as file:
                        for line in file:
                            # Check if the value in column 1 contains "23604" or "28881"
                            if re.search(r'23604|28881', line):
                            # Split the line into columns
                                columns = line.split('\t')
                                # Extract the desired values
                                pos = pd.to_numeric(columns[1], errors='coerce') 
                                ref = columns[2]
                                depth_ill = pd.to_numeric(columns[3], errors='coerce') 
                                alt_bases = re.findall(r'[A-Za-z]', columns[4].upper()) 
                                alt = set(re.findall(r'[A-Za-z]', columns[4].upper())) # Convert the string to a set of individual characters
                                # Append the extracted values to the data list
                                data.append([pos, ref, depth_ill, alt_bases, alt])
                
                # Create a DataFrame from the extracted data
                snp_ill = pd.DataFrame(data, columns=['POS', 'REF', 'DEPTH_ill','ALT_BASES', 'ALT'])
                # Count occurrences of each character in ALT and store in a new column
                snp_ill['ALT_DP'] = snp_ill['ALT_BASES'].apply(lambda x: {i: x.count(i) for i in set(x)})
                
                #Create individual snp rows and modify to fit existing varscan table
                #Create an empty DataFrame to store the new rows
                new_rows = []
                # Iterate over the rows of msnp_ill DataFrame
                for index, row in snp_ill.iterrows():
                    pos = row["POS"]
                    ref = row["REF"]
                    depth_ill = row ["DEPTH_ill"]
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
                            "ALT_F": alt_DP[alt]/ depth_ill * 100.0,
                            "DEPTH_ill": depth_ill
                            }   
                        new_rows.append(new_row)
                df_snps = pd.DataFrame(data=new_rows)
            
                ill = ill[~ill["POS"].isin([23604])]
                ill = ill[~ill["POS"].isin([28881])]
                pos_23604_row = df_snps.loc[df_snps['POS'] == 23604]
                pos_28881_row = df_snps.loc[df_snps['POS'] == 28881]
                
                        # Check if the POS value in msnp_ill is either 23604 or 28881
                if 23604 in pos_values_in_ill.tolist():
                    frames = [ill, pos_23604_row]
                    ill = pd.concat(frames).reset_index(drop=True)
                
                if 28881 in pos_values_in_ill.tolist():
                    frames = [ill, pos_28881_row]
                    ill = pd.concat(frames).reset_index(drop=True)
            else:
                pass

                
            ill_depth = pd.merge(ill, ill_depth, how= "left", on=["POS"])
            if "DEPTH_ill_x" in ill_depth.columns: 
                ill_depth["DEPTH_ill"] = ill_depth["DEPTH_ill_x"].combine_first(ill_depth["DEPTH_ill_y"])
                ill_depth.drop("DEPTH_ill_x",axis=1, inplace=True)
                ill_depth.drop("DEPTH_ill_y",axis=1, inplace=True)
            ill_depth_true = pd.merge(control_for_sample, ill_depth, on=["POS","REF","ALT","lineage"], how="outer", indicator=True)
            FP_O = ill_depth_true[(ill_depth_true["_merge"] != "both") & (ill_depth_true["_merge"] != "left_only")]
            FP_ill = FP_O.drop("_merge", axis=1)
            FP_ill["ALT_DP"] = pd.to_numeric(FP_ill["ALT_DP"], errors="coerce")
            try: 
                if threshold == "optimal":
                    #FP_ill = FP_ill[FP_ill["ALT_DP"] >= 4]
                    FP_ill = FP_ill[(FP_ill['POS'] == 23604) | (FP_ill['POS'] == 28881) & (FP_ill["ALT_DP"] >= 4)]
            except IndexError:
                pass

            # create new columnin the FN table to contain the directory name(sample name)
            FP_ont["source"] = dir
            FP_ill["source"] = dir
            # Loop through rows in "a_conc" DataFrame to extract the sample concentration 
            for index, row in a_conc.iterrows():
                if row["plate_well"] == dir.split('-')[0]:
                    alpha_conc = row["Alpha-C15"]
                    alpha_exp_frequency = row["alpha_exp_frequency"]
                    beta_conc = row["Beta-C16"]
                    beta_exp_frequency = row["beta_exp_frequency"]
                    delta_conc = row["Delta-C23"]
                    delta_exp_frequency = row["deltaC23_exp_frequency"]
                    deltaAY2_conc = row["DeltaAY.2-C29"]
                    deltaAY2_exp_frequency = row["deltaAY2_exp_frequency"]
                    omicron_conc = row["Omicron-C48"]
                    omicron_frequency = row["omicron_exp_frequency"]
            
                    # Create new columns in FP_ont table and assign values
                    FP_ont["alpha_conc"] = alpha_conc
                    FP_ont["beta_conc"] = beta_conc
                    FP_ont["delta_conc"] = delta_conc
                    FP_ont["deltaAY2_conc"] = deltaAY2_conc
                    FP_ont["omicron_conc"] = omicron_conc
                    FP_ont["alpha_exp_frequency"] = alpha_exp_frequency
                    FP_ont["beta_exp_frequency"] = beta_exp_frequency
                    FP_ont["delta_exp_frequency"] = delta_exp_frequency
                    FP_ont["deltaAY2_exp_frequency"] = deltaAY2_exp_frequency
                    FP_ont["omicron_frequency"] = omicron_frequency
            
                    # Fill the entire columns with the assigned values
                    FP_ont["alpha_conc"].fillna(alpha_conc, inplace=True)
                    FP_ont["beta_conc"].fillna(beta_conc, inplace=True)
                    FP_ont["delta_conc"].fillna(delta_conc, inplace=True)
                    FP_ont["deltaAY2_conc"].fillna(deltaAY2_conc, inplace=True)
                    FP_ont["omicron_conc"].fillna(omicron_conc, inplace=True)
                    FP_ont["alpha_exp_frequency"].fillna(alpha_exp_frequency, inplace=True)
                    FP_ont["beta_exp_frequency"].fillna(beta_exp_frequency, inplace=True)
                    FP_ont["delta_exp_frequency"].fillna(delta_exp_frequency, inplace=True)
                    FP_ont["deltaAY2_exp_frequency"].fillna(deltaAY2_exp_frequency, inplace=True)
                    FP_ont["omicron_frequency"].fillna(omicron_frequency, inplace=True)
        
                    # Create new columns in FP_ill table and assign values
                    FP_ill["alpha_conc"] = alpha_conc
                    FP_ill["beta_conc"] = beta_conc
                    FP_ill["delta_conc"] = delta_conc
                    FP_ill["deltaAY2_conc"] = deltaAY2_conc
                    FP_ill["omicron_conc"] = omicron_conc
                    FP_ill["alpha_exp_frequency"] = alpha_exp_frequency
                    FP_ill["beta_exp_frequency"] = beta_exp_frequency
                    FP_ill["delta_exp_frequency"] = delta_exp_frequency
                    FP_ill["deltaAY2_exp_frequency"] = deltaAY2_exp_frequency
                    FP_ill["omicron_frequency"] = omicron_frequency
            
                    # Fill the entire columns with the assigned values
                    FP_ill["alpha_conc"].fillna(alpha_conc, inplace=True)
                    FP_ill["beta_conc"].fillna(beta_conc, inplace=True)
                    FP_ill["delta_conc"].fillna(delta_conc, inplace=True)
                    FP_ill["deltaAY2_conc"].fillna(deltaAY2_conc, inplace=True)
                    FP_ill["omicron_conc"].fillna(omicron_conc, inplace=True)
                    FP_ill["alpha_exp_frequency"].fillna(alpha_exp_frequency, inplace=True)
                    FP_ill["beta_exp_frequency"].fillna(beta_exp_frequency, inplace=True)
                    FP_ill["delta_exp_frequency"].fillna(delta_exp_frequency, inplace=True)
                    FP_ill["deltaAY2_exp_frequency"].fillna(deltaAY2_exp_frequency, inplace=True)
                    FP_ill["omicron_frequency"].fillna(omicron_frequency, inplace=True)
            
            FP_ont_table.append(FP_ont) # Append the df to the FP_ont_table list
            FP_ill_table.append(FP_ill) # Append the df to the FP_ill_table list
    else:
        pass 
        
FP_ONT_table = pd.concat(FP_ont_table).reset_index(drop=True) # Concatenate all DataFrames in the FP_ont_table list into a single DataFrame 
FP_ILL_table = pd.concat(FP_ill_table).reset_index(drop=True) # Concatenate all DataFrames in the FP_ill_table list into a single DataFrame 
# merge both tables       
FP_table = pd.merge(FP_ONT_table, FP_ILL_table, on=["POS","REF","ALT","lineage","source","alpha_conc","beta_conc","delta_conc","deltaAY2_conc","omicron_conc","alpha_exp_frequency","beta_exp_frequency","delta_exp_frequency","deltaAY2_exp_frequency","omicron_frequency"], how="outer", suffixes=("_ont", "_ill"))
# Save the final DataFrame to a TSV file       
FP_table.to_csv(f'{directory}/{table}.txt',sep="\t",index = False)
        
        
      
        
        
        








































# import os
# import pandas as pd
# from pathlib import Path
# import glob

# ## nostrandbiasfilter
# directory = '/home/mbxha18/allcombinations/New_illumina/exeter/exeter_ill_vs_exeter_ont_NT001/'

# All_control = pd.read_csv(f'/home/mbxha18/OvsI_prom/ONT_match_prom/All_control.tsv',sep="\t")
# lineage_control = pd.read_csv(f'/home/mbxha18/OvsI_prom/ONT_match_prom/lineage_control.txt',sep="\t")
# a_conc = pd.read_csv("/home/mbxha18/allcombinations/controls_old/TWIST_synthetic_mixes_sample_sheet_150222_edited.txt", sep="\t",index_col=0)

# FP_ont_table=[]
# FP_ill_table=[]
# all_subdirs = glob.glob(f"{directory}/*/")

# for subdir in all_subdirs:
#     dir = subdir.split('/')[7]
#     print(f"{dir}")
    
#     os.chdir(f'{directory}/{dir}')
#     sample_ont = f"{dir}_ont"
#     sample_ill = f"{dir}_ill"


#     ont_files = glob.glob(f"{directory}/{dir}/*.ont.mincov1.mpileup2snp.varscan.tsv") 
#     ont = pd.read_csv(ont_files[0], sep="\t")
        
#     ont_depth_files = glob.glob(f"{directory}/{dir}/*.ont.samtools.cov")
#     ont_depth = pd.read_csv(ont_depth_files[0], sep="\t", header = None, names = ["REGION","POS","DEPTH"], usecols = ["POS","DEPTH"])
    
#     if ont.empty == False:
#         ont[['Cons', 'Cov', 'Reads1', 'Reads2', 'Freq', 'P-value']] = ont['Cons:Cov:Reads1:Reads2:Freq:P-value'].str.split(':', expand=True)
#     else:
#         ont[['Reads1','Reads2','Freq']] = " "

#     ont = ont.rename(columns={"Var": "ALT", "Position": "POS", "Ref": "REF", "Reads1": "REF_DP", "Reads2": "ALT_DP", "Freq" : "ALT_F", "Chrom" : "REGION"})
#     ont['ALT'] = ont['ALT'].map(lambda x: x.lstrip('+-'))
#     ont["ALT_F"] = ont["ALT_F"].str.rstrip('%').astype('float')
#     ont = ont.sort_values('ALT_F').drop_duplicates(subset='POS', keep='last', ignore_index= True).sort_values('POS').reset_index(drop=True)
#     ont = ont[["POS","REF","ALT","REF_DP","ALT_DP", "ALT_F"]]
    
#     if ont.empty == True:
#         new_row = pd.DataFrame(data= {'POS':0, 'REF':0, 'ALT':0, 'REF_DP':0, "ALT_DP":0, "ALT_F":0 }, index=[0])
#         ont = pd.concat([ont, new_row]).reset_index(drop=True)

#     try: 
#         for index, row in a_conc.iterrows():
#             if (row["plate_well"]  == dir.split('-')[0]):
#                 print(row["lineage"])
#                 ont["lineage"] = row["lineage"]
#                 sample_lineage = ont.lineage.values.tolist()[0]
#     except IndexError:
#         pass

#     control_for_sample = All_control[All_control["lineage"].eq(sample_lineage)].drop_duplicates(ignore_index= True)

#     ont_depth = pd.merge(ont, ont_depth, how= "left", on=["POS"])

#     ont_depth_true = pd.merge(control_for_sample, ont_depth, on=["POS","REF","ALT","lineage"], how="outer", indicator=True)
#     FP_O = ont_depth_true[(ont_depth_true["_merge"] != "both") & (ont_depth_true["_merge"] != "left_only")]
#     FP_ont = FP_O.drop("_merge", axis=1)



#     ## nostrandbiasfilter
#     ill_files = glob.glob(f"{directory}/{dir}/*.ill.mincov1.mpileup2snp.nostrandbiasfilter.varscan.tsv")
#     ill = pd.read_csv(ill_files[0], sep="\t")
        
#     ill_depth_files = glob.glob(f"{directory}/{dir}/*.ill.samtools.cov")
#     ill_depth = pd.read_csv(ill_depth_files[0], sep="\t", header = None, names = ["REGION","POS","DEPTH"], usecols = ["POS","DEPTH"])
    
#     if ill.empty == False:
#         ill[['Cons', 'Cov', 'Reads1', 'Reads2', 'Freq', 'P-value']] = ill['Cons:Cov:Reads1:Reads2:Freq:P-value'].str.split(':', expand=True)
#     else:
#         ill[['Reads1','Reads2','Freq']] = " "
    
#     ill = ill.rename(columns={"Var": "ALT", "Position": "POS", "Ref": "REF", "Reads1": "REF_DP", "Reads2": "ALT_DP", "Freq" : "ALT_F"})
#     ill['ALT'] = ill['ALT'].map(lambda x: x.lstrip('+-'))
#     ill["ALT_F"] = ill["ALT_F"].str.rstrip('%').astype('float')
#     ill = ill.sort_values('ALT_F').drop_duplicates(subset='POS', keep='last', ignore_index= True).sort_values('POS').reset_index(drop=True)
#     ill = ill[["POS","REF","ALT","REF_DP","ALT_DP", "ALT_F"]]
    
#     if ill.empty == True:
#         new_row = pd.DataFrame(data= {'POS':0, 'REF':0, 'ALT':0, 'REF_DP':0, "ALT_DP":0, "ALT_F":0 }, index=[0])
#         ill = pd.concat([ill, new_row]).reset_index(drop=True)

#     try: 
#         for index, row in a_conc.iterrows():
#             if (row["plate_well"]  == dir.split('-')[0]):
#                 ill["lineage"] = row["lineage"]
#                 sample_lineage = ill.lineage.values.tolist()[0]
#     except IndexError:
#         pass
    
#     control_for_sample = All_control[All_control["lineage"].eq(sample_lineage)].drop_duplicates(ignore_index= True)

#     ill_depth = pd.merge(ill, ill_depth, how= "left", on=["POS"])

#     ill_depth_true = pd.merge(control_for_sample, ill_depth, on=["POS","REF","ALT","lineage"], how="outer", indicator=True)
#     FP_I = ill_depth_true[(ill_depth_true["_merge"] != "both") & (ill_depth_true["_merge"] != "left_only")]
#     FP_ill = FP_I.drop("_merge", axis=1)

    
#     FP_ont["source"] = dir
#     FP_ill["source"] = dir

#     for index, row in a_conc.iterrows():
#         if row["plate_well"] == dir.split('-')[0]:
#             alpha_conc = row["Alpha-C15"]
#             alpha_exp_frequency = row["alpha_exp_frequency"]
#             beta_conc = row["Beta-C16"]
#             beta_exp_frequency = row["beta_exp_frequency"]
#             delta_conc = row["Delta-C23"]
#             delta_exp_frequency = row["deltaC23_exp_frequency"]
#             deltaAY2_conc = row["DeltaAY.2-C29"]
#             deltaAY2_exp_frequency = row["deltaAY2_exp_frequency"]
#             omicron_conc = row["Omicron-C48"]
#             omicron_frequency = row["omicron_exp_frequency"]
    
#             # Create new columns in FP_ont table and assign values
#             FP_ont["alpha_conc"] = alpha_conc
#             FP_ont["beta_conc"] = beta_conc
#             FP_ont["delta_conc"] = delta_conc
#             FP_ont["deltaAY2_conc"] = deltaAY2_conc
#             FP_ont["omicron_conc"] = omicron_conc
#             FP_ont["alpha_exp_frequency"] = alpha_exp_frequency
#             FP_ont["beta_exp_frequency"] = beta_exp_frequency
#             FP_ont["delta_exp_frequency"] = delta_exp_frequency
#             FP_ont["deltaAY2_exp_frequency"] = deltaAY2_exp_frequency
#             FP_ont["omicron_frequency"] = omicron_frequency
    
#             # Fill the entire columns with the assigned values
#             FP_ont["alpha_conc"].fillna(alpha_conc, inplace=True)
#             FP_ont["beta_conc"].fillna(beta_conc, inplace=True)
#             FP_ont["delta_conc"].fillna(delta_conc, inplace=True)
#             FP_ont["deltaAY2_conc"].fillna(deltaAY2_conc, inplace=True)
#             FP_ont["omicron_conc"].fillna(omicron_conc, inplace=True)
#             FP_ont["alpha_exp_frequency"].fillna(alpha_exp_frequency, inplace=True)
#             FP_ont["beta_exp_frequency"].fillna(beta_exp_frequency, inplace=True)
#             FP_ont["delta_exp_frequency"].fillna(delta_exp_frequency, inplace=True)
#             FP_ont["deltaAY2_exp_frequency"].fillna(deltaAY2_exp_frequency, inplace=True)
#             FP_ont["omicron_frequency"].fillna(omicron_frequency, inplace=True)

#             # Create new columns in FP_ill table and assign values
#             FP_ill["alpha_conc"] = alpha_conc
#             FP_ill["beta_conc"] = beta_conc
#             FP_ill["delta_conc"] = delta_conc
#             FP_ill["deltaAY2_conc"] = deltaAY2_conc
#             FP_ill["omicron_conc"] = omicron_conc
#             FP_ill["alpha_exp_frequency"] = alpha_exp_frequency
#             FP_ill["beta_exp_frequency"] = beta_exp_frequency
#             FP_ill["delta_exp_frequency"] = delta_exp_frequency
#             FP_ill["deltaAY2_exp_frequency"] = deltaAY2_exp_frequency
#             FP_ill["omicron_frequency"] = omicron_frequency
    
#             # Fill the entire columns with the assigned values
#             FP_ill["alpha_conc"].fillna(alpha_conc, inplace=True)
#             FP_ill["beta_conc"].fillna(beta_conc, inplace=True)
#             FP_ill["delta_conc"].fillna(delta_conc, inplace=True)
#             FP_ill["deltaAY2_conc"].fillna(deltaAY2_conc, inplace=True)
#             FP_ill["omicron_conc"].fillna(omicron_conc, inplace=True)
#             FP_ill["alpha_exp_frequency"].fillna(alpha_exp_frequency, inplace=True)
#             FP_ill["beta_exp_frequency"].fillna(beta_exp_frequency, inplace=True)
#             FP_ill["delta_exp_frequency"].fillna(delta_exp_frequency, inplace=True)
#             FP_ill["deltaAY2_exp_frequency"].fillna(deltaAY2_exp_frequency, inplace=True)
#             FP_ill["omicron_frequency"].fillna(omicron_frequency, inplace=True)
    
#     FP_ont_table.append(FP_ont)
#     FP_ill_table.append(FP_ill)
    
# FP_ONT_table = pd.concat(FP_ont_table).reset_index(drop=True)
# FP_ILL_table = pd.concat(FP_ill_table).reset_index(drop=True)

# FP_table = pd.merge(FP_ONT_table, FP_ILL_table, on=["POS","REF","ALT","lineage","source","alpha_conc","beta_conc","delta_conc","deltaAY2_conc","omicron_conc","alpha_exp_frequency","beta_exp_frequency","delta_exp_frequency","deltaAY2_exp_frequency","omicron_frequency"], how="outer", suffixes=("_ont", "_ill"))

# FP_table.to_csv(f'{directory}/FP_table_NT001.txt',sep="\t",index = False)
# #FP_ILL_table.to_csv(f'{directory}/FP_ILL_table_NT002',sep="\t",index = False)



# import os
# import pandas as pd
# from pathlib import Path
# import glob

# directory = '/home/mbxha18/allcombinations/New_illumina/exeter/exeter_ill_vs_exeter_ont_NT002/'
# comparison = "exeter_ill_vs_exeter_ont_NT002.mincov1.mpileup2snp.varscan"
# all_subdirs = glob.glob(f"{directory}/*/")

# for subdir in all_subdirs:
#     dir = subdir.split('/')[7]
#     print(f"{dir}")

#     os.chdir(f'{directory}/{dir}')
#     files = glob.glob(f"{directory}/{dir}.variant.{comparison}.varscan.tsv") 
#     file = pd.read_csv(files[0], sep="\t")



