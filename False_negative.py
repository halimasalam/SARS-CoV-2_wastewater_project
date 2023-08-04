# .ill.mpileup2snp.mincov4.minvar0.06.nostrandbiasfilter.varscan.tsv
# .ont.mpileup.mincov4.minvar0.06.varscan.tsv

#
# .ill.mincov1.mpileup2snp.nostrandbiasfilter.varscan.tsv
# .ont.mincov1.mpileup2snp.nostrandbiasfilter.varscan.tsv
# import required libraries
import os
import pandas as pd
from pathlib import Path
import glob
import numpy as np
import re

# Set the directory path that contains each samples folder 
directory = '/home/mbxha18/allcombinations/New_illumina/exeter/exeter_ill_vs_exeter_ont_NT002/'
table = "FN_all_table_NT002" # final table name
# Read data from the scontrol file
All_control = pd.read_csv(f'/home/mbxha18/OvsI_prom/ONT_match_prom/All_control_with_indels.tsv',sep="\t")
# remove indel positions
Control_remove_indel_positions = All_control[~All_control['POS'].str.startswith(('del', 'ins'))].astype({'POS': int})
# remove misalignment_errors_positions 
misalignment_errors_positions = [21766, 22034, 6515, 21987, 21988, 22196, 28363]
Control = Control_remove_indel_positions[~Control_remove_indel_positions['POS'].isin(misalignment_errors_positions)]
# Read data from the control file without covmix
lineage_control = pd.read_csv(f'/home/mbxha18/OvsI_prom/ONT_match_prom/single_lineage_control.txt',sep="\t")

# Read data from the TWIST_synthetic_mixes_sample_sheet that contains the sample information
a_conc = pd.read_csv("/home/mbxha18/allcombinations/controls_old/TWIST_synthetic_mixes_sample_sheet_150222_edited.txt", sep="\t",index_col=0)
# Initialize an empty list to store data for each sample
FN_table=[]

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
        elif file.endswith(".ont.mincov1.mpileup2snp.varscan.tsv"):
            ont_file_present = True
    
    if ill_file_present and ont_file_present:
        print("Run")
        # Set the working directory to the current subdirectory
        os.chdir(f'{directory}/{dir}')
        sample_ill = f"{dir}_ill"
        sample_ont = f"{dir}_ont"

        # Get a list of needed files and read them in dfs
        ont_files = glob.glob(f"{directory}/{dir}/*.ont.mincov1.mpileup2snp.nostrandbiasfilter.varscan.tsv") 
        ont = pd.read_csv(ont_files[0], sep="\t")
        ont_pileup = glob.glob(f"{directory}/{dir}/*.ont.pileup")
        av_cov_ont_files = glob.glob(f"{directory}/{dir}/*.ont.mosdepth.summary.txt")
        av_cov_ont = pd.read_csv(av_cov_ont_files[0], sep="\t")
        ont_depth_files = glob.glob(f"{directory}/{dir}/*.ont.samtools.cov")
        ont_depth = pd.read_csv(ont_depth_files[0], sep="\t", header = None, names = ["REGION","POS","DEPTH_ont"], usecols = ["POS","DEPTH_ont"])

        ill_files = glob.glob(f"{directory}/{dir}/*.ill.mincov1.mpileup2snp.varscan.tsv") 
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

            # merge the varscan df with the control for this sample and depth file
            ont_true = pd.merge(control_for_sample, ont, how= "left", on=["POS","REF","ALT","lineage"])
            ont_true_depth = pd.merge(ont_true, ont_depth, how= "left", on=["POS"])
            if "DEPTH_ont_x" in ont_true_depth.columns: 
                ont_true_depth["DEPTH_ont"] = ont_true_depth["DEPTH_ont_x"].combine_first(ont_true_depth["DEPTH_ont_y"])
                ont_true_depth.drop("DEPTH_ont_x",axis=1, inplace=True)
                ont_true_depth.drop("DEPTH_ont_y",axis=1, inplace=True)
            ont_table = ont_true_depth.fillna(0)
            A = ont_table

                 
 

            
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
                
            ill_true = pd.merge(control_for_sample, ill, how= "left", on=["POS","REF","ALT","lineage"])
            ill_true_depth = pd.merge(ill_true, ill_depth, how= "left", on=["POS"])
            if "DEPTH_ill_x" in ill_true_depth.columns: 
                ill_true_depth["DEPTH_ill"] = ill_true_depth["DEPTH_ill_x"].combine_first(ill_true_depth["DEPTH_ill_y"])
                ill_true_depth.drop("DEPTH_ill_x",axis=1, inplace=True)
                ill_true_depth.drop("DEPTH_ill_y",axis=1, inplace=True)
            ill_table = ill_true_depth.fillna(0)
            B = ill_table

 

            All = pd.merge(A, B, how="outer", on=["POS","REF","ALT","lineage"], suffixes=("_ont", "_ill"))
            
            def categorise_ill(row): 
                if row['lineage'] == "neg":
                    return 'neg'
                if row['ALT_DP_ill'] == 0 and row["DEPTH_ill"] <= 0: 
                    return 'no_coverage'
                elif row['ALT_DP_ill']== 0 and row["DEPTH_ill"] >= 1 :
                    return 'not_called'
                else: 
                    return 'called'
            All["REASON_ill"] = All.apply(lambda row: categorise_ill(row), axis=1)
                    
            def categorise_ont(row): 
                if row['lineage'] == "neg":
                    return 'neg'
                if row['ALT_DP_ont'] == 0 and row["DEPTH_ont"] <= 0: 
                    return 'no_coverage'
                elif row['ALT_DP_ont']== 0 and row["DEPTH_ont"] >= 1 :
                    return 'not_called'
                else: 
                    return 'called'
            All["REASON_ont"] = All.apply(lambda row: categorise_ont(row), axis=1)
   
            # extract ont FN (true snps not called)
            FN = (All[All["REASON_ont"].eq("not_called")]) 
            # create new columnin the FN table to contain the directory name(sample name)
            FN["source"] = dir
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
                    FN["alpha_conc"] = alpha_conc
                    FN["beta_conc"] = beta_conc
                    FN["delta_conc"] = delta_conc
                    FN["deltaAY2_conc"] = deltaAY2_conc
                    FN["omicron_conc"] = omicron_conc
                    FN["alpha_exp_frequency"] = alpha_exp_frequency
                    FN["beta_exp_frequency"] = beta_exp_frequency
                    FN["delta_exp_frequency"] = delta_exp_frequency
                    FN["deltaAY2_exp_frequency"] = deltaAY2_exp_frequency
                    FN["omicron_frequency"] = omicron_frequency
            
                    FN["alpha_conc"].fillna(alpha_conc, inplace=True)
                    FN["beta_conc"].fillna(beta_conc, inplace=True)
                    FN["delta_conc"].fillna(delta_conc, inplace=True)
                    FN["deltaAY2_conc"].fillna(deltaAY2_conc, inplace=True)
                    FN["omicron_conc"].fillna(omicron_conc, inplace=True)
                    FN["alpha_exp_frequency"].fillna(alpha_exp_frequency, inplace=True)
                    FN["beta_exp_frequency"].fillna(beta_exp_frequency, inplace=True)
                    FN["delta_exp_frequency"].fillna(delta_exp_frequency, inplace=True)
                    FN["deltaAY2_exp_frequency"].fillna(deltaAY2_exp_frequency, inplace=True)
                    FN["omicron_frequency"].fillna(omicron_frequency, inplace=True)
            
            FN_table.append(FN) # Append the df to the FN_table list
    else:
        pass ## skip to the next iteration (sample) if the files are not present in the folder
    
FN_all_table = pd.concat(FN_table).reset_index(drop=True) # Concatenate all DataFrames in the FN_table list into a single DataFrame 

# To assign the snps to their respective lineage, Merge the FN_all_table df to the lineage_control based on the position column
merged_df = pd.merge(FN_all_table, lineage_control, on="POS", how="left")  
# Group the merged DataFrame by position and aggregate the lineages into a list
grouped_df = merged_df.groupby("POS")["source_lineage"].agg(lambda x: ', '.join(x))    
# Create a new column "lineage" in FN_all_table df and fill it with the aggregated lineages
FN_all_table["source_lineage"] = FN_all_table["POS"].map(grouped_df)
# Iterate over each value in the "lineage" column
for i, cell_data in enumerate(FN_all_table["source_lineage"]):
    values = cell_data.split(", ")
    unique_values = list(set(values))
    result = ", ".join(unique_values)
    FN_all_table.at[i, "source_lineage"] = result

# Save the final DataFrame to a TSV file
FN_all_table.to_csv(f'{directory}/{table}.txt',sep="\t",index = False)

