import pandas as pd

# NOTE: Make sure you run the project from the DatesetParser folder

# Path to the dataset folder
basePath: str = './Dataset/' 
# Path to settings
settingsPath: str = basePath + 'Settings/'
# Path to gdc data files
inputPath: str = basePath + 'OriginalFiles/'
# Path to the parsed files
outputPath: str = basePath + 'ProcessedFiles/'
# Path to the metadata file
metadataPath: str = settingsPath + 'metadata.cart.2025-05-22.json'

# Gene-model:
geneModel = "# gene-model: GENCODE v36\n"

def updateFileExtension(path: str, newExtension: str) -> str:
    # Splitting the path into segments
    segments = str.split(path, '.')
    # Removing the last segment (the current extension)
    segments = segments[0:-1]
    # Adding the new extension
    segments.append(newExtension)
    # Joining the remaining segments back together
    return '.'.join(segments)

def processGeneData(originalPath: str) -> tuple[pd.DataFrame | None, str]:
    with open(originalPath) as file:
        header = file.readline()

        # Verifying if the first line has the correct gene model
        if header == geneModel:
            # Updating the header variable
            header = file.readline()

            # Reading the file into a pandas dataframe
            dataframe = pd.read_csv(originalPath, delimiter="\t", skiprows=1)

            # Filtering the dataframe to only include lncRNA
            dataframe = dataframe[dataframe['gene_type'] == 'lncRNA']

            # Creating the columns for unstranded and tpm values
            dataframe['unstranded_col'] = dataframe['gene_name'] + '_unstranded'
            dataframe['tpm_col'] = dataframe['gene_name'] + '_tpm_unstranded'

            # Build a dictionary with new column names and their values
            unstranded_vals = dict(zip(dataframe['unstranded_col'], dataframe['unstranded']))
            tpm_vals = dict(zip(dataframe['tpm_col'], dataframe['tpm_unstranded']))

            # Merge the two into one row using a dictionary
            combined = {**unstranded_vals, **tpm_vals}

            # Create a new single-row DataFrame from the combined dict
            result_df = pd.DataFrame([combined])

        return result_df, header    

def processClinicalData(originalPath: str) -> pd.DataFrame | None:
    # Loading the clinical data from the original file
    dataFrame = pd.read_xml(originalPath)

    # Flattening the XML data
    firstRow = dataFrame.loc[0]
    for item in dataFrame.columns:
        # Checking if the column is empty and removing it
        if firstRow[item] is None or firstRow[item] == '':
            dataFrame.loc[0, item] = dataFrame.loc[1, item]  # Copying the value from the second row

    # Filtering the clinical data to only include the relevant columns
    columnsToKeep = ['histological_type', 'icd_o_3_histology']

    # Checking if the columns exist in the dataframe
    if not set(columnsToKeep).issubset(dataFrame.columns):
        print("Some columns are missing in the clinical data. Please check the file.")
        return None
    
    # The xml file is not parsed flat but in two rows, so we need to flatten it
    dataFrame.drop(index=1, inplace=True)
  
    return dataFrame

# Reading the metadata to create folders for each case
def readMetadataFile(path: str) -> dict:
    import json
    with open(path) as f:
        d = json.load(f)
        return d
     
def mergeCaseData(metadataPath: str, inputPath: str, outputPath: str, storeSubfiles:bool = True) -> pd.DataFrame:
    import os

    # Initializing the two main data frames
    dataFrameColumns = ['case_id']
    clinicalDataFrame = pd.DataFrame(columns=dataFrameColumns)
    geneDataFrame = pd.DataFrame(columns=dataFrameColumns)

    fileCount = 0

    data = readMetadataFile(metadataPath)
    for file in data:
        # Creating the output folder for the case
        caseId = file['associated_entities'][0]['case_id']

        # Retrieving the file name and folder name
        fileName = file['file_name']
        folderName = file['file_id']
        fileFormat = file['data_format']

        fileCount += 1
        print(f"Processing file {fileCount}: {fileName} for case {caseId}")

        # Creating the storage and output paths
        dataFile = inputPath + folderName + '/' + fileName
        outputFile = updateFileExtension(outputPath + caseId + '/' + fileName, "csv")
  
        # Checking if the file exists (Data set may be corrupted or incomplete)
        if not os.path.isfile(dataFile):
            print("File not found: " + dataFile)
            continue

         # Handling the different file types
        dataFrame = None
        if fileFormat == 'TSV':
            # continue
            (dataFrame, _) = processGeneData(dataFile) 
        else:  
            dataFrame = processClinicalData(dataFile)

        # Adding the data to the main dataframe and storing the file
        if dataFrame is not None and not dataFrame.empty:
            dataFrame['case_id'] = caseId
         
            # Adding the loaded data frame to the main data frame
            if fileFormat == 'TSV':
                geneDataFrame = pd.concat([geneDataFrame, dataFrame])
            else:  
                clinicalDataFrame = pd.concat([clinicalDataFrame, dataFrame])
             
            # Only storing the subfiles if the flag is set
            if storeSubfiles:
                dataFrame.to_csv(outputFile, index=False, header=True)


    return pd.merge(clinicalDataFrame, geneDataFrame, on="case_id", how="inner")      

mainDataFrame = mergeCaseData(metadataPath, inputPath, outputPath, True)
print(mainDataFrame)

# Storing the main data frame to a file
mainDataFrame.to_csv(outputPath + 'merged_data.csv', index=False, header=True)