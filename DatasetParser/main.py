from os import listdir
from os.path import  isfile, join
import pandas as pd
from typing import Callable

# NOTE: Make sure you run the project from the DatesetParser folder

# Path to the dataset folder
basePath: str = './Dataset/' 
# Path to settings
settingsPath: str = basePath + 'Settings/'
# Path to the gene models
geneModelsPath: str = basePath + 'Gene Models'
# Path to the clinical data
clinicalPath: str = basePath + 'clinical.tsv'
# Path to gdc data files
inputPath: str = basePath + 'OriginalFiles/'
# Path to the parsed files
outputPath: str = basePath + 'ProcessedFiles/'
# Path to the metadata file
metadataPath: str = settingsPath + 'metadata.cart.2025-05-22.json'

# Gene-model:
geneModel = "# gene-model: GENCODE v36\n"

# Tracking the different heafers and file lengths
headers: dict = dict()
lengths: set = set()

# Running code for each sub folder in a given folder
def forEachDataSetFolder(path: str, func: Callable[[str], None]) -> None:
    # Looping through each file in the directory
    for file in listdir(path):
        folder = join(path, file)
        # If the path is not a file, it is a folder
        if not isfile(folder):       
            func(folder)

# Running code for the first file in each sub folder
def forFirstFileWithExtension(path: str, func: Callable[[str], None], extension = '.tsv') -> None:
    # Looping through each file in the directory
    for file in listdir(path):
        filePath = join(path, file)
        # If the path is a file trigger the function and break
        if isfile(filePath) and filePath.endswith(extension):
            func(filePath)
            break  

def removeSegmentsFromPath(path: str, count: int = 1, hasTrailingSlash: bool = True) -> str:
    # Splitting the path into segments
    segments = str.split(path, '\\')
    # Removing the last count segments
    segments = segments[0:-count]
    # Joining the remaining segments back together
    path = '\\'.join(segments)
    # Adding a trailing slash if needed
    return path + ('\\' if hasTrailingSlash else '')

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
            frame = pd.read_csv(originalPath, delimiter="\t", skiprows=1)

            # Filtering the dataframe to only include lncRNA
            frame = frame[frame['gene_type'] == 'lncRNA']

            flattendFrame = {}
            for _, row in frame.iterrows():
                name = row['gene_name']
                flattendFrame[name + '_unstranded'] = row['unstranded']
                flattendFrame[name + '_tpm_unstranded'] = row['tpm_unstranded']   

            frame = pd.DataFrame(flattendFrame, index=[0, len(flattendFrame)])
        return frame, header    

def processClinicalData(originalPath: str) -> pd.DataFrame | None:
    # Loading the clinical data from the original file
    clinicalData = pd.read_xml(originalPath)

    # Filtering the clinical data to only include the relevant columns
    columnsToKeep = ['histological_type', 'icd_o_3_histology']

    # Checking if the columns exist in the dataframe
    if not set(columnsToKeep).issubset(clinicalData.columns):
        print("Some columns are missing in the clinical data. Please check the file.")
        return None
    
    # Removing e nb mpty row, that happens to be added when reading the xml file
    clinicalData.drop(index=1, inplace=True)
  
    return clinicalData

# Reading the metadata to create folders for each case
def readMetadataFile(path: str) -> dict:
    import json
    with open(path) as f:
        d = json.load(f)
        return d
     
def mergeCaseData(metadataPath: str, mainDataFrame: pd.DataFrame, storeSubfiles:bool = True) -> pd.DataFrame:
    global outputPath, inputPath
    import os
    data = readMetadataFile(metadataPath)
    for file in data:
        # Creating the output folder for the case
        caseId = file['associated_entities'][0]['case_id']

        # Retrieving the file name and folder name
        fileName = file['file_name']
        folderName = file['file_id']
        fileFormat = file['data_format']

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
            dataFrame['case_id'] = pd.Series([caseId])
            if 'case_id' in mainDataFrame.columns and not mainDataFrame[mainDataFrame['case_id'] == caseId].empty:
                pd.merge(dataFrame, mainDataFrame, on="case_id")  
            else: 
                mainDataFrame = pd.concat([mainDataFrame, dataFrame], ignore_index=True)      
            # Only storing the subfiles if the flag is set
            if storeSubfiles:
                dataFrame.to_csv(outputFile, index=False, header=True)

        print("Processed: " + dataFile + " -> " + outputFile)    

        # print(caseId)    
    return mainDataFrame     

mainDataFrame = mergeCaseData(metadataPath, pd.DataFrame())
print(mainDataFrame)

# Storing the main data frame to a file
mainDataFrame.to_csv(outputPath + 'merged_data.csv', index=False, header=True)