import csv
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

def filterGeneModel(originalPath: str, storePath: str) -> tuple[pd.DataFrame, str]:
    with open(originalPath) as f:
        header = f.readline()

        # Verifying if the first line has the correct gene model
        if header == geneModel:
            # Updating the header variable
            header = f.readline()

            # Reading the file into a pandas dataframe
            frame = pd.read_csv(originalPath, delimiter="\t", skiprows=1)

            # Filtering the dataframe to only include lncRNA
            frame = frame[frame['gene_type'] == 'lncRNA']

            # Filtering the clinical data to only include the relevant columns
            columnsToKeep = ['gene_id', 'gene_name', 'unstranded', 'tpm_unstranded']

            flattendFrame = {}
            for _, row in frame.iterrows():
                name = row['gene_name']
                flattendFrame[name + '_unstranded'] = row['unstranded']
                flattendFrame[name + '_tpm_unstranded'] = row['tpm_unstranded']   
            # Storing the filterd file
            pd.DataFrame(flattendFrame, index=[0, len(flattendFrame)]).to_csv(storePath, index=False, header=True)
        return frame, header    

def filterClinicalData(originalPath: str, storePath: str, caseId: str) -> pd.DataFrame:
    # Loading the clinical data from the original file
    clinicalData = pd.read_xml(originalPath)

    # Filtering the clinical data to only include the relevant columns
    columnsToKeep = ['histological_type', 'icd_o_3_histology']

    # Checking if the columns exist in the dataframe
    if not set(columnsToKeep).issubset(clinicalData.columns):
        print("Some columns are missing in the clinical data. Please check the file.")
        return pd.DataFrame()

    clinicalData['case_id'] = caseId
    clinicalData.to_csv(storePath, index=False, header=True)

    return clinicalData

def readTsvFile(filePath: str) -> None:
    (frame, header) = filterGeneModel(filePath)
    clinicalFrame = filterClinicalData(filePath)

    # if header in headers.keys():
    #     headers[header] += 1
    # else:
    #     headers[header] = 1

    # length = sum(1 for _ in open(filePath))
    # lengths.add(length)

# Reading the metadata to create folders for each case
def readMetadataFile(path: str) -> dict:
    import json
    with open(path) as f:
        d = json.load(f)
        return d
     
def mergeCaseData(metadataPath: str, mainDataFrame: pd.DataFrame) -> pd.DataFrame:
    global outputPath, inputPath
    import os
    # import shutil

    data = readMetadataFile(metadataPath)
    for file in data:
        # Creating the output folder for the case
        caseId = file['associated_entities'][0]['case_id']
        os.makedirs(outputPath + caseId, exist_ok=True)

        # Getting the file name and folder name
        fileName = file['file_name']
        folderName = file['file_id']
        fileFormat = file['data_format']
        dataFile = inputPath + folderName + '/' + fileName
        outputFile = updateFileExtension(outputPath + caseId + '/' + fileName, "csv")
  
        # Checking if the file exists (Data set may be corrupted or incomplete)
        if not os.path.isfile(dataFile):
            print("File not found: " + dataFile)
            continue

         # Handling the different file types
        dataFrame = None
        if fileFormat == 'TSV':
            continue
            filterGeneModel(dataFile, outputFile) 
        else:  
            dataFrame = filterClinicalData(dataFile, outputFile, caseId)

        # Adding the data to the main dataframe
        if dataFrame is not None and not dataFrame.empty:
            mainDataFrame = pd.concat([mainDataFrame, dataFrame], ignore_index=True)    

        print("Processed: " + dataFile + " -> " + outputFile)    

        # print(caseId)    
    return mainDataFrame     

mainDataFrame = mergeCaseData(metadataPath, pd.DataFrame())

# Storing the main data frame to a file
mainDataFrame.to_csv(outputPath + 'merged_data.csv', index=False, header=True)