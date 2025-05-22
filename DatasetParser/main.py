import csv
from os import listdir
from os.path import  isfile, join
import pandas as pd
from typing import Callable

# Path to the dataset folder
basePath: str = '.\\Dataset\\' 
# Path to the gene models
geneModelsPath: str = basePath + 'Gene Models'
# Path to the clinical data
clinicalPath: str = basePath + 'clinical.tsv'

# Gene-model:
geneModel = "# gene-model: GENCODE v36\n"

# Gene output file
geneOutputFile = 'lncRna.csv'
clinicalOutputFile = 'clinical.csv'

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

# Reading the clinical data
clinicalData = pd.read_csv(clinicalPath, delimiter="\t")

def filterGeneModel(path: str) -> tuple[pd.DataFrame, str]:
    with open(path) as f:
        header = f.readline()

        # Verifying if the first line has the correct gene model
        if header == geneModel:
            # Updating the header variable
            header = f.readline()

            # Reading the file into a pandas dataframe
            frame = pd.read_csv(path, delimiter="\t", skiprows=1)

            # Filtering the dataframe to only include lncRNA
            frame = frame[frame['gene_type'] == 'lncRNA']

            # Saving the filtered data to a new file
            outputPath = removeSegmentsFromPath(path) + geneOutputFile
            frame.to_csv(outputPath, index=False, header=True)
        return frame, header    

def filterClinicalData(path: str) -> pd.DataFrame:
    global clinicalData
    caseId = path.split('\\')[-2]

    # Filtering the dataframe to only include lncRNA
    clinicalData = clinicalData[clinicalData['cases.case_id'] == caseId]

    # Saving the filtered data to a new file
    outputPath = removeSegmentsFromPath(path) + clinicalOutputFile
    clinicalData.to_csv(outputPath, index=False, header=True)

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

forEachDataSetFolder(geneModelsPath, lambda folder: forFirstFileWithExtension(folder, readTsvFile))

# Print the headers
print("Headers:" + str(len(headers)))
for key, value in headers.items():
    print(key, value)

# Print the lengths
print("\nLengths:" + str(len(lengths)))
for length in lengths:
    print(length)