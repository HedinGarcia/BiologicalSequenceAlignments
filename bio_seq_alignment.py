# An implementation of the Needleman-Wunsch - Global Alignment Algorithm
import sys
import csv

def scoring_matrix(compSequenceA, compSequenceB):
    if (compSequenceA == compSequenceB): return 1  # Reward on match
    return -1  #Penalty for mismatch

def print_matrix(matrix):
    numRows = len(matrix)
    for i in range(0,numRows):
        print(matrix[i])
    return

def needleman_wunsch_matrix(sequence1, sequence2):
    numColumns = len(sequence1) + 2
    numRows = len(sequence2) + 2
    matrix = [[0 for x in range(numColumns)] for x in range(numRows)]  # Create grid for sequence
    d = -2  # Gap penalty

    seqIndex, initIndex= 0, 0
    for i in range(0, 2):
        for j in range(0, numColumns):
            if( i == 0 and j >= 2):
                matrix[i][j] = sequence1[seqIndex]  # Fill the first sequence on top of the matrix starting from the third column
                seqIndex += 1
            if( i == 1 and j < numColumns-1):
                matrix[i][j + 1] = d * j  # Initialization Step
                initIndex += 1

    seqIndex, initIndex = 0, 0
    for j in range(0,2):
        for i in range(0, numRows):
            if( j == 0 and i >= 2):
                matrix[i][j] = sequence2[seqIndex]   # Fill the second sequence on left side of the matrix starting from the third row
                seqIndex += 1
            if( j == 1 and i < numRows-1):
                matrix[i + 1][j] = d * i  # Initialization Step
                initIndex += 1
    iIndex = 2
    jIndex = 2
    for i in range(iIndex,numRows):
        for j in range(jIndex,numColumns):
            matrix[i][j] = max(matrix[i-1][j-1]+scoring_matrix(matrix[i][0],matrix[0][j]), matrix[i][j-1] + d, matrix[i-1][j] + d)
    return matrix

# Read input CSV file from command line
paramList = sys.argv
inputfile = paramList[1]
if ".csv" in inputfile:
    print("The input CSV file is", inputfile)
    print("File content:")
    with open(inputfile, newline='', encoding='utf-8-sig') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        next(reader)  # skip header row
        for row in reader:
            print("Evaluating: " + row[0] + "," + row[1])
            print("Number of rows: " + str(len(row[1]) + 2) + ", " + "Number of columns: " + str(len(row[0]) + 2))
            print_matrix(needleman_wunsch_matrix(row[0], row[1]))
else: print("Input is not a .csv file")
