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

def needleman_wunsch_sequence_alignment(sequence1, sequence2):
    numColumns = len(sequence1) + 2
    numRows = len(sequence2) + 2
    matrix = [[0 for x in range(numColumns)] for x in range(numRows)]  # Create grid for sequence
    backtrackingMatrix = [[0 for x in range(numColumns)] for x in range(numRows)]  # Tracking chosen paths
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
    # Calculate score of sequences
    iIndex = 2
    jIndex = 2
    for i in range(iIndex,numRows):
        for j in range(jIndex,numColumns):
            diagonalPath = matrix[i-1][j-1] + scoring_matrix(matrix[i][0],matrix[0][j])
            leftPath, topPath =  matrix[i][j-1] + d, matrix[i-1][j] + d
            highestScore = max(diagonalPath, topPath, leftPath)
            matrix[i][j] = highestScore
            # Track the path chosen for score of pairing
            if (highestScore == diagonalPath and highestScore == leftPath):
                backtrackingMatrix[i][j] = "L"  # Priority over this option
            elif (highestScore == diagonalPath and highestScore == topPath):
                backtrackingMatrix[i][j] = "T"
            elif (highestScore == leftPath and highestScore == topPath):
                backtrackingMatrix[i][j] = "L"
            elif (highestScore == diagonalPath):
                backtrackingMatrix[i][j] = "D"
            elif (highestScore == leftPath):
                backtrackingMatrix[i][j] = "L"
            elif (highestScore == topPath):
                backtrackingMatrix[i][j] = "T"
    # Backtrack to align sequences
    lastRow = numRows -1
    lastColumn = numColumns -1
    score = matrix[lastRow][lastColumn]
    while(lastRow !=1 or lastColumn !=1):
        path = backtrackingMatrix[lastRow][lastColumn]
        if (path == "D"):
            lastRow -= 1
            lastColumn -= 1
        elif (path == "T"):
            lastRow -= 1
            sequence1 = sequence1[:lastColumn-1] + "-" + sequence1[lastColumn-1:]
        elif (path == "L"):
            lastColumn -= 1
            sequence2 = sequence2[:lastRow-1] + "-" + sequence2[lastRow-1:]
        elif (lastRow == 1 and lastColumn != 1):  # When the path is in the init values
            lastColumn -= 1
            sequence2 = "-" + sequence2
        elif (lastColumn == 1 and lastRow != 1):  # Likewise when the path is in the init values
            lastRow -= 1
            sequence1 = "-" + sequence1
    return matrix, backtrackingMatrix, sequence1, sequence2, score

# Read input CSV file from command line
paramList = sys.argv
inputfile = paramList[1]
if ".csv" in inputfile:
    inputfile = sys.argv[1]
    with open(inputfile, newline='', encoding='utf-8-sig') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        next(reader)  # skip header row
        for row in reader:
            print("Evaluating: " + row[0] + "," + row[1])
            print("Number of rows: " + str(len(row[1]) + 2) + ", " + "Number of columns: " + str(len(row[0]) + 2))
            matrix, backtrackingMatrix, seq1, seq2, score = needleman_wunsch_sequence_alignment(row[0], row[1])
            print_matrix(backtrackingMatrix)
            print("\n")
            print_matrix(matrix)
            print("Result: " + seq1 + " " + seq2 + " " + str(score) + "\n")
else: print("Input is not a .csv file")
