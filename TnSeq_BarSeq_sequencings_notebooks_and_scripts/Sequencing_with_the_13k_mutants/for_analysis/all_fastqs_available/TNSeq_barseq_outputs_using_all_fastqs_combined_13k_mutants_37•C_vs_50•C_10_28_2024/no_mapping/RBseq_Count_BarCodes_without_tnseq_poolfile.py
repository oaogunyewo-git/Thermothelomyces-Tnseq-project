import argparse
import json
import sys
import numpy as np
import pandas as pd
from datetime import datetime

# Assume extract_barcode is defined elsewhere or provide a simple placeholder
def extract_barcode(line):
    # Placeholder: Replace with logic for extracting a barcode from a line
    return line.strip()

# Define OffByOneList function at the top of your script
def OffByOneList(barcode):
    bases = ['A', 'T', 'C', 'G']
    barcode_variants = []
    for i, base in enumerate(barcode):
        for new_base in bases:
            if new_base != base:  # Skip if it's the same base
                new_barcode = barcode[:i] + new_base + barcode[i+1:]
                barcode_variants.append(new_barcode)
    return barcode_variants
# Test call to ensure function works
print(OffByOneList("ATCG"))
# Main function
def main(argv):
    # Initial setup and argument parsing
    timestamp = datetime.now().strftime('%Y%m%H%M')
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--metafile", dest="metafile", help="Metadata file for BarSeq runs. A tab-delimited file with columns titled Fastq,SampleName,BarSeqModel,OutputDir.", default="metafile.txt")
    parser.add_argument("-l", "--logFile", dest="logFile", help="File to write run log to. Default is Count_TIMESTAMP.log", default="Count_" + timestamp + ".log")
    parser.add_argument("-q", "--qual", dest="minQual", help="Minimum quality score for the barcode region for a read to counted. Default is 10", default=10, type=int)
    parser.add_argument("-b", "--matchBefore", dest="matchBefore", help="Number of bases before the barcode to match. Default is 6", default=6, type=int)
    parser.add_argument("-a", "--matchAfter", dest="matchAfter", help="Number of bases after the barcode to match. Default is 6", default=6, type=int)
    parser.add_argument("-Q", "--quietMode", action='store_true', dest="quietMode", help="Give fewer details about barcodes found during the run and in the log file. Summary stats will still be reported in summaryStats.txt", default=False)
    options = parser.parse_args()

    # Assume metaFrame and fileIndexes are initialized from metafile
    metaFrame = pd.read_csv(options.metafile, sep='\t')
    fileIndexes = metaFrame.index

    # Initialize sumStats dictionary
    sumStats = {
        'FileIndex': [], 'MostAbundantBarcode': [], 'CountsForMostAbundantBarcode': [], 
        'ErrorRatePercent': [], 'OneRead': [], 'TwoReads': [], 'ThreeReadsOrMore': [],
        'EstimatedPopulationSize': []
    }

    # Placeholder for printUpdate function
    def printUpdate(logFile, message):
        print(message)

    barcodeCounts = {}  # Define barcodeCounts in the main scope

    for fileIndex in fileIndexes:
        sumStats['FileIndex'].append(fileIndex)
        indexRow = metaFrame.loc[fileIndex]
        
        barcodeCounts.clear()  # Reset barcodeCounts for each fileIndex

        if indexRow['UsePrecounted']:
            # Load precounted barcodes
            statusUpdate = '  Loading previously counted reads for ' + indexRow['OutputDir']
            printUpdate(options.logFile, statusUpdate)
            fileToOpen = indexRow['OutputDir'] + "/countsFiles/" + str(fileIndex) + '.counts'
            try:
                with open(fileToOpen, 'r') as file:
                    barcodeCounts = json.load(file)
            except IOError:
                statusUpdate = " Could not read file:" + fileToOpen + " ...exiting."
                printUpdate(options.logFile, statusUpdate)
                sys.exit()
        else:
            # Count barcodes directly from the input file if not using precounted data
            input_file = indexRow['Fastq']  # Assuming 'Fastq' column has the file path
            statusUpdate = '  Counting barcodes in ' + input_file
            printUpdate(options.logFile, statusUpdate)
            
            try:
                with open(input_file, 'r') as file:
                    for line in file:
                        barcode = extract_barcode(line)
                        if barcode:
                            if barcode in barcodeCounts:
                                barcodeCounts[barcode] += 1
                            else:
                                barcodeCounts[barcode] = 1
            except IOError:
                statusUpdate = " Could not read file:" + input_file + " ...exiting."
                printUpdate(options.logFile, statusUpdate)
                sys.exit()

        # Create barcodeCountsFrame from barcodeCounts within the loop
        barcodeCountsFrame = pd.DataFrame.from_dict(barcodeCounts, orient='index', columns=[fileIndex], dtype=int)

        # Find the most abundant barcode and count
        mostAbundantBarcode = barcodeCountsFrame[fileIndex].idxmax()
        mostAbundantCounts = barcodeCountsFrame.loc[mostAbundantBarcode][fileIndex]

        sumStats['MostAbundantBarcode'].append(mostAbundantBarcode)
        sumStats['CountsForMostAbundantBarcode'].append(mostAbundantCounts)

        # Sequencing error rate calculation
        if mostAbundantCounts > 1000:
            errorList = OffByOneList(mostAbundantBarcode)  # Placeholder for OffByOneList logic
            errorsSeen = [barcode for barcode in errorList if barcode in barcodeCountsFrame.index]
            errorCounts = barcodeCountsFrame.loc[errorsSeen][fileIndex].sum()
            errorRate = max(.005, float(errorCounts) / (mostAbundantCounts + errorCounts)) * 100
        else:
            errorRate = 1.0

        sumStats['ErrorRatePercent'].append(errorRate)

        ones = (barcodeCountsFrame[fileIndex] == 1).sum()
        twos = (barcodeCountsFrame[fileIndex] == 2).sum()
        threes = len(barcodeCountsFrame) - ones - twos

        sumStats['OneRead'].append(ones)
        sumStats['TwoReads'].append(twos)
        sumStats['ThreeReadsOrMore'].append(threes)

        # Chao population size estimate
        Nch = ones**2 / (2 * twos) if twos > 0 else np.nan
        sumStats['EstimatedPopulationSize'].append(int(round(Nch, 2)) if not np.isnan(Nch) else np.nan)

        # Print status update for non-quiet mode
        if not options.quietMode:
            printUpdate(options.logFile, f"  Estimated sequencing error rate for barcodes: {errorRate}%")
            printUpdate(options.logFile, f"  Barcodes seen once: {ones}, twice: {twos}, three times or more: {threes}")
            printUpdate(options.logFile, f"  Chao estimate of population size: {Nch}")

    # Save the summary statistics to a file
    sumStatFrame = pd.DataFrame.from_dict(sumStats)
    outputDir = metaFrame.loc[0, 'OutputDir']  # Use output directory from metadata
    fileToSave = outputDir + "/fastqSummaryStats.txt"
    printUpdate(options.logFile, f"Saving summary statistics for fastqs to: {fileToSave}")
    sumStatFrame.to_csv(fileToSave, sep='\t', index=None)

if __name__ == "__main__":
    main(sys.argv[1:])
