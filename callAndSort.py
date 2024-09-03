import argparse
import gzip
import time
import os

scriptVersion = "Bug Fix Post LARGE-PD Annual Meeting"

#=========================================================== Utils ====================================================
def createFolder(folderName, logFile):
    try:
        os.mkdir(folderName)
    except OSError as error:
        print(error)
    logFile.write(f"Folder {folderName} created using os library")

    return folderName

def execute(commandLine, logFile, execute = True):
    start = time.time()
    print(commandLine)

    if execute:
        os.system(commandLine)
    end = time.time()
    diff = end-start

    logFile.write(f'{commandLine}\nTime: {diff} seconds\n\n')


#=========================================================== QC ====================================================
def gtc2VCF(bcftools, bpm, egt, csv, folderGTC, allGTCFolder, genomeReference, vcfFolder, threads, outputName, mem, logFile):
    for file in os.listdir(folderGTC):
        line = f"cp {folderGTC}/{file} {allGTCFolder}"
        #execute(line, logFile)

    line = f"{bcftools} +gtc2vcf --no-version -Oz --bpm {bpm} --egt {egt} --csv {csv} --gtcs {allGTCFolder} " \
           f"--fasta-ref {genomeReference} --output {vcfFolder}/{outputName}.vcf.gz --threads {threads} " \
           f"--use-gtc-sample-names --adjust-clusters"
    #line = f"{bcftools} +gtc2vcf --no-version -Oz --bpm {bpm} --egt {egt} --csv {csv} --gtcs {allGTCFolder} " \
    #       f"--output {vcfFolder}/{outputName}.vcf.gz --threads {threads} " \
    #       f"--use-gtc-sample-names --adjust-clusters"

    execute(line, logFile)

    line = f"{bcftools} sort {vcfFolder}/{outputName}.vcf.gz -Oz -o {vcfFolder}/{outputName}_Sort.vcf.gz -m {mem}"
    execute(line, logFile)

    line = f"{bcftools} index {vcfFolder}/{outputName}_Sort.vcf.gz"
    execute(line, logFile)

    return f"{vcfFolder}/{outputName}_Sort.vcf.gz"



def generateGTC(iaap, bpm, egt, folder, outFolder, outputName, batchList, threads, logFile):
    threadsTest = 8

    inputFile = open(batchList)
    print(f"Creating the file {folder}/{outputName}_SampleSheet.csv")
    outputFile = open(f"{folder}/{outputName}_SampleSheet.csv", 'w')
    outputFile.write("[Data]\n")
    outputFile.write("Sample_ID,SentrixBarcode_A,SentrixPosition_A,Path\n")

    dictBatch = {}
    listID = []

    for line in inputFile:
        illuminaCode2External, sampleID2Sentrix, path, ID = line.strip().split()

        dictCodes = {}

        fileCodes = open(illuminaCode2External)
        header = True
        for line in fileCodes:
            if header:
                if "IlluminaID" in line:
                    header = False
            else:
                illuminaID, externalID = line.strip().replace(" ", "_").split()

                if externalID == "#N/A":
                    externalID = "NonLARGE"

                if externalID not in listID:
                    listID.append(externalID)
                else:
                    count = 1
                    baseID = externalID
                    externalID = f"{baseID}_{count}"

                    while externalID in listID:
                        count = count + 1
                        externalID = f"{baseID}_{count}"

                    listID.append(externalID)


                if illuminaID not in dictCodes:
                    dictCodes[illuminaID] = {}
                    dictCodes[illuminaID]["ExternalID"] = externalID
                    dictCodes[illuminaID]["Batch"] = ID
        fileCodes.close()

        fileSentrix = open(sampleID2Sentrix)
        header = True
        for line in fileSentrix:
            if header:
                if "Sample_ID" in line:
                    header = False
            else:
                illuminaID, sentrixBarcode, sentrixPosition = line.strip().split(",")[0:3]
                if illuminaID in dictCodes:
                    dictCodes[illuminaID]['sentrixBarcode'] = sentrixBarcode
                    dictCodes[illuminaID]['sentrixPosition'] = sentrixPosition
                else:
                    input(f"Error: the ID {illuminaID} was not in external code")

        for illuminaID in dictCodes:
            sentrixBarcode = dictCodes[illuminaID]['sentrixBarcode']
            sentrixPosition = dictCodes[illuminaID]['sentrixPosition']
            pathToLook = path + "/" + sentrixBarcode
            externalID = dictCodes[illuminaID]["ExternalID"] + "_" + dictCodes[illuminaID]["Batch"]
            dictBatch[f"{sentrixBarcode}_{sentrixPosition}"] = externalID

            if os.path.exists(pathToLook):
                outputFile.write(f"{externalID},{sentrixBarcode},{sentrixPosition},{pathToLook}\n")
            else:
                input(f"Not found: {pathToLook}")

        fileSentrix.close()

    outputFile.close()

    line = f"{iaap} gencall {bpm} {egt} -s {folder}/{outputName}_SampleSheet.csv --output-gtc {outFolder} " \
           f"-t {threadsTest} -c 0.2"
           #f"--gender-estimate-call-rate-threshold -0.1 -t {threadsTest} -c 0.2"
    execute(line, logFile)

    return f"{folder}/{outputName}_SampleSheet.csv"

def removeLiftOverProblems(vcfFileName, errors, bcftools, outputName, outputFolder):
    dictID = {}
    fileError = open(errors)
    for line in fileError:
        split = line.split()
        dictID[split[0]] = ""

    vcfFile = gzip.open(vcfFileName)
    vcfOut = open(f"{outputFolder}/{outputName}_NoLiftProblems.vcf", "w")

    count = 0
    removed = 0
    header = True
    for line in vcfFile:
        line = line.decode("utf-8")
        if header:
            vcfOut.write(line)
            if "#CHROM" in line:
                header = False
        else:
            count = count + 1
            split = line.split()

            if split[2] not in dictID:
                vcfOut.write(f"{line}")
            else:
                removed = removed + 1
    vcfOut.close()
    print(f"The VCF has {count} variants. We removed {removed} from this")
    #os.system(f"{bcftools} view -T ^{outputFolder}/errorToRemoveLift.txt -Oz -o {outputFolder}/{outputName}_NoLiftProblems.vcf.gz {vcfFileName}")
    os.system(f"bgzip {outputFolder}/{outputName}_NoLiftProblems.vcf")
    os.system(f"{bcftools} index {outputFolder}/{outputName}_NoLiftProblems.vcf.gz")

    return f"{outputFolder}/{outputName}_NoLiftProblems.vcf.gz"

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Call and sort alleles')

    requiredGeneral = parser.add_argument_group("Required arguments for all steps")
    requiredGeneral.add_argument('-o', '--outputFolder', help='Name of output folder', required=True)
    requiredGeneral.add_argument('-O', '--outputName', help='Name of output name', required=True)

    requiredCalling = parser.add_argument_group("Required arguments for variant calling")
    requiredCalling.add_argument('-b', '--bpm', help='BPM File from Illumina', required=True)
    requiredCalling.add_argument('-c', '--csv', help='CSV manifest file', required=True)
    requiredCalling.add_argument('-e', '--egt', help='EGT File from Illumina', required=True)
    requiredCalling.add_argument('-r', '--genomeReference', help='Human genome reference', required=True)
    requiredCalling.add_argument('-l', '--listError', help='List of problematic variants', required=True)
    requiredCalling.add_argument('-L', '--batchList',
                                 help='File with four columns: \"IlluminaID_ExternalID, SampleID_SentrixBarcode_SentrixPosition,'
                                      'Path to RawFiles,ID\" separated by tab', required=True)

    programs = parser.add_argument_group("Programs")
    programs.add_argument('-B', '--bcftools', help='bcftools', required=False, default="bcftools")
    programs.add_argument('-I', '--iaap', help='Illumina Array Analysis Platform path ', required=True)


    optional = parser.add_argument_group("Optional arguments")
    optional.add_argument('-t', '--threads', help='Number of threads', default=96, type=int, required=False)
    optional.add_argument('-m', '--mem', help='Memmory to bcftools sort', default="75G", required=False)
    optional.add_argument('-a', '--analysis', help='Generate the analysis data', default=False, action="store_true")
    args = parser.parse_args()

    folder = args.outputFolder
    os.system(f"mkdir {folder}")
    logFile = open(f"{args.outputFolder}/{args.outputName}.log", 'w')

    folderGTC = createFolder(folder+"/GTC/", logFile)
    sampleSheet = generateGTC(args.iaap, args.bpm, args.egt, folder, folderGTC, args.outputName, args.batchList, args.threads, logFile)

    allGTCFolder = createFolder(folder+"/GTCAll/", logFile)
    vcfFolder = createFolder(folder+"/VCFs/", logFile)
    VCF = gtc2VCF(args.bcftools, args.bpm, args.egt, args.csv, folderGTC, allGTCFolder, args.genomeReference, vcfFolder,
            args.threads, args.outputName, args.mem, logFile)

    VCFNoLift = removeLiftOverProblems(VCF, args.listError, args.bcftools, args.outputName, vcfFolder)

    if args.analysis:
        execute(f"{args.bcftools} filter -e 'INFO/GenTrain_Score<=0.7' {VCFNoLift} | "
                f"{args.bcftools} annotate -x 'INFO' | "
                f"{args.bcftools} annotate -x 'FORMAT' -Oz -o {vcfFolder}/{args.outputName}_GenTrain_noAnnotation.vcf.gz",
                logFile)
        execute(f"{args.bcftools} index {vcfFolder}/{args.outputName}_GenTrain_noAnnotation.vcf.gz", logFile)

        sheetFile = open(sampleSheet)
        headerFile = open(f"{vcfFolder}/toChangeID.txt", "w")
        header = True
        for line in sampleSheet:
            if header:
                if "Sample_ID,SentrixBarcode_A,SentrixPosition_A,Path" in line:
                    header = False
            else:
                ID, barcode, position, path = line.strip().split()
                IDWithoutBatch = ID.split("_B")[0]
                oldID = f"{barcode}_{position}"
                headerFile.write(f"{oldID} {IDWithoutBatch}\n")
        headerFile.close()
        execute(f"{args.bcftools} reheader -s {vcfFolder}/toChangeID.txt -o {vcfFolder}/{args.outputName}_toAnalysis.vcf.gz "
                  f"{vcfFolder}/{args.outputName}_GenTrain_noAnnotation.vcf.gz", logFile)
        execute(f"{args.bcftools} index {vcfFolder}/{args.outputName}_toAnalysis.vcf.gz", logFile)
    logFile.close()