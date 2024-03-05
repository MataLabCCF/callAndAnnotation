import os
import gzip
import argparse

def removeLiftOverProblemsNew(vcfFileName, errors, bcftools, outputName, outputFolder):
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
    os.system(f"bcftools index {outputFolder}/{outputName}_NoLiftProblems.vcf.gz")

    return f"{outputFolder}/{outputName}_NoLiftProblems.vcf.gz"

def removeLiftOverProblems(vcfFileName, errors, bcftools, outputName, outputFolder):
    dictID = {}
    vcfFile = gzip.open(vcfFileName)

    var = 0

    header = True
    for line in vcfFile:
        line = line.decode("utf-8")
        if header:
            if "#CHROM" in line:
                header = False
        else:
            var = var+1
            lineToSplit = line[0:150]
            split = lineToSplit.split()
            dictID[split[2]] = f"{split[0]}\t{split[1]}\n"
    vcfFile.close()

    print(f"The VCF has {var} variants")

    fileOut = open(f"{outputFolder}/errorToRemoveLift.txt", "w")
    fileError = open(errors)
    for line in fileError:
        split = line.split()
        if split[0] in dictID:
            fileOut.write(dictID[split[0]])
    fileOut.close()

    os.system(f"{bcftools} view -T ^{outputFolder}/errorToRemoveLift.txt -Oz -o {outputFolder}/{outputName}_NoLiftProblems.vcf.gz {vcfFileName}")

    return f"{outputFolder}/{outputName}_NoLiftProblems.vcf.gz"

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Prepare Annotation')

    requiredGeneral = parser.add_argument_group("Required arguments for all steps")
    requiredGeneral.add_argument('-v', '--vcf', help='VCF input file', required=True)
    requiredGeneral.add_argument('-O', '--outputFolder', help='Name of output path', required=True)
    requiredGeneral.add_argument('-o', '--outputName', help='Name of output name', required=True)
    requiredGeneral.add_argument('-e', '--errors', help='File with LiftOver problems and unsolved tri-allelics', required=True)
    #requiredGeneral.add_argument('-r', '--reference', help='VCF reference to be used in eagle', required=True)


    programs = parser.add_argument_group("Programs")
    programs.add_argument('-B', '--bcftools', help='bcftools path', required=False, default="bcftools")
    #programs.add_argument('-E', '--eagle', help='eagle path ', required=True, default = "eagle")


    args = parser.parse_args()
    file = removeLiftOverProblemsNew(args.vcf, args.errors, args.bcftools, args.outputName, args.outputFolder)

    print(f"File to annotate is in {file}")