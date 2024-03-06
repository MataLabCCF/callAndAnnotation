import os
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Gentrain and reheader')

    requiredGeneral = parser.add_argument_group("Required arguments for all steps")
    requiredGeneral.add_argument('-o', '--outputFolder', help='Name of output folder', required=True)
    requiredGeneral.add_argument('-O', '--outputName', help='Name of output name', required=True)
    requiredGeneral.add_argument('-s', '--sampleSheet', help='Sample sheet used to perform variant calling', required=True)

    optional = parser.add_argument_group("Optional arguments")
    optional.add_argument('-r', '--removeBatch', help='Flag used to remove batch ID', default=False, action="store_true")
    args = parser.parse_args()


    removeBatch = args.removeBatch
    outFolder = args.outputFolder
    outName = args.outputName
    sampleSheetName = args.sampleSheet


    os.system(f'bcftools filter --threads 45 -e "INFO/GenTrain_Score<=0.7 " {outFolder}/{outName}_NoLiftProblems.vcf.gz | '
               f'bcftools annotate --threads 45 -x "INFO" | bcftools annotate --threads 45 -x "FORMAT" -Oz -o '
               f'{outFolder}/{outName}_GenTrain.vcf.gz')



    sampleSheet = open(sampleSheetName)
    fileToChangeID = open(f"{outFolder}/{outName}_convertIDs.txt", "w")

    header = True
    for line in sampleSheet:
        if header:
            if  "Sample_ID,SentrixBarcode_A,SentrixPosition_A,Path" in line:
                header = False
        else:
            ID, barcode, position, path = line.strip().split(",")
            if removeBatch:
                ID = ID.split("_B")[0]

            VCFCode = f"{barcode}_{position}"

            fileToChangeID.write(f"{VCFCode} {ID}\n")

    fileToChangeID.close()

    os.system(f"bcftools reheader -s {outFolder}/{outName}_convertIDs.txt -o {outFolder}/{outName}_Header.vcf.gz {outFolder}/{outName}_GenTrain.vcf.gz")
