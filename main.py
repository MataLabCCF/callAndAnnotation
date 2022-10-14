import myvariant
import argparse
import gzip
import time
import os
import numpy as np


mv = myvariant.MyVariantInfo()

#=========================================================== Utils ====================================================
def createFolder(folderName, logFile):
    line = f"mkdir {folderName}"
    execute(line, logFile)

    return folderName

def execute(commandLine, logFile):
    start = time.time()
    print(commandLine)
    os.system(commandLine)
    end = time.time()
    diff = end-start

    logFile.write(f'{commandLine}\nTime: {diff} seconds\n\n')


#=========================================================== QC ====================================================
def gtc2VCF(bcftools, bpm, egt, csv, folderGTC, allGTCFolder, genomeReference, vcfFolder, threads, outputName, logFile):
    for folder in os.listdir(folderGTC):
        line = f"cp {folderGTC}/{folder}/* {allGTCFolder}"
        execute(line, logFile)

    line = f"{bcftools} +gtc2vcf --no-version -Oz --bpm {bpm} --egt {egt} --csv {csv} --gtcs {allGTCFolder} " \
           f"--fasta-ref {genomeReference} --output {vcfFolder}/{outputName}.vcf.gz --threads {threads} " \
           f"--use-gtc-sample-names --adjust-clusters"
    execute(line, logFile)

    line = f"{bcftools} sort -Oz -o {vcfFolder}/{outputName}_Sort.vcf.gz {vcfFolder}/{outputName}.vcf.gz"
    execute(line, logFile)

    line = f"{bcftools} norm -Oz -c x -f {genomeReference} --threads {threads}" \
           f" -o {vcfFolder}/{outputName}_Norm.vcf.gz {vcfFolder}/{outputName}_Sort.vcf.gz"
    execute(line, logFile)

    line = f"{bcftools} index {vcfFolder}/{outputName}_Norm.vcf.gz --threads {threads}"
    execute(line, logFile)

    return f"{vcfFolder}/{outputName}_Norm.vcf.gz"

def generateGTC(iaap, bpm, egt, folder, outFolder, threads, logFile):
    threadsTest = 4

    line = f"{iaap} gencall {bpm} {egt} -f {folder} --output-gtc {outFolder} --gender-estimate-call-rate-threshold " \
           f"-0.1 -t {threadsTest}"
    execute(line, logFile)

#def generatePED(iaap, bpm, egt, folder, outFolder, threads):
#    threadsTest = 4

#    line = f"{iaap} gencall {bpm} {egt} -f {folder} --output-gtc {outFolder} --gender-estimate-call-rate-threshold -0.1 -t {threadsTest}"
#    execute(line)

def readFileQC(qcFile):
    dictQC = {}

    if qcFile == "":
        #INCLUDE, not exclude
        #dictQC["F_MISSING"] = " < 0.05"  # VCF
        dictQC["GenTrain_Score"] = " >= 0.67"  # VCF
        dictQC["Cluster_Sep"] = " > 0.425"  # VCF
        dictQC["devTHETA_AA"] = " <= 0.05"  # VCF
        dictQC["devTHETA_BB"] = " <= 0.05"  # VCF
        dictQC["meanR_AA"] = " >= 0.2"  # VCF
        dictQC["meanR_BB"] = " >= 0.2"  # VCF
        dictQC["meanR_AB"] = " >= 0.2"  # VCF
    else:
        pass
        #TODO

    return dictQC


def basicVCFQC(VCF, parametersQC, bcftools, vcfFolder, vcfName, logFile):
    line = f'{bcftools}  view -i \"'

    first = True

    for parameter in parametersQC:
        if first:
            line = f'{line}INFO/{parameter}{parametersQC[parameter]}'
            first = False
        else:
            line = f'{line} && INFO/{parameter}{parametersQC[parameter]}'

    line = f'{line}\" -Oz -o {vcfFolder}/{vcfName}_QC.vcf.gz {VCF}'
    execute(line, logFile)

    return f"{vcfFolder}/{vcfName}_QC.vcf.gz"


def basicGenotypingQC(VCF, folder, name, plink2, logFile):
    folderPGEN = f"{folder}/PGEN"
    execute(f"mkdir {folderPGEN}", logFile)

    #Convert
    execute(f"{plink2} --vcf {VCF} --make-pfile --out {folderPGEN}/{name}", logFile)

    #Remove missing
    execute(f"{plink2} --pfile {folderPGEN}/{name} --geno 0.05 --mind 0.05 --make-pgen --out {folderPGEN}/{name}_missing", logFile)

    #Remove big deviation on HWE
    execute(f"{plink2} --pfile {folderPGEN}/{name}_missing --hwe 1e-10 --make-pgen --out {folderPGEN}/{name}_HWE", logFile)

    #PGEN to VCF
    execute(f"{plink2} --pfile {folderPGEN}/{name}_HWE --recode vcf-iid --out {folderPGEN}/{name}_QC", logFile)
    return f"{folderPGEN}/{name}_QC.vcf"

#=========================================================== Anotation ====================================================


def readPreviousSearchFile(MyVariantFile):
    file = open(MyVariantFile)
    dictInfo = {}

    header = True
    for line in file:
        if header:
            headerLine = line.strip().split("\t")
            for i in range(0, len(headerLine)):
                if headerLine[i] == "Query":
                    keyIndex = i
            header = False
        else:
            split = line.strip().split("\t")
            ID = split[keyIndex]
            dictInfo[ID] = {}

            for i in range(0, len(split)):
                if i != keyIndex:
                    dictInfo[ID][headerLine[i]] = ""
                    dictInfo[ID][headerLine[i]] = split[i]

    return dictInfo


def getDataToPrint(data):
    if isinstance(data, list):
        if len(data) < 1:
            toReturn = ""
        else:
            toReturn = ('; ').join(data)
    else:
        toReturn = data

    return toReturn

def readVCFs(vcfName, variantInfo):

    print(f"Open VCF file ({vcfName})")
    toReturn = []
    searches = {}

    if vcfName[-2:] == ".gz":
        file = gzip.open(vcfName)
    else:
        file = open(vcfName)

    header = True
    countToSearch = 0
    count = 0
    for line in file:
        if vcfName[-2:] == ".gz":
            line = line.decode("utf-8")

        if header:
            if "#CHROM" in line:
                header = False
        else:
            split = line.strip().split()
            prefix = f"chr{split[0]}:g.{split[1]}"
            key = f"{prefix}{split[3]}>{split[4]}"
            if key not in variantInfo:
                if key not in searches:
                    searches[key] = {}
                    searches[key]["Allele 1"] = split[3]
                    searches[key]["Allele 2"] = split[4]
                    searches[key]["prefix"] = prefix
                    searches[key]["alleles"] = split[-1]
                    countToSearch = countToSearch + 1
            count = count + 1
    print(f"\tThe VCF has {count} variants, which {countToSearch} is new \n")
    return searches


def dbsnpData(dictData):
        alleles = ""
        allelesFreq = ""

        if 'rsid' in dictData:
            rsID = dictData['rsid']
        else:
            rsID = "NA"

        if 'ref' in dictData:
            ref = dictData['ref']
        else:
            ref = "NA"

        if 'alleles' in dictData:
            for observation in dictData['alleles']:
                allele = observation['allele']
                if 'topmed' not in observation['freq']:
                    freq = "NA"
                else:
                    freq = observation['freq']['topmed']

                if alleles == "":
                    alleles = f'{allele}'
                else:
                    alleles = f'{alleles},{allele}'

                if allelesFreq == "":
                    allelesFreq = f'{freq}'
                else:
                    allelesFreq = f'{allelesFreq},{freq}'
        else:
            alleles = "NA"
            allelesFreq = "NA"

        return alleles, allelesFreq, rsID, ref


def caddData(dictData):
    polyphen = ""
    sift = ""

    geneID = ""
    geneName = ""

    annotype = dictData['annotype']
    consdetail = dictData['consdetail']
    consequence = dictData['consequence']
    consscore = str(dictData['consscore'])
    phredLike = dictData['phred']

    if 'gene' in dictData:

        if isinstance(dictData['gene'], list):
            geneID = dictData['gene'][0]['gene_id']
            geneName = dictData['gene'][0]['genename']

            for i in range(1, len(dictData['gene'])):
                geneID = geneID + "; " + dictData['gene'][i]['gene_id']
                geneName = geneName + "; " + dictData['gene'][i]['genename']
        else:
            geneID = dictData['gene']['gene_id']
            geneName = dictData['gene']['genename']
    else:
        geneID = ""
        geneName = ""

    if 'polyphen' in dictData:
        if isinstance(dictData['polyphen'], list):
            sift = dictData['polyphen'][0]['cat']
            for i in range(len(dictData['polyphen'])):
                sift = sift + "; " + dictData['polyphen'][i]['cat']
        else:
            polyphen = dictData['polyphen']['cat']

    if 'sift' in dictData:
        if isinstance(dictData['sift'], list):
            sift = dictData['sift'][0]['cat']
            for i in range(len(dictData['sift'])):
                sift = sift + "; " + dictData['sift'][i]['cat']
        else:
            sift = dictData['sift']['cat']

    return annotype, consdetail, consequence, consscore, geneID, geneName, polyphen, sift, phredLike


def getConditionName(data):
    name = ""
    if isinstance(data, list):
        for i in range(0, len(data)):
            # print(data[i])
            if 'name' in data[i]:
                if name == "":
                    name = data[i]['name']
                else:
                    name = name + "; " + data[i]['name']
    else:
        if 'name' in data:
            name = data['name']

    # print(name)
    return name


def getConditionIdentifiers(data):
    mondo = []
    medgen = []
    omim = []
    orphanet = []
    HPO = []

    # print(data)

    if isinstance(data, list):
        # print("Is list")
        for i in range(len(data)):
            if 'identifiers' in data[i]:
                if 'mondo' in data[i]['identifiers']:
                    mondo.append(data[i]['identifiers']['mondo'])
                if 'medgen' in data[i]['identifiers']:
                    medgen.append(data[i]['identifiers']['medgen'])
                if 'omim' in data[i]['identifiers']:
                    omim.append(data[i]['identifiers']['omim'])
                if 'orphanet' in data[i]['identifiers']:
                    orphanet.append(data[i]['identifiers']['orphanet'])
                if 'human_phenotype_ontology' in data[i]['identifiers']:
                    HPO.append(data[i]['identifiers']['human_phenotype_ontology'])
    else:
        # print("Is not list")
        if 'identifiers' in data:
            if 'mondo' in data['identifiers']:
                mondo.append(data['identifiers']['mondo'])
            if 'medgen' in data['identifiers']:
                medgen.append(data['identifiers']['medgen'])
            if 'omim' in data['identifiers']:
                omim.append(data['identifiers']['omim'])
            if 'orphanet' in data['identifiers']:
                orphanet.append(data['identifiers']['orphanet'])
            if 'human_phenotype_ontology' in data['identifiers']:
                HPO.append(data['identifiers']['human_phenotype_ontology'])

    return medgen, mondo, omim, orphanet, HPO


def getRCVData(data):
    omimList = []
    mondoList = []
    medgenList = []
    orphanetList = []
    HPOList = []

    preferredName = ""
    name = ""
    # NAME
    if isinstance(data, list):
        for i in range(0, len(data)):
            returnedName = getConditionName(data[i]['conditions'])
            medgen, mondo, omim, orphanet, HPO = getConditionIdentifiers(data[i]['conditions'])

            omimList = omimList + omim
            mondoList = mondoList + mondo
            medgenList = medgenList + medgen
            orphanetList = orphanetList + orphanet
            HPOList = HPOList + HPO

            if name == "":
                name = returnedName
            else:
                name = name + "; " + returnedName

            if 'preferred_name' in data[i]:
                if preferredName == "":
                    preferredName = data[i]['preferred_name']
                else:
                    preferredName = preferredName + "; " + data[i]['preferred_name']

    else:
        name = getConditionName(data['conditions'])
        medgen, mondo, omim, orphanet, HPO = getConditionIdentifiers(data['conditions'])

        omimList = omimList + omim
        mondoList = mondoList + mondo
        medgenList = medgenList + medgen
        orphanetList = orphanetList + orphanet
        HPOList = HPOList + HPO

        if 'preferred_name' in data:
            preferredName = data['preferred_name']

    return name, preferredName, omimList, mondoList, medgenList, orphanetList, HPOList


def makeSearch(toSearch, variantInfo, fields):
    print("\tQuery")
    notFound = []
    for data in mv.getvariants(toSearch, fields=fields, returnall=True):
        if "notfound" in data:
            notFound.append(data['query'])
        else:
            if isinstance(data, dict):
                query = data["query"]

                if 'cadd' in data:
                    annotype, consdetail, consequence, consscore, geneID, geneName, polyphen, sift, phred = caddData(
                        data["cadd"])
                else:
                    annotype = ""
                    consdetail = ""
                    consequence = ""
                    consscore = ""
                    geneID = ""
                    geneName = ""
                    polyphen = ""
                    sift = ""
                    phred = ""

                if 'dbsnp' in data:
                    alleles, freqs, rs, ref = dbsnpData(data["dbsnp"])
                else:
                    rs = "NA"
                    alleles = "NA"
                    freqs = "NA"
                    ref = "NA"

                if 'clinvar' in data:
                    name, preferredName, omimList, mondoList, medgenList, orphanetList, humanPhenotypeOntologyList = getRCVData(
                        data["clinvar"]["rcv"])
                else:
                    name = ""
                    preferredName = ""
                    omimList = []
                    mondoList = []
                    medgenList = []
                    orphanetList = []
                    humanPhenotypeOntologyList = []

                omimOut = getDataToPrint(omimList)
                mondoOut = getDataToPrint(mondoList)
                medgenOut = getDataToPrint(medgenList)
                orphanetOut = getDataToPrint(orphanetList)
                humanPhenotypeOntologyOut = getDataToPrint(humanPhenotypeOntologyList)

                annotypeOut = getDataToPrint(annotype)
                consdetailOut = getDataToPrint(consdetail)
                consequenceOut = getDataToPrint(consequence)
                consscoreOut = getDataToPrint(consscore)

                keys = ['rs', 'Alleles', 'Ref', 'Freqs', 'CADD Phred', 'CADD Annotype',
                        'CADD Consequence Detail', 'CADD Consequence', 'CADD Consequence Score', 'CADD Gene ID',
                        'CADD Gene Name', 'Sift Cathegory', 'Polyphen Cathegory', 'CLINVAR Name',
                        'CLINVAR Preferred Name',
                        'CLINVAR Preferred Name', 'OMIM IDs', 'MONDO IDs', 'MEDGEN IDs', 'Orphanet IDs',
                        'Human Phenotype Ontology']

                variantInfo[query] = {}
                for key in keys:
                    variantInfo[query][key] = ""

                variantInfo[query]['rs'] = rs
                variantInfo[query]['Freqs'] = freqs
                variantInfo[query]['Alleles'] = alleles
                variantInfo[query]['Ref'] = ref


                variantInfo[query]['CADD Phred'] = phred
                variantInfo[query]['CADD Annotype'] = annotypeOut
                variantInfo[query]['CADD Consequence Detail'] = consdetailOut
                variantInfo[query]['CADD Consequence'] = consequenceOut
                variantInfo[query]['CADD Consequence Score'] = consscoreOut
                variantInfo[query]['CADD Gene ID'] = geneID
                variantInfo[query]['CADD Gene Name'] = geneName
                variantInfo[query]['Sift Cathegory'] = sift
                variantInfo[query]['Polyphen Cathegory'] = polyphen
                variantInfo[query]['CLINVAR Name'] = name
                variantInfo[query]['CLINVAR Preferred Name'] = preferredName
                variantInfo[query]['OMIM IDs'] = omimOut
                variantInfo[query]['MONDO IDs'] = mondoOut
                variantInfo[query]['MEDGEN IDs'] = medgenOut
                variantInfo[query]['Orphanet IDs'] = orphanetOut
                variantInfo[query]['Human Phenotype Ontology'] = humanPhenotypeOntologyOut

    return notFound, variantInfo

def searchSNPs(variantInfo, toSearchDict):
    print("Get annotation")

    fields = ['clinvar.rcv.conditions.identifiers', 'clinvar.rcv.conditions.name', 'clinvar.rcv.preferred_name',
              'dbsnp.alleles', 'dbsnp.rsid', 'dbsnp.ref', 'cadd.gene.gene_id', 'cadd.gene.genename', 'cadd.annotype',
              'cadd.consdetail', 'cadd.consequence', 'cadd.consscore', 'cadd.sift.cat', 'cadd.polyphen.cat',
              'cadd.phred']

    toSearch = list(toSearchDict.keys())
    notFound, variantInfo = makeSearch(toSearch, variantInfo, fields)

    notFoundToSearch = {}
    for ID in notFound:
        A2 = toSearchDict[ID]["Allele 1"]
        A1 = toSearchDict[ID]["Allele 2"]
        prefix = toSearchDict[ID]["prefix"]

        notFoundToSearch[f"{prefix}{A1}>{A2}"] = ID

    toSearchNotFound = list(notFoundToSearch.keys())
    notFound, variantInfo = makeSearch(toSearchNotFound, variantInfo, fields)

    return variantInfo

def isNotNull(data):
    if data != "NA" and data != "":
        return True
    return False

def getAFs(data):
    toReturn = []
    split = data.strip().split(",")

    for AF in split:
        if AF == "NA" or AF == "":
            return split, False
        else:
            toReturn.append(float(AF))
    return toReturn, True

def saveSearch(fileOutName):
    header = True
    fileReference = open(fileOutName, "w")
    for query in variantInfo:
        if header:
            fileReference.write("Query")
            for field in variantInfo[query]:
                fileReference.write(f"\t{field}")
            fileReference.write("\n")
            header = False
        fileReference.write(f"{query}")
        for field in variantInfo[query]:
            fileReference.write(f"\t{variantInfo[query][field]}")
        fileReference.write("\n")

def openCorrespondenceFile(fileInput):
    fileCorrespondence = open(fileInput)
    dictCorrespondence = {}
    for line in fileCorrespondence:
        split = line.strip().split("\t")
        dictCorrespondence[split[0]] = split[1]
    return dictCorrespondence

def createAnnotationTables(line, correspondence, dictCorrespondence, folderAnot):
    headerIDsTemp = line.strip().split()
    headerIDs = []

    dictOut = {}

    for i in range(len(headerIDsTemp)):
        IDTemp = headerIDsTemp[i]
        if i < 9:
            headerIDs.append(IDTemp)
        else:
            #Change to correspondence list if exists
            if correspondence != "":
                if IDTemp in dictCorrespondence:
                    ID = dictCorrespondence[IDTemp]
                else:
                    ID = f"NoCorrespondence_{IDTemp}"
            else:
                ID = IDTemp

            #Create a unique file to each ID, even if the name is the same
            if ID not in dictOut:
                dictOut[ID] = open(f"{folderAnot}/Annot_{ID}.tsv", 'w')
                headerIDs.append(ID)
            else:
                count = 1
                while f"{ID}_{count}" in headerIDs:
                    count = count + 1

                headerIDs.append(f"{ID}_{count}")
                dictOut[f"{ID}_{count}"] = open(f"{folderAnot}/Annot_{ID}_{count}.tsv", 'w')

    return dictOut, headerIDs

def openVCFAndSearch(vcfName, variantInfo, folderAnot, name, correspondence):
    toSearch = readVCFs(vcfName, variantInfo)
    variantInfo = searchSNPs(variantInfo, toSearch)
    saveSearch(f"{folderAnot}/{name}_AnnotTable.tsv")

    if correspondence != "":
        dictCorrespondence = openCorrespondenceFile(correspondence)

    print(f"Open VCF file again({vcfName}) ")

    if vcfName[-2:] == ".gz":
        file = gzip.open(vcfName)
    else:
        file = open(vcfName)

    headerOnOutput = False
    header = True
    for line in file:
        if vcfName[-2:] == ".gz":
            line = line.decode("utf-8")
        if header:
            if "#CHROM" in line:
                dictOut, headerIDs = createAnnotationTables(line, correspondence, dictCorrespondence, folderAnot)
                header = False
        else:
            split = line.strip().split()
            prefix = f"chr{split[0]}:g.{split[1]}"
            key = f"{prefix}{split[3]}>{split[4]}"

            if key in variantInfo:
                query = key
                A0 = split[3]
                A1 = split[4]
            else:
                key = f"{prefix}{split[4]}>{split[3]}"
                if key in variantInfo:
                    query = key
                    A0 = split[4]
                    A1 = split[3]
                else:
                    query = "notFound"

            if query != "notFound":
                if not headerOnOutput:
                    headerOnOutput = True
                    for ind in dictOut:
                        dictOut[ind].write("Query\tGenotype")
                        for item in variantInfo[query]:
                            dictOut[ind].write(f"\t{item}")
                        dictOut[ind].write("\n")


                #Select CADD Phred if it is valid
                validPhred = False
                if isNotNull(variantInfo[query]['CADD Phred']):
                    phred = float(variantInfo[query]['CADD Phred'])
                    validPhred = True

                # Select MAF if it is valid (2nd higher AF)
                AFList, validAF = getAFs(variantInfo[query]['Freqs'])

                if validAF:
                    if len(AFList) == 1:
                        validAF = False
                    if len(AFList) == 2:
                        maf = min(AFList)
                    if len(AFList) > 2:
                        sortedAF = np.sort(AFList)
                        maf = sortedAF[-2]

                if validPhred and validAF:
                    if variantInfo[query]['CADD Phred'] > 15 and maf < 0.05:
                        for i in range(9, len(split)):
                            ID = headerIDs[i]
                            information = split[i].split(":")
                            if "1" in information[0]:
                                information[0].replace("0", A0).replace("1", A1)
                                dictOut[ID].write(f"{query}\t{information[0]}")
                                for item in variantInfo[query]:
                                    dictOut[ID].write(f"\t{variantInfo[query][item]}")
                                dictOut[ID].write(f"\n")
    for ID in dictOut:
        dictOut[ID].close()






if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='From Illumina to annotation')

    requiredGeneral = parser.add_argument_group("Required arguments for all steps")
    requiredGeneral.add_argument('-o', '--outputFolder', help='Name of output folder', required=True)
    requiredGeneral.add_argument('-O', '--outputName', help='Name of output name', required=True)

    requiredCalling = parser.add_argument_group("Required arguments for variant calling")
    requiredCalling.add_argument('-b', '--bpm', help='BPM File from Illumina', required=True)
    requiredCalling.add_argument('-c', '--csv', help='CSV manifest file', required=True)
    requiredCalling.add_argument('-e', '--egt', help='EGT File from Illumina', required=True)
    requiredCalling.add_argument('-f', '--folder', help='Folder with files from illumina', required=True)
    requiredCalling.add_argument('-r', '--genomeReference', help='Human genome reference', required=True)
    requiredCalling.add_argument('-l', '--locusSummary', help='Human genome reference', required=True)

    programs = parser.add_argument_group("Programs")
    programs.add_argument('-B', '--bcftools', help='bcftools', required=False, default="bcftools")
    programs.add_argument('-i', '--iaap', help='Illumina Array Analysis Platform path ', required=True)
    programs.add_argument('-p', '--plink2', help='Plink2 path ', required=True)

    optional = parser.add_argument_group("Optional arguments")
    optional.add_argument('-t', '--threads', help='Number of threads', default=96, type=int, required=False)
    optional.add_argument('-q', '--qc', help='QC parameters to VCF', required=False, default="")
    optional.add_argument('-P', '--previousSearch', help='File with MyVariant previus search', required=False,
                          default="")
    optional.add_argument('-C', '--correspondenceList', help='List with the correspondence Illumina to patientes ID',
                          required=False, default="")
    args = parser.parse_args()

    logFile = open(f"{args.outputFolder}/{args.outputName}.log", 'w')
    folder = createFolder(args.outputFolder, logFile)



    folderGTC = createFolder(folder+"/GTC/", logFile)
    generateGTC(args.iaap, args.bpm, args.egt, args.folder, folderGTC, args.threads, logFile)
    #generatePED(args.iaap, args.bpm, args.egt, args.folder, folderGTC, args.threads)

    allGTCFolder = createFolder(folder+"/GTCAll/", logFile)
    vcfFolder = createFolder(folder+"/VCFs/", logFile)
    VCF = gtc2VCF(args.bcftools, args.bpm, args.egt, args.csv, folderGTC, allGTCFolder, args.genomeReference, vcfFolder, args.threads, args.outputName, logFile)
    parametersQC = readFileQC(args.qc)

    VCFQC = basicVCFQC(VCF, parametersQC, args.bcftools, vcfFolder, args.outputName, logFile)

    VCFFinal = basicGenotypingQC(VCFQC, folder, args.outputName, args.plink2, logFile)

    anotFolder = createFolder(folder + "/Anot/", logFile)

    start = time.time()

    variantInfo = {}
    if args.previousSearch != "":
        variantInfo = readPreviousSearchFile(args.previousSearch)

    variantInfo = openVCFAndSearch(VCFFinal, variantInfo, anotFolder, args.outputName, args.correspondenceList)

    end = time.time()
    diff = end - start

    logFile.write(f'Annotation using MyVariants\nTime: {diff} seconds\n\n')
    logFile.close()