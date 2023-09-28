import myvariant
import argparse
import gzip
import time
import sys
import os
import re
import csv
import numpy as np
from xlsxwriter.workbook import Workbook



mv = myvariant.MyVariantInfo()

#=========================================================== Utils ====================================================
def createFolder(folderName, logFile):
    line = f"mkdir {folderName}"
    execute(line, logFile)

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
def gtc2VCF(bcftools, bpm, egt, csv, folderGTC, allGTCFolder, genomeReference, vcfFolder, threads, outputName, logFile):
    for file in os.listdir(folderGTC):
        line = f"cp {folderGTC}/{file} {allGTCFolder}"
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

def generateGTC(iaap, bpm, egt, folder, outFolder, outputName, batchList, threads, logFile):
    threadsTest = 4

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
            IDIllumina, number = illuminaID.split("-")
            pathToLook = path + "/" + sentrixBarcode
            externalID = dictCodes[illuminaID]["ExternalID"] + "_" + dictCodes[illuminaID]["Batch"]

            if os.path.exists(pathToLook):
                outputFile.write(f"{externalID},{sentrixBarcode},{sentrixPosition},{pathToLook}\n")
            else:
                input(f"Not found: {pathToLook}")

        fileSentrix.close()

    outputFile.close()

    line = f"{iaap} gencall {bpm} {egt} -s {folder}/{outputName}_SampleSheet.csv --output-gtc {outFolder} --gender-estimate-call-rate-threshold " \
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

def getInfo(line):
    infoDict = {}
    split = line.strip().split()[7].split(";")
    for fieldVCF in split:
        if "=" in fieldVCF:
            field,value = fieldVCF.split("=")
            infoDict[field]=value
    return infoDict

def wasAccepted(value, parameter):
    split = parameter.split(" ")
    cutoff = float(split[-1])

    if ">=" in parameter:
        if value >= cutoff:
            return "ok"
        return "fail"
    elif "<=" in parameter:
        if value <= cutoff:
            return "ok"
        return "fail"
    elif ">" in parameter:
        if value > cutoff:
            return "ok"
        return "fail"
    elif "<" in parameter:
        if value < cutoff:
            return "ok"
        return "fail"
    else:
        sys.exit(f"Unknown parameter: {parameter}")


def basicVCFQCToIncludeAsInfo(VCF, parametersQC, outFolder):

    dictInfo= {}
    if ".gz" in VCF:
        file = gzip.open(VCF)
    else:
        file = open(VCF)

    header = True
    for line in file:
        if ".gz" in VCF:
            line = line.decode("utf-8")

        if header:
            if "#CHROM" in line:
                header = False
        else:
            infoDict = getInfo(line)
            ID = line.strip().split()[2]

            failQC = ""

            for field in parametersQC:
                if field in infoDict:
                    if wasAccepted(float(infoDict[field]),parametersQC[field]) == "fail":
                        if failQC == "":
                            failQC = field
                        else:
                            failQC = f"{failQC}, {field}"

            if failQC == "":
                dictInfo[ID] = "PASS"
            else:
                dictInfo[ID] = failQC

    return dictInfo

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
                geneID = f"{geneID}; {dictData['gene'][i]['gene_id']}"
                geneName = f"{geneName}; {dictData['gene'][i]['genename']}"
        else:
            geneID = dictData['gene']['gene_id']
            geneName = dictData['gene']['genename']
    else:
        geneID = ""
        geneName = ""

    if 'polyphen' in dictData:
        if isinstance(dictData['polyphen'], list):
            polyphen = dictData['polyphen'][0]['cat']
            for i in range(len(dictData['polyphen'])):
                polyphen = f"{polyphen}; {dictData['polyphen'][i]['cat']}"
        else:
            polyphen = dictData['polyphen']['cat']

    if 'sift' in dictData:
        if isinstance(dictData['sift'], list):
            sift = dictData['sift'][0]['cat']
            for i in range(len(dictData['sift'])):
                sift = f"{sift}; {dictData['sift'][i]['cat']}"
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
                    name = f"{name}; {data[i]['name']}"
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
                name = f"{name}; {returnedName}"

            if 'preferred_name' in data[i]:
                if preferredName == "":
                    preferredName = data[i]['preferred_name']
                else:
                    preferredName = f"{preferredName}; {data[i]['preferred_name']}"

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

def fixAlleles(original, new):
    dataOriginal = re.split(r'(\d+)', original)
    dataNew = re.split(r'(\d+)', new)

    newId = ""
    for i in range(0, len(dataNew)-1):
        newId = newId+dataNew[i]

    newId = newId+dataOriginal[-1]
    return newId



def getIDsNonHG38(listVar):
    print(f"\tConverting HG38 to HG19")

    listVarHg19 = []
    print(f"We will try to convert {len(listVar)} variants")

    for var in listVar:
        # preparation to future
        dbSNP = ""
        info = re.split(r'(\d+)', var)
        alleles = info[-1].split(">")
        flag = 0



        search = mv.getvariant(var, fields="dbsnp.rsid", assembly="hg38", returnall=True)
        #print(search)
        if search and ('dbsnp' in search):
            dbSNP = search['dbsnp']['rsid']
        else:
            newVar = ""
            for i in range(len(info) - 1):
                newVar = newVar + info[i]
            newVar = newVar + f"{alleles[1]}>{alleles[0]}"
            search = mv.getvariant(newVar, fields="dbsnp.rsid", assembly="hg38", returnall=True)
            if search:
                dbSNP = search['dbsnp']['rsid']
            else:
                print(f"\tWe could not find {var} or {newVar}")

        if dbSNP:
            search = mv.getvariant(dbSNP, returnall=True)
            if isinstance(search, list):
                for answer in search:
                    allelesAnswer = re.split(r'(\d+)', answer['_id'])[-1].split(">")
                    if (allelesAnswer[0] == alleles[0] and allelesAnswer[1] == alleles[1]) or (
                            allelesAnswer[0] == alleles[1] and allelesAnswer[1] == alleles[0]):
                        listVarHg19.append(answer['_id'])
                        flag = 1
            else:
                allelesAnswer = re.split(r'(\d+)', search['_id'])[-1].split(">")
                if (allelesAnswer[0] == alleles[0] and allelesAnswer[1] == alleles[1]) or (
                        allelesAnswer[0] == alleles[1] and allelesAnswer[1] == alleles[0]):
                    listVarHg19.append(search['_id'])
                    flag = 1
        if flag != 1:
            print(f"\t\tUnable to convert {var}")

    print(f"We got {len(listVarHg19)} variants")
    return listVarHg19

def makeSearch(toSearch, variantInfo, fields):
    print("Query")
    isHg38 = False

    notFound = []

    if isHg38:
        #I will not use because some variants are with hg19 after calling with hg38 bpm and norm with hg38
        notFound, toSearch19, dict19To38 = getIDsNonHG38(toSearch)
    else:
        toSearch19 = toSearch

    print(f"\tSearching {len(toSearch19)}")
    for var in toSearch19:

        #print(f"Looking for {var}")

        found = True
        data = mv.getvariant(var, fields=fields, returnall=True)

        if not data:
            info = re.split(r'(\d+)', var)
            alleles = info[-1].split(">")
            newVar = ""
            for i in range(len(info) - 1):
                newVar = newVar + info[i]
            newVar = newVar + f"{alleles[1]}>{alleles[0]}"

            data = mv.getvariant(newVar, fields=fields, returnall=True)
            if not data:
                print(f"Unable to find {var} or {newVar}")
                found = False
        if found:
            #print(type(data))
            if isinstance(data, dict):
                query = data["_id"]

                if 'cadd' in data:
                    #print(f"\tWe have cadd for {query} -> {data['cadd']}")
                    annotype, consdetail, consequence, consscore, geneID, geneName, polyphen, sift, phred = caddData(
                        data["cadd"])
                else:
                    print(f"\tWe do not have cadd for {query}")
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
                    if 'hg19' in data["clinvar"]:
                        posHG19 = data["clinvar"]['hg19']["start"]
                    if 'hg38' in data["clinvar"]:
                        posHG38 = data["clinvar"]['hg38']["start"]
                    if 'chrom' in data['clinvar']:
                        chrom = data['clinvar']['chrom']

                else:
                    name = ""
                    preferredName = ""
                    posHG19 = ""
                    posHG38 = ""
                    chrom = ""
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

                keys = ['rs', 'Alleles', 'Ref', 'Chrom' ,'Pos hg19','Pos hg38', 'Freqs', 'CADD Phred', 'CADD Annotype',
                        'CADD Consequence Detail', 'CADD Consequence', 'CADD Consequence Score', 'CADD Gene ID',
                        'CADD Gene Name', 'Sift Cathegory', 'Polyphen Cathegory', 'CLINVAR Name',
                        'CLINVAR Preferred Name','CLINVAR Preferred Name', 'OMIM IDs', 'MONDO IDs', 'MEDGEN IDs',
                        'Orphanet IDs', 'Human Phenotype Ontology']

                variantInfo[query] = {}
                for key in keys:
                    variantInfo[query][key] = ""

                variantInfo[query]['rs'] = rs
                variantInfo[query]['Freqs'] = freqs
                variantInfo[query]['Alleles'] = alleles
                variantInfo[query]['Ref'] = ref
                variantInfo[query]['Chrom'] = chrom
                variantInfo[query]['Pos hg19'] = posHG19
                variantInfo[query]['Pos hg38'] = posHG38

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
        else:
            notFound.append(var)

    print(f"\t#Not Found {len(notFound)}")
    return variantInfo, notFound

def searchSNPs(variantInfo, toSearchDict, folderAnot):
    print("Get annotation")

    fields = ['clinvar.rcv.conditions.identifiers', 'clinvar.rcv.conditions.name', 'clinvar.rcv.preferred_name',
              'dbsnp.alleles', 'dbsnp.rsid', 'dbsnp.ref', 'cadd.gene.gene_id', 'cadd.gene.genename', 'cadd.annotype',
              'cadd.consdetail', 'cadd.consequence', 'cadd.consscore', 'cadd.sift.cat', 'cadd.polyphen.cat','cadd.phred',
              'clinvar.hg38.start', 'clinvar.hg19.start', 'clinvar.chrom']

    toSearch = list(toSearchDict.keys())
    variantInfo, notFound = makeSearch(toSearch, variantInfo, fields)

    return variantInfo, notFound

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
        if line.strip() != "":
            split = line.strip().split("\t")
            dictCorrespondence[split[0]] = split[1]
    return dictCorrespondence

def createAnnotationTables(line, folderAnot):
    headerIDsTemp = line.strip().split()
    headerIDs = []

    dictOut = {}

    for i in range(len(headerIDsTemp)):
        ID = headerIDsTemp[i]
        if i < 9:
            headerIDs.append(ID)
        else:

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


def readIlluminaSentrix(fileName):
    file = open(fileName)

    dictReturn = {}

    header = True
    for line in file:
        if header:
            if "Sample_ID" in line:
                header = False
        else:
            if line != "" and line != " " and line != "\n":
                illuminaID, sentrixBarcode, sentrixPosition = line.strip().split(",")[0:3]
                sentrixID = f"{sentrixBarcode}_{sentrixPosition}"

                dictReturn[sentrixID] = illuminaID

    return dictReturn

def readIlluminaExternal(fileName):
    file = open(fileName)

    dictReturn = {}

    header = True
    for line in file:
        if header:
            if "IlluminaID" in line:
                header = False
        else:
            if line != "" and line != " " and line != "\n":
                illuminaID, externalID = line.strip().split("\t")
                dictReturn[illuminaID] = externalID.replace(" ", "_").replace("#N/A", "NonLARGE")

    return dictReturn

def extractGenes(genes, interval, VCF, vcfFolder, outputName, bcftools, batchList, logFile):
    dictGenes = {}
    file = open(genes)
    fileToExtract = open(f"{vcfFolder}/toExtract.txt", "w")


    for line in file:
        geneName, chrom, begin, end = line.strip().split()
        beginInterval = int(begin) - int(interval)
        endInterval = int(end) + int(interval)
        fileToExtract.write(f"{chrom}\t{beginInterval}\t{endInterval}\n")
        dictGenes[geneName]= {}
        dictGenes[geneName]["Begin"] = beginInterval
        dictGenes[geneName]["End"] = endInterval
        dictGenes[geneName]["Chrom"] = int(chrom)

    file.close()
    fileToExtract.close()

    file = open(batchList)
    fileChangeName = open(f"{vcfFolder}/toChangeID.txt", "w")
    externalList = []
    for line in file:
        IlluminaExternal, IlluminaSentrix, folder, ID = line.strip().split()
        dictIlluminaExternal = readIlluminaExternal(IlluminaExternal)
        dictIllumiaSentrix = readIlluminaSentrix(IlluminaSentrix)

        for sentrixID in dictIllumiaSentrix:
            illuminaID = dictIllumiaSentrix[sentrixID]
            externalID = dictIlluminaExternal[illuminaID]
            changed = False
            if externalID not in externalList:
                externalList.append(externalID)
            else:
                changed = True


                count = 1
                externalIDBase = externalID
                externalID = f"{externalIDBase}_{count}"

                while externalID in externalList:
                    count = count+1
                    externalID = f"{externalIDBase}_{count}"

                externalList.append(externalID)

            if changed:
                print(f"The external ID {externalIDBase} ({illuminaID}/{sentrixID}) was changed to {externalID}")
            else:
                print(f"We used external ID {externalID} ({illuminaID}/{sentrixID})")
            fileChangeName.write(f"{sentrixID} {externalID}\n")

    fileChangeName.close()

    execute(f"{bcftools} view -R {vcfFolder}/toExtract.txt -Oz -o {vcfFolder}/{outputName}_extracted.vcf.gz "
            f"{VCF}", logFile)

    execute(f"{bcftools} reheader -s {vcfFolder}/toChangeID.txt -o {vcfFolder}/{outputName}_externalID.vcf.gz "
            f"{vcfFolder}/{outputName}_extracted.vcf.gz", logFile)

    return f"{vcfFolder}/{outputName}_externalID.vcf.gz", dictGenes


def openVCFAndSearch(vcfName, variantInfo, folderAnot, name, dictGenes, dictFilter, logFile):
    toSearch = readVCFs(vcfName, variantInfo)
    variantInfo, notFound = searchSNPs(variantInfo, toSearch, folderAnot)
    saveSearch(f"{folderAnot}/{name}_AnnotTable.tsv")

    fileNot = open(f"{folderAnot}/{name}_NotFound.tsv", "w")
    for var in notFound:
        fileNot.write(f"{var}\n")
    fileNot.close()


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
                dictOut, headerIDs = createAnnotationTables(line, folderAnot)
                header = False
        else:
            split = line.strip().split()
            prefix = f"chr{split[0]}:g.{split[1]}"
            key = f"{prefix}{split[3]}>{split[4]}"

            A0 = split[3]
            A1 = split[4]
            if key in variantInfo:
                query = key
                turn = False
            else:
                key = f"{prefix}{split[4]}>{split[3]}"
                if key in variantInfo:
                    query = key
                    turn = True
                else:
                    query = "notFound"

            if query != "notFound":
                if not headerOnOutput:
                    headerOnOutput = True
                    for ind in dictOut:
                        dictOut[ind].write("Query\tID on VCF\tHomozygous\tGenotype")
                        for item in variantInfo[query]:
                            dictOut[ind].write(f"\t{item}")

                        dictOut[ind].write("\t\tFILTER\n")



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
                    #if maf < 0.05:
                    if phred > 15 and maf < 0.05:
                        IDVCF = split[2]

                        for i in range(9, len(split)):
                            ID = headerIDs[i]
                            information = split[i].split(":")
                            if not turn:
                                if "1" in information[0]:
                                    information[0] = information[0].replace("0", A0).replace("1", A1)

                                    alleles = information[0].split("/")
                                    print(alleles)
                                    hom = "FALSE"
                                    if alleles[0] == alleles[1]:
                                        hom = "TRUE"

                                    dictOut[ID].write(f"{query}\t{IDVCF}\t{hom}\t{information[0]}")
                                    for item in variantInfo[query]:
                                        dictOut[ID].write(f"\t{variantInfo[query][item]}")
                                    dictOut[ID].write(f"\t\t{dictFilter[IDVCF]}\n")
                            else:
                                if "0" in information[0]:
                                    information[0] = information[0].replace("1", A1).replace("0", A0)
                                    dictOut[ID].write(f"{query}\t{IDVCF}\t{information[0]}")
                                    for item in variantInfo[query]:
                                        dictOut[ID].write(f"\t{variantInfo[query][item]}")
                                    dictOut[ID].write(f"\t\t{dictFilter[IDVCF]}\n")
    for ID in dictOut:
        dictOut[ID].close()


def createExcelOnToSendFolder(TSV, folder, ind):
    #code adapted from  https://www.geeksforgeeks.org/convert-a-tsv-file-to-excel-using-python/

    workbook = Workbook(f"{folder}/{ind}.xlsx")
    worksheet = workbook.add_worksheet()

    readTSV = csv.reader(open(TSV, 'r', encoding='utf-8'), delimiter='\t')
    for row, data in enumerate(readTSV):
        worksheet.write_row(row, 0, data)
    workbook.close()



def makeSummary(folderAnot, logFile):
    #Prepare to send -> Summary table

    toSendFolder =  f"{folderAnot}/toSend/"
    execute(f"mkdir {toSendFolder}", logFile, True)
    execute(f"wc -l {folderAnot}/*.tsv > {folderAnot}/myTSV.txt", logFile, True)

    summaryDict = {}
    idRS = {}


    fileWC = open(f"{folderAnot}/myTSV.txt")
    for line in fileWC:
        split = line.strip().split()
        if 'Annot_' in split[1]:
            if int(split[0]) > 1: #Have at least one variant
                info = split[1].split("Annot_")
                if "_" in info[1]:
                    pop = info[1].split("_")[0]
                else:
                    pop = info[1][0:2]

                if pop not in summaryDict:
                    execute(f"mkdir {toSendFolder}/{pop}", logFile, True)
                    summaryDict[pop] = {}


                #execute(f"cp {split[1]} {toSendFolder}/{pop}", logFile, True)
                ind = info[1].replace(".tsv", "")
                createExcelOnToSendFolder(split[1], f"{toSendFolder}/{pop}", ind)


                idRS[ind] = {}

                fileTSV = open(split[1])
                header = True

                for line in fileTSV:
                    if header:
                        header = False
                    else:
                        split = line.strip().split("\t")
                        rs = split[3]

                        fields =  re.split(r'(\d+)', split[0])
                        chrom = int(fields[1])
                        pos = int(fields[3])


                        if rs not in idRS[ind]:
                            idRS[ind][rs] = ""
                            for gene in dictGenes:
                                if dictGenes[gene]["Chrom"] == chrom:
                                    if dictGenes[gene]["Begin"] <= pos <= dictGenes[gene]["End"]:
                                        if ind not in summaryDict[pop]:
                                            summaryDict[pop][ind] = {}

                                        if gene not in summaryDict[pop][ind]:
                                            summaryDict[pop][ind][gene] = 0
                                        summaryDict[pop][ind][gene] = summaryDict[pop][ind][gene]+1

    print(f"Creating the summary table")
    for pop in summaryDict:
        summaryTable = open(f"{toSendFolder}/{pop}/Summary.tsv", 'w')
        summaryTable.write("ID")
        for gene in dictGenes:
            summaryTable.write(f"\t{gene}")
        summaryTable.write("\n")

        for ind in summaryDict[pop]:
            summaryTable.write(f"{ind}")
            for gene in dictGenes:
                if gene in summaryDict[pop][ind]:
                    summaryTable.write(f"\t{summaryDict[pop][ind][gene]}")
                else:
                    summaryTable.write(f"\t0")
            summaryTable.write(f"\n")
        summaryTable.close()

        createExcelOnToSendFolder(f"{toSendFolder}/{pop}/Summary.tsv", f"{toSendFolder}/{pop}/", "Summary")

        execute(f"rm {toSendFolder}/{pop}/Summary.tsv", logFile, True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='From Illumina to annotation')

    requiredGeneral = parser.add_argument_group("Required arguments for all steps")
    requiredGeneral.add_argument('-o', '--outputFolder', help='Name of output folder', required=True)
    requiredGeneral.add_argument('-O', '--outputName', help='Name of output name', required=True)

    requiredCalling = parser.add_argument_group("Required arguments for variant calling")
    requiredCalling.add_argument('-b', '--bpm', help='BPM File from Illumina', required=True)
    requiredCalling.add_argument('-c', '--csv', help='CSV manifest file', required=True)
    requiredCalling.add_argument('-e', '--egt', help='EGT File from Illumina', required=True)
    #requiredCalling.add_argument('-f', '--folder', help='Folder with files from illumina', required=True)
    requiredCalling.add_argument('-r', '--genomeReference', help='Human genome reference', required=True)
    requiredCalling.add_argument('-L', '--batchList',
                                 help='File with four columns: \"IlluminaID_ExternalID, SampleID_SentrixBarcode_SentrixPosition,'
                                      'Path to RawFiles,ID\" separated by tab', required=True)

    requiredAnnotation = parser.add_argument_group("Required arguments for Annotation")
    requiredAnnotation.add_argument('-g', '--genes', help='File with gene list to be annotated. Four columns required '
                                                       '(separated by tab):Gene Name, Chromosome, Begin of the gene (bp), '
                                                       'End of the gene (bp).', required=True)
    requiredAnnotation.add_argument('-i', '--interval', help='size (in bp) of the region flanking the genes to be '
                                                             'included. Default = 1000', required=True, default = 1000)

    programs = parser.add_argument_group("Programs")
    programs.add_argument('-B', '--bcftools', help='bcftools', required=False, default="bcftools")
    programs.add_argument('-I', '--iaap', help='Illumina Array Analysis Platform path ', required=True)
    programs.add_argument('-P', '--plink2', help='Plink2 path ', required=True)

    optional = parser.add_argument_group("Optional arguments")
    optional.add_argument('-t', '--threads', help='Number of threads', default=96, type=int, required=False)
    optional.add_argument('-q', '--qc', help='QC parameters to VCF', required=False, default="")
    optional.add_argument('-p', '--previousSearch', help='File with MyVariant previous search', required=False,
                          default="")
    optional.add_argument('-C', '--correspondenceList', help='List with the correspondence Illumina to patientes ID',
                          required=False, default="")
    args = parser.parse_args()

    folder = args.outputFolder
    os.system(f"mkdir {folder}")
    logFile = open(f"{args.outputFolder}/{args.outputName}.log", 'w')
    #folder = createFolder(args.outputFolder, logFile)



    folderGTC = createFolder(folder+"/GTC/", logFile)
    generateGTC(args.iaap, args.bpm, args.egt, folder, folderGTC, args.outputName, args.batchList, args.threads, logFile)
    #generatePED(args.iaap, args.bpm, args.egt, args.folder, folderGTC, args.threads)

    allGTCFolder = createFolder(folder+"/GTCAll/", logFile)
    vcfFolder = createFolder(folder+"/VCFs/", logFile)
    VCF = gtc2VCF(args.bcftools, args.bpm, args.egt, args.csv, folderGTC, allGTCFolder, args.genomeReference, vcfFolder, args.threads, args.outputName, logFile)
    parametersQC = readFileQC(args.qc)

    #Change 9/13 -> We do not exclude based on VCF parameters
    VCFQC = basicVCFQC(VCF, parametersQC, args.bcftools, vcfFolder, args.outputName, logFile)


    VCFRegion, dictGenes = extractGenes(args.genes, args.interval, VCF, vcfFolder, args.outputName, args.bcftools,
                                        args.batchList, logFile)
    anotFolder = createFolder(folder + "/Anot/", logFile)
    infoQC = basicVCFQCToIncludeAsInfo(VCFRegion, parametersQC, anotFolder)

    VCFFinal = basicGenotypingQC(VCFRegion, folder, args.outputName, args.plink2, logFile)




    start = time.time()

    variantInfo = {}
    if args.previousSearch != "":
        variantInfo = readPreviousSearchFile(args.previousSearch)

    variantInfo = openVCFAndSearch(VCFFinal, variantInfo, anotFolder, args.outputName, dictGenes, infoQC, logFile)
    makeSummary(anotFolder, logFile)

    end = time.time()
    diff = end - start



    logFile.write(f'Annotation using MyVariants\nTime: {diff} seconds\n\n')
    logFile.close()