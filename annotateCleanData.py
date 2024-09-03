from datetime import date
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
scriptVersion = "Version 2"

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


def checkFam(famFile, allSamples, label, indQC):
    someSamples = []
    fam = open(famFile)
    for line in fam:
        split = line.strip().split()
        someSamples.append(split[1])
    fam.close()

    for ind in allSamples:
        if ind not in someSamples:
            if ind in infoQC:
                indQC[ind] = infoQC[ind] + f",{label}"
            else:
                indQC[ind] = f"{label}"
    return indQC

def checkBim(bimFile, allVariants, label, infoQC):
    someVariants = []
    bim = open(bimFile)
    for line in bim:
        split = line.strip().split()
        someVariants.append(split[1])
    bim.close()

    for variant in allVariants:
        if variant not in someVariants:
            if variant in infoQC:
                infoQC[variant] = infoQC[variant] + f", {label}"
            else:
                infoQC[variant] = infoQC[variant] + f"{label}"
    return infoQC


def checkPlinkFiles(prefixName,  allVariants, allSamples, label, infoQC, indQC):
    infoQC = checkBim(f"{prefixName}.bim", allVariants, label, infoQC)
    indQC = checkFam(f"{prefixName}.fam", allSamples, label, indQC)

    return infoQC, indQC

def basicGenotypingQC(VCF, folder, name, plink2, infoQC, logFile):
    indQC ={}
    folderPGEN = f"{folder}/PGEN"
    execute(f"mkdir {folderPGEN}", logFile)

    #Convert
    execute(f"{plink2} --vcf {VCF} --make-bed --out {folderPGEN}/{name}", logFile)

    allVariants = []
    allSamples = []

    bim = open(f"{folderPGEN}/{name}.bim")
    for line in bim:
        split = line.strip().split()
        allVariants.append(split[1])
    bim.close()

    fam = open(f"{folderPGEN}/{name}.fam")
    for line in fam:
        split = line.strip().split()
        allSamples.append(split[1])
    fam.close()


    #Remove missing
    execute(f"{plink2} --bfile {folderPGEN}/{name} --geno 0.05 --mind 0.05 --make-bed --out {folderPGEN}/{name}_missing",
            logFile)


    infoQC = checkBim(f"{folderPGEN}/{name}_missing.bim", allVariants,"Geno", infoQC)
    indQC = checkFam(f"{folderPGEN}/{name}_missing.fam", allSamples, "Mind", indQC)

    #Remove big deviation on HWE
    execute(f"{plink2} --bfile {folderPGEN}/{name}_missing --hwe 1e-10 --make-bed --out {folderPGEN}/{name}_HWE",
            logFile)

    infoQC = checkBim(f"{folderPGEN}/{name}_HWE.bim", allVariants,"HWE", infoQC)

    #PGEN to VCF
    #execute(f"{plink2} --pfile {folderPGEN}/{name}_HWE --recode vcf-iid --out {folderPGEN}/{name}_QC", logFile)

    return infoQC, indQC

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
            toReturn = (';').join(data)
    else:
        toReturn = data

    return toReturn

def readVCFs(vcfName, variantInfo, infoQC):

    print(f"Open VCF file ({vcfName})")
    toReturn = []
    searches = {}
    if vcfName[-2:] == "gz":
        file = gzip.open(vcfName)
    else:
        file = open(vcfName)

    header = True
    countToSearch = 0
    count = 0
    for line in file:
        if vcfName[-2:] == "gz":
            line = line.decode("utf-8")

        if header:
            if "#CHROM" in line:
                header = False
        else:
            split = line.strip().split()

            if infoQC[split[2]] == "PASS":

                if len(split[3]) == 1 and len(split[4]) == 1:
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
                else:
                    print(f"We will not look the {prefix} because {split[3]} or {split[4]} is indel")
            else:
                print(f"Remove from search because {split[2]} failed on {infoQC[split[2]]}")
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
        #print(dictData['gene'])
        if isinstance(dictData['gene'], list):
            for i in range(0, len(dictData['gene'])):
                if 'gene_id' in dictData['gene'][i]:
                    if geneID == "":
                        geneID = dictData['gene'][i]['gene_id']
                        geneName = dictData['gene'][i]['genename']
                    else:
                        geneID = f"{geneID};{dictData['gene'][i]['gene_id']}"
                        geneName = f"{geneName};{dictData['gene'][i]['genename']}"
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
                polyphen = f"{polyphen};{dictData['polyphen'][i]['cat']}"
        else:
            polyphen = dictData['polyphen']['cat']

    if 'sift' in dictData:
        if isinstance(dictData['sift'], list):
            sift = dictData['sift'][0]['cat']
            for i in range(len(dictData['sift'])):
                sift = f"{sift};{dictData['sift'][i]['cat']}"
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
                    name = f"{name};{data[i]['name']}"
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
    clinicalSignificance = ""
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
                name = f"{name};{returnedName}"

            if 'preferred_name' in data[i]:
                if preferredName == "":
                    preferredName = data[i]['preferred_name']
                else:
                    preferredName = f"{preferredName};{data[i]['preferred_name']}"
            if 'clinical_significance' in data[i]:
                if clinicalSignificance == "":
                    clinicalSignificance = data[i]['clinical_significance']
                else:
                    clinicalSignificance = clinicalSignificance+";"+data[i]['clinical_significance']
            else:
                if clinicalSignificance == "":
                    clinicalSignificance=";"
                else:
                    clinicalSignificance = clinicalSignificance+";"

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
        if 'clinical_significance' in data:
            clinicalSignificance = data['clinical_significance']

    return name, preferredName, clinicalSignificance, omimList, mondoList, medgenList, orphanetList, HPOList

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

def checkAndGetData(data, firstKey, secondKey):
    if firstKey in data:
        if secondKey == '':
            toReturn = data[firstKey]
        elif secondKey in data[firstKey]:
            toReturn = data[firstKey][secondKey]
        else:
            toReturn = ""
    else:
        toReturn = ""

    return getDataToPrint(toReturn)

def checkDataSNPEff(data):
    effectsSNPEff = putativeImpactSNPEff = featureIDSNPEff = lofGeneName = lofPercent = lofNumberTranscriptGene = ""
    if "snpeff" in data:
        if "ann" in data["snpeff"]:
            if isinstance(data["snpeff"]["ann"], list):
                for i in range(len(data["snpeff"]["ann"])):
                    if effectsSNPEff == "":
                        effectsSNPEff =data["snpeff"]["ann"][i]["effect"]
                    else:
                        effectsSNPEff = effectsSNPEff+";"+data["snpeff"]["ann"][i]["effect"]

                    if putativeImpactSNPEff == "":
                        putativeImpactSNPEff = data["snpeff"]["ann"][i]["putative_impact"]
                    else:
                        putativeImpactSNPEff = putativeImpactSNPEff+";"+data["snpeff"]["ann"][i]["putative_impact"]

                    if featureIDSNPEff == "":
                        featureIDSNPEff = data["snpeff"]["ann"][i]["feature_id"]
                    else:
                        featureIDSNPEff = featureIDSNPEff+";"+data["snpeff"]["ann"][i]["feature_id"]
            else:
                effectsSNPEff = data["snpeff"]["ann"]["effect"]
                putativeImpactSNPEff = data["snpeff"]["ann"]["putative_impact"]
                featureIDSNPEff = data["snpeff"]["ann"]["feature_id"]

        if "lof" in data["snpeff"]:
            if isinstance(data["snpeff"]["lof"], list):
                for i in range(len(data["snpeff"]["lof"])):
                    if lofGeneName == "":
                        lofGeneName = data["snpeff"]["lof"][i]["genename"]
                    else:
                        lofGeneName = lofGeneName + ";" + data["snpeff"]["lof"][i]["genename"]

                    if lofNumberTranscriptGene == "":
                        lofNumberTranscriptGene = data["snpeff"]["lof"][i]["number_of_transcripts_in_gene"]
                    else:
                        lofNumberTranscriptGene = lofNumberTranscriptGene+ ";"+ data["snpeff"]["lof"][i]["number_of_transcripts_in_gene"]

                    if lofPercent == "":
                        lofPercent = data["snpeff"]["lof"][i]["percent_of_transcripts_affected"]
                    else:
                        lofPercent = lofPercent + ";" + data["snpeff"]["lof"][i]["percent_of_transcripts_affected"]

    return effectsSNPEff, putativeImpactSNPEff, featureIDSNPEff, lofGeneName, lofPercent, lofNumberTranscriptGene



def dbnsfpData(data):
    #Basic info
    #print(data)
    chrom = hg19 = hg38 = bayesdelPred = clinpred = clinvarSignificance = clinvarID = fathmmXf = caddPred = \
        polyphen2HVARPred = polyphen2HDIVPred = metaSVMPred = primateai = sift = freqs1000G = ""

    chrom = checkAndGetData(data, 'chrom', '')
    hg19 = checkAndGetData(data, 'hg19', 'start')
    hg38 = checkAndGetData(data, 'hg38', 'start')
    freqs1000G = checkAndGetData(data, '1000gp3', 'af')

    if 'bayesdel' in data:
        bayesdelPred = checkAndGetData(data['bayesdel'], 'no_af','pred')
    clinpred = checkAndGetData(data, 'clinpred', 'pred')
    clinvarSignificance = checkAndGetData(data, 'clinvar', 'clinsig')
    clinvarID = checkAndGetData(data, 'clinvar', 'clinvar_id')
    fathmmXf = checkAndGetData(data, 'fathmm-xf', 'coding_pred')
    caddPred = checkAndGetData(data, 'cadd', 'pred')
    if 'polyphen2' in data:
        polyphen2HDIVPred = checkAndGetData(data['polyphen2'], 'hdiv', 'pred')
        polyphen2HVARPred = checkAndGetData(data['polyphen2'], 'hvar', 'pred')
    metaSVMPred = checkAndGetData(data, 'metasvm', 'pred')
    primateai = checkAndGetData(data, 'primateai', 'pred')
    sift = checkAndGetData(data, 'sift', 'pred')

    return chrom, hg19, hg38, bayesdelPred, clinpred, clinvarSignificance, clinvarID, fathmmXf, caddPred, polyphen2HVARPred, polyphen2HDIVPred, metaSVMPred, primateai, sift, freqs1000G


def checkDataDBNSFP(data):
    chrom = freqs1000G = hg19 = hg38 = bayesdelPred = clinpred = clinvarSignificance = clinvarID = fathmmXf = caddPred = polyphen2HVARPred = polyphen2HDIVPred = metaSVMPred = primateai = sift = ""

    if 'dbnsfp' in data:
        return dbnsfpData(data['dbnsfp'])

    return chrom, hg19, hg38, bayesdelPred, clinpred, clinvarSignificance, clinvarID, fathmmXf, caddPred, polyphen2HVARPred, polyphen2HDIVPred, metaSVMPred, primateai, sift, freqs1000G

def checkDataCadd(data):
    annotype, consdetail, consequence, consscore, geneID, geneName, polyphen, sift, phred = "", "", "", "", "", "", "", "", ""
    if 'cadd' in data:
        annotype, consdetail, consequence, consscore, geneID, geneName, polyphen, sift, phred = caddData(data["cadd"])

    annotypeOut = getDataToPrint(annotype)
    consdetailOut = getDataToPrint(consdetail)
    consequenceOut = getDataToPrint(consequence)
    consscoreOut = getDataToPrint(consscore)

    return annotypeOut, consdetailOut, consequenceOut, consscoreOut, geneID, geneName, polyphen, sift, phred


def checkDataDbSNP(data):
    if 'dbsnp' in data:
        return dbsnpData(data["dbsnp"])
    return "NA", "NA", "NA", "NA"

def checkDataClinvar(data):
    name, preferredName, clinicalSignificance, omimList, mondoList, medgenList, orphanetList, humanPhenotypeOntologyList = "", "", "", [], [], [], [], []
    if 'clinvar' in data:
        name, preferredName, clinicalSignificance, omimList, mondoList, medgenList, orphanetList, humanPhenotypeOntologyList = getRCVData(data["clinvar"]["rcv"])

    omimOut = getDataToPrint(omimList)
    mondoOut = getDataToPrint(mondoList)
    medgenOut = getDataToPrint(medgenList)
    orphanetOut = getDataToPrint(orphanetList)
    humanPhenotypeOntologyOut = getDataToPrint(humanPhenotypeOntologyList)

    return name, preferredName, clinicalSignificance, omimOut, mondoOut, medgenOut, orphanetOut, humanPhenotypeOntologyOut

def readPreviousSearch(pathPrevious):
    file = open(pathPrevious)
    variantInfo = {}

    header = True
    for line in file:
        if header:
            keys = line.strip().split("\t")
            header = False
        else:

            split = line.split("\t")

            variantInfo[split[0]] = {}
            for i in range(1,len(keys)):
                variantInfo[split[0]][keys[i]] = split[i].strip()

    return variantInfo, []

def makeSearch(toSearch, variantInfo, fields, dictGenes, pathPrevious):
    keys = ['rs', 'Alleles', 'Ref', 'Chrom', 'Pos hg19', 'Pos hg38', 'Freqs DBSNP', 'Freqs 1000G Phase 3', 'CADD Phred', 'CADD Annotype',
            'CADD Consequence Detail', 'CADD Consequence', 'CADD Consequence Score', 'CADD Gene ID',
            'CADD Gene Name', 'Sift Cathegory', 'Polyphen Cathegory', 'CLINVAR Name',
            'CLINVAR Preferred Name', 'CLINVAR Clinical Significance', 'OMIM IDs', 'MONDO IDs', 'MEDGEN IDs',
            'Orphanet IDs', 'Human Phenotype Ontology']
    pred = ["Bayesdel", "Clinpred", "CLINVAR Significance", "CLINVAR ID", "FATHMM-XF", "Polyphen2 HVAR",
            "Polyphen2 HDIV", "PrimateAI", "SIFT", "DBNSFP SVM", "SnpEff Effect", "SnpEff Feature ID",
            "SnpEff Putative Impact", "SnpEff LoF Gene", "SnpEff Number of Transcripts in Gene",
            "SnpEff Percent Of Transcripts Affected"]

    if pathPrevious:
        variantInfo, notFound = readPreviousSearch(pathPrevious)
        return variantInfo, notFound, pred


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

                annotypeOut, consdetailOut, consequenceOut, consscoreOut, geneID, geneName, polyphen, sift, phred = checkDataCadd(data)

                alleles, freqs, rs, ref = checkDataDbSNP(data)

                name, preferredName,clinicalSignificance, omimOut, mondoOut, medgenOut, orphanetOut, humanPhenotypeOntologyOut = checkDataClinvar(data)

                chrom, posHG19, posHG38, bayesdelPred, clinpred, clinvarSignificance, clinvarID, fathmmXf, caddPred, \
                polyphen2HVARPred, polyphen2HDIVPred, metaSVMPred, primateai, sift, freqs1000G = checkDataDBNSFP(data)

                effectsSNPEff, putativeImpactSNPEff, featureIDSNPEff, lofGeneName, lofPercent, lofNumberTranscriptGene = checkDataSNPEff(data)

                dataSplit = re.split(r'(\d+)', query)

                gene, inOrOut = getGeneBasedInPos(int(dataSplit[1]), int(dataSplit[3]), dictGenes)

                if gene != "Error" and inOrOut != "Error":
                    variantInfo[query] = {}
                    for key in keys:
                        variantInfo[query][key] = ""

                    for key in pred:
                        variantInfo[query][key] = ""

                    variantInfo[query]['rs'] = rs
                    variantInfo[query]['Freqs DBSNP'] = freqs
                    variantInfo[query]['Freqs 1000G Phase 3'] = freqs1000G
                    variantInfo[query]['Alleles'] = alleles
                    variantInfo[query]['Ref'] = ref
                    variantInfo[query]['Chrom'] = chrom
                    variantInfo[query]['Pos hg19'] = posHG19
                    variantInfo[query]['Pos hg38'] = posHG38

                    variantInfo[query]['Gene'] = gene
                    variantInfo[query]['Inside/Outside'] = inOrOut

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
                    variantInfo[query]['CLINVAR Clinical Significance'] = clinicalSignificance
                    if clinvarSignificance != "":
                        variantInfo[query]['CLINVAR Clinical Significance'] = variantInfo[query]['CLINVAR Clinical Significance']+f";{clinvarSignificance}"
                        variantInfo[query]['CLINVAR Preferred Name'] = variantInfo[query]['CLINVAR Preferred Name']+f";{clinvarID}"

                    variantInfo[query]['OMIM IDs'] = omimOut
                    variantInfo[query]['MONDO IDs'] = mondoOut
                    variantInfo[query]['MEDGEN IDs'] = medgenOut
                    variantInfo[query]['Orphanet IDs'] = orphanetOut
                    variantInfo[query]['Human Phenotype Ontology'] = humanPhenotypeOntologyOut

                    variantInfo[query]['Bayesdel'] = bayesdelPred
                    variantInfo[query]['Clinpred'] = clinpred
                    variantInfo[query]['CLINVAR Significance'] = clinvarSignificance
                    variantInfo[query]['CLINVAR ID'] = clinvarID
                    variantInfo[query]["FATHMM-XF"] = fathmmXf
                    #variantInfo[query]["CADD"] = caddPred
                    variantInfo[query]["Polyphen2 HVAR"] = polyphen2HVARPred
                    variantInfo[query]["Polyphen2 HDIV"] = polyphen2HDIVPred
                    variantInfo[query]["PrimateAI"] = primateai
                    variantInfo[query]["SIFT"] = sift
                    variantInfo[query]["DBNSFP SVM"] = metaSVMPred
                    variantInfo[query]["SnpEff Effect"] = effectsSNPEff
                    variantInfo[query]["SnpEff Feature ID"] = featureIDSNPEff
                    variantInfo[query]["SnpEff Putative Impact"] = putativeImpactSNPEff
                    variantInfo[query]["SnpEff LoF Gene"] = lofGeneName
                    variantInfo[query]["SnpEff Number of Transcripts in Gene"] = lofNumberTranscriptGene
                    variantInfo[query]["SnpEff Percent Of Transcripts Affected"] = lofPercent
                else:
                    print(f"-> {query} : {chrom}, {posHG19} -> REMOVED for window constraint")

        else:
            notFound.append(var)

    print(f"\t#Not Found {len(notFound)}")
    return variantInfo, notFound, pred

def searchSNPs(variantInfo, toSearchDict, folderAnot, dictGenes, previousSearch):
    print(" annotation")

    #Old version
    #fields = ['clinvar.rcv.conditions.identifiers', 'clinvar.rcv.conditions.name', 'clinvar.rcv.preferred_name',
    #          'dbsnp.alleles', 'dbsnp.rsid', 'dbsnp.ref', 'cadd.gene.gene_id', 'cadd.gene.genename', 'cadd.annotype',
    #          'cadd.consdetail', 'cadd.consequence', 'cadd.consscore', 'cadd.sift.cat', 'cadd.polyphen.cat','cadd.phred',
    #          'clinvar.hg38.start', 'clinvar.hg19.start', 'clinvar.chrom']

    #New version
    fields = ['clinvar', 'dbsnp', 'cadd', 'dbnsfp', 'snpeff']

    toSearch = list(toSearchDict.keys())
    variantInfo, notFound, pred = makeSearch(toSearch, variantInfo, fields, dictGenes, previousSearch)

    return variantInfo, notFound, pred

def isClinvarPathogenic(dictData):
    if 'CLINVAR Clinical Significance' in dictData:
        if "Pathogenic" in dictData['CLINVAR Clinical Significance'] or "Likely_pathogenic" in dictData['CLINVAR Clinical Significance'] or "risk_factor" in dictData['CLINVAR Clinical Significance']:
            #input(f"Returning TRUE {dictData['CLINVAR Significance']}")
            return True
    #input(f"Returning FALSE {dictData['CLINVAR Significance']}")
    return False

def isPathogenic(dictData):
    #D - deleterious, T - tolarated, N- neutral, B - benign, P - Pathogenic,  U - Uncertain
    predictors = ['Bayesdel', 'Clinpred', 'FATHMM-XF', 'Polyphen2 HVAR', 'Polyphen2 HDIV',
                  'PrimateAI', 'SIFT']

    inferences = ['D', 'P', 'N', 'T', 'U', 'B']
    dictPrediction = {}
    for infer in inferences:
        dictPrediction[infer] = 0
    dictPrediction['total'] = 0

    for predictor in predictors:
        data = ""
        if ";" in dictData[predictor]:
            dataList = dictData[predictor].strip().split(";")
            if 'D' in dataList:
                data = 'D'
            elif 'P' in dataList:
                data = 'P'
            else:
                data = dataList[0]
        else:
            data = dictData[predictor]


        if data != "":
            dictPrediction[data] = dictPrediction[data] + 1
            dictPrediction['total'] = dictPrediction['total']+1

    if dictPrediction['total']/2 <= (dictPrediction["D"]+dictPrediction["P"]) and dictPrediction['total'] > 0:
        return True
    else:
        return False




def isValidCADDPhredAndMAF(dictData, phredCADDMin, cutoffMAF):
    if 'CADD Phred' in dictData:
        if isNotNull(dictData['CADD Phred']):
            #print("notNull")
            phred = float(dictData['CADD Phred'])
            AFList, validAF = getAFs(dictData['Freqs DBSNP'])
            #print(f"AFList DBSNP : {AFList} -> {validAF}")
            if not validAF:
                AFList, validAF = getAFs(dictData['Freqs 1000G Phase 3'])
                #print(f"AFList 1000G : {AFList} -> {validAF}")
            if validAF:
                sortedAF = np.sort(AFList)
                #print(f"AF sorted: {sortedAF}")
                if len(AFList) == 1:
                    maf = sortedAF[0]
                elif len(AFList) == 2:
                    maf = sortedAF[0]
                elif len(AFList) == 3:
                    maf = sortedAF[1]
                else:
                    maf = sortedAF[2]

                #print(f"{phred} > {phredCADDMin} and {maf} < {cutoffMAF}")
                if phred > phredCADDMin and maf < cutoffMAF:
                    return True
                else:
                    return False
            else:
                return False
        else:
            return False
    else:
        return False

def isNotNull(data):
    if data != "NA" and data != "":
        return True
    return False

def getAFs(data):
    toReturn = []
    data = str(data)
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

def createAnnotationTables(line, folderAnot, variantInfo, pred):
    headerIDsTemp = line.strip().split()
    headerIDs = []

    dictOut = {}
    headerToInclude = ""
    for var in variantInfo:
        for item in variantInfo[var]:
            if item not in pred:
                headerToInclude = headerToInclude + f"\t{item}"
        break
    headerToInclude = headerToInclude+"\t\tFILTER\t"
    for item in pred:
        headerToInclude = headerToInclude + f"\t{item}"
    headerToInclude = headerToInclude + "\n"


    for i in range(len(headerIDsTemp)):
        ID = headerIDsTemp[i]
        if i < 9:
            headerIDs.append(ID)
        else:
            workID = ID
            #Create a unique file to each ID, even if the name is the same
            if ID not in dictOut:
                dictOut[ID] = open(f"{folderAnot}/Annot_{ID}.tsv", 'w')
                headerIDs.append(ID)
            else:
                count = 1
                while f"{ID}_{count}" in headerIDs:
                    count = count + 1

                workID = f"{ID}_{count}"
                headerIDs.append(workID)
                dictOut[workID] = open(f"{folderAnot}/Annot_{ID}_{count}.tsv", 'w')
            dictOut[workID].write(f"Accepted CADD+MAF\tAccepted Clinvar\tAccepted Predictors (>50%)\t\tQuery\tID on VCF\tHomozygous\tGenotype\t{headerToInclude}")


    return dictOut, headerIDs

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

def extractGenes(genes, interval, VCF, vcfFolder, outputName, bcftools, sampleSheet, logFile):
    dictGenes = {}
    file = open(genes)
    fileToExtract = open(f"{vcfFolder}/toExtract.txt", "w")


    for line in file:
        geneName, chrom, begin, end = line.strip().split()
        beginInterval = int(begin) - int(interval)
        endInterval = int(end) + int(interval)
        fileToExtract.write(f"{chrom}\t{beginInterval}\t{endInterval}\n")
        dictGenes[geneName]= {}
        dictGenes[geneName]["BeginOriginal"] = int(begin)
        dictGenes[geneName]["EndOriginal"] = int(end)
        dictGenes[geneName]["Begin"] = beginInterval
        dictGenes[geneName]["End"] = endInterval
        dictGenes[geneName]["Chrom"] = int(chrom)

    file.close()
    fileToExtract.close()

    file = open(sampleSheet)
    fileChangeName = open(f"{vcfFolder}/toChangeID.txt", "w")

    header = True
    for line in file:
        if header:
            if "Sample_ID,SentrixBarcode_A,SentrixPosition_A,Path" in line:
                header = False
        else:
            split = line.split(",")
            externalID = split[0]
            sentrixID = f"{split[1]}_{split[2]}"

            fileChangeName.write(f"{sentrixID} {externalID}\n")

    fileChangeName.close()

    # First extract
    execute(f"{bcftools} view -R {vcfFolder}/toExtract.txt -Oz -o {vcfFolder}/{outputName}_extracted.vcf.gz "
            f"{VCF}", logFile)
    execute(f"{bcftools} index {vcfFolder}/{outputName}_extracted.vcf.gz", logFile)

    # Change header
    execute(f"{bcftools} reheader -s {vcfFolder}/toChangeID.txt -o {vcfFolder}/{outputName}_externalID.vcf.gz "
            f"{vcfFolder}/{outputName}_extracted.vcf.gz", logFile)
    execute(f"{bcftools} index {vcfFolder}/{outputName}_externalID.vcf.gz", logFile)


    return f"{vcfFolder}/{outputName}_externalID.vcf.gz", dictGenes


def openVCFAndSearch(vcfName, variantInfo, folderAnot, name, infoQC, dictGenes, previousSearch = ""):
    toSearch = readVCFs(vcfName, variantInfo, infoQC)
    variantInfo, notFound, pred = searchSNPs(variantInfo, toSearch, folderAnot, dictGenes, previousSearch)
    saveSearch(f"{folderAnot}/{name}_AnnotTable.tsv")

    fileNot = open(f"{folderAnot}/{name}_NotFound.tsv", "w")
    for var in notFound:
        fileNot.write(f"{var}\n")
    fileNot.close()

    return variantInfo, pred

def readSampleSheet(inputFile):
    sampleSheetFile = open(inputFile)
    header = True
    dictSamples = {}
    for line in sampleSheetFile:
        if header:
            if "Sample_ID,SentrixBarcode_A,SentrixPosition_A,Path" in line:
                header = False
        else:
            infoID = line.strip().split(',')[0]
            split = infoID.split('_')


            for i in range(0, len(split)):
                if i == 0:
                    ID = split[i]
                else:
                    ID = ID + f"_{split[i]}"

            #ID and batch
            batch = split[-1]
            dictSamples[ID] = {}

            #Letter code -> get the ID without any underline
            if "_" in ID:
                countryID = ID.split("_")[0]
            else:
                countryID = ID

            #Special cases -> HIHG, Nuytemans, X code that are weird, one letter code
            if "HIHG" in ID:
                countryIDOut = "HIHG"
            elif "NonLARGE" in ID:
                countryIDOut = "NonLARGE"
            elif "X00" in ID:
                countryIDOut = "X00"
            elif len(countryID) > 2:
                countryIDOut = countryID[0:2]
            else:
                countryIDOut = countryID

            dictSamples[ID]["Code"] = countryIDOut
            dictSamples[ID]["Batch"] = batch

    return dictSamples

def checkIfVariantWasFound(line, variantInfo):
    split = line.strip().split()
    prefix = f"chr{split[0]}:g.{split[1]}"
    key = f"{prefix}{split[3]}>{split[4]}"

    if key in variantInfo:
        query = key
    else:
        key = f"{prefix}{split[4]}>{split[3]}"
        if key in variantInfo:
            query = key
            #input(f"Inverse -> {prefix}")
        else:
            query = "notFound"
    return query, split

def getGeneBasedInPos(chrom, pos, dictGenes):
    for gene in dictGenes:
        if dictGenes[gene]["Chrom"] == chrom:
            if dictGenes[gene]["Begin"] <= pos <= dictGenes[gene]["End"]:
                if dictGenes[gene]["BeginOriginal"] <= pos <= dictGenes[gene]["EndOriginal"]:
                    return gene, "Inside"
                else:
                    return gene, "Outside"
    return "Error", "Error"

def writeFiles(vcfName, variantInfo, pred, folderAnot, variantQC, indQC, dictGenes, sampleSheet, logFile, phredCADDMin = 15, cutoffMAF = 0.05):
    print("Reading sample sheet")
    dictSampleSheet = readSampleSheet(sampleSheet)

    dictSummary = {}

    print(f"Open VCF file again({vcfName}) ")

    if vcfName[-2:] == "gz":
        file = gzip.open(vcfName)
    else:
        file = open(vcfName)

    header = True

    for line in file:
        if vcfName[-2:] == "gz":
            line = line.decode("utf-8")

        if header:
            if "#CHROM" in line:
                dictOut, headerIDs = createAnnotationTables(line, folderAnot, variantInfo, pred)
                header = False
        else:
            query, split = checkIfVariantWasFound(line, variantInfo)

            chrom, pos, IDVCF, A0, A1 = split[0:5]

            if query != "notFound":

                CADDMAF = isValidCADDPhredAndMAF(variantInfo[query], float(phredCADDMin), float(cutoffMAF))
                #input(f"CADDMAF : {CADDMAF}")
                preditorsPathogenic = isPathogenic(variantInfo[query])
                clinvarPathogenic = isClinvarPathogenic(variantInfo[query])
                if CADDMAF or preditorsPathogenic or clinvarPathogenic:
                    #Individual table

                    for i in range(9, len(split)):
                        ID = headerIDs[i]
                        information = split[i].split(":")
                        GT = information[0]
                        #At least one risk allele
                        if "1" in GT:
                            GT = GT.replace("0", A0).replace("1", A1)

                            alleles = GT.split("/")
                            hom = "FALSE"
                            if alleles[0] == alleles[1]:
                                hom = "TRUE"

                            #Individual table
                            begin = f"{CADDMAF}\t{clinvarPathogenic}\t{preditorsPathogenic}\t"
                            dictOut[ID].write(f"{begin}\t{query}\t{IDVCF}\t{hom}\t{GT}\t")
                            for item in variantInfo[query]:
                                if item not in pred:
                                    dictOut[ID].write(f"\t{variantInfo[query][item]}")
                            dictOut[ID].write(f"\t\t{variantQC[IDVCF]}\t")

                            for item in pred:
                                dictOut[ID].write(f"\t{variantInfo[query][item]}")
                            dictOut[ID].write(f"\n")

                            pop = dictSampleSheet[ID]["Code"]

                            print(f"{variantQC[IDVCF]}, {clinvarPathogenic}, {CADDMAF}, {preditorsPathogenic}")
                            color = getColor(variantQC[IDVCF], clinvarPathogenic, CADDMAF, preditorsPathogenic)


                            gene,inOrOut = getGeneBasedInPos(int(chrom), int(pos), dictGenes)
                            if pop not in dictSummary:
                                dictSummary[pop] = {}

                            if ID not in dictSummary[pop]:
                                dictSummary[pop][ID] = {}

                            if gene != "Error":
                                if gene not in dictSummary[pop][ID]:
                                    dictSummary[pop][ID][gene] = {}
                                    dictSummary[pop][ID][gene]["red"] = []
                                    dictSummary[pop][ID][gene]["green"] = []
                                    dictSummary[pop][ID][gene]["blue"] = []
                                    dictSummary[pop][ID][gene]["grey"] = []

                                if pos not in dictSummary[pop][ID][gene][color]:
                                    dictSummary[pop][ID][gene][color].append(pos)
                            else:
                                print(f"Very weird {chrom} {pos} {dictGenes}")



    for ID in dictOut:
        dictOut[ID].close()

    today = date.today()

    execute(f"mkdir {folderAnot}/toSend/", logFile)

    for pop in dictSummary:
        toSendFolder = f"{folderAnot}/toSend/{pop}"
        execute(f"mkdir {toSendFolder}", logFile)

        fileSummary = open(f"{folderAnot}/{pop}_Summary.tsv", "w")
        fileSummary.write("SampleID\tBatch")
        for gene in dictGenes:
            fileSummary.write(f"\t{gene}\t{gene}_Green|Blue|Grey")
        fileSummary.write("\tQC\t\tAnnotationDate\tScriptVersion\n")

        for ind in dictSummary[pop]:
            batch = dictSampleSheet[ID]["Batch"]

            fileSummary.write(f"{ind}\t{batch}")

            toConvert = False
            for gene in dictGenes:
                if gene not in dictSummary[pop][ind]:
                    total = numVarR = numVarG = numVarB = numVarGrey = 0
                else:
                    total = 0
                    toConvert = True
                    #numVarR = len(dictSummary[pop][ind][gene]["red"])
                    numVarG = len(dictSummary[pop][ind][gene]["green"])
                    total = total + numVarG
                    numVarB = len(dictSummary[pop][ind][gene]["blue"])
                    total = total + numVarB
                    numVarGrey = len(dictSummary[pop][ind][gene]["grey"])
                    total = total + numVarGrey

                if total == 0:
                    total = ""

                fileSummary.write(f"\t{total}\t{numVarG}|{numVarB}|{numVarGrey}")
            if ind not in indQC:
                indInfo = ""
            else:
                indInfo = indQC[ind]
            fileSummary.write(f"\t{indInfo}\t\t{today}\t{scriptVersion}\n")

            if toConvert:
                createExcelOnToSendFolderData(f"{folderAnot}/Annot_{ind}.tsv", toSendFolder, ind)
        fileSummary.close()
        createExcelOnToSendFolder(f"{folderAnot}/{pop}_Summary.tsv", toSendFolder, f"{pop}_summary")

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


#=========================================================== Main ====================================================

def getColor(QC, Clinvar, CADD, InSilico):
    if QC != "PASS": #35
        return "red"
    elif Clinvar == True or Clinvar == "True" or Clinvar == "TRUE": #1
        return "green"
    elif CADD == True or CADD == "True" or CADD == "TRUE": #0
        return "blue"
    elif InSilico == True or InSilico == "True" or InSilico == "TRUE": #2
        return "grey"


def createExcelOnToSendFolderData(TSV, folder, ind):
    #code adapted from  https://www.geeksforgeeks.org/convert-a-tsv-file-to-excel-using-python/

    dictData = {}

    fileTSV = open(TSV)
    header = True
    for line in fileTSV:
        if header:
            headerLine = line.strip().split("\t")
            header = False
        else:
            split = line.strip().split("\t")
            ID = split[5]

            if ID not in dictData:
                dictData[ID] = {}
                for i in range(len(split)):
                    field = headerLine[i]
                    dictData[ID][field] = split[i]
            dictData[ID]["color"] = getColor(split[37], split[1], split[0], split[2])

    summaryFields = ["ID on VCF", "Query", "", "Gene", "Inside/Outside",
                     "Chrom", "Pos hg19", "Pos hg38", "",
                     "Alleles", "Ref", "Genotype", "",
                     "rs","Freqs DBSNP", "Freqs 1000G Phase 3","",
                     "FILTER", "",
                     "Accepted CADD+MAF", "Accepted Clinvar", "Accepted Predictors (>50%)"]
    infoFromCLINVAR = ["ID on VCF", "CLINVAR Name", "CLINVAR Preferred Name", "CLINVAR Clinical Significance", "",
                       "OMIM IDs","MONDO IDs", "MEDGEN IDs", "Orphanet IDs", "Human Phenotype Ontology"]
    infoFromCADD = ["ID on VCF", "CADD Phred", "CADD Annotype", "CADD Consequence Detail", "CADD Consequence",
                    "CADD Consequence Score", "CADD Gene ID", "CADD Gene Name", "",
                    "Sift Cathegory", "Polyphen Cathegory"]
    infoFromInSilico = ["ID on VCF", "CLINVAR Significance", "CLINVAR ID","",
                        "Bayesdel", "Clinpred", "FATHMM-XF", "Polyphen2 HVAR", "Polyphen2 HDIV", "PrimateAI", "SIFT", "",
                        "SnpEff Effect", "SnpEff Feature ID", "SnpEff Putative Impact", "SnpEff LoF Gene",
                        "SnpEff Number of Transcripts in Gene", "SnpEff Percent Of Transcripts Affected"]


    workbook = Workbook(f"{folder}/{ind}.xlsx")
    bold = workbook.add_format({'bold': True})

    red = workbook.add_format()
    red.set_font_color("#FF0000")

    green = workbook.add_format()
    green.set_bg_color("#49E468")

    blue = workbook.add_format()
    blue.set_bg_color("#31C0CF")

    grey = workbook.add_format()
    grey.set_bg_color("#B8BFBA")

    worksheetSummary = workbook.add_worksheet("Summary")

    line = 0
    for i in range(len(summaryFields)):
        if summaryFields[i] == "ID on VCF":
            worksheetSummary.write(line, i, "VCF ID", bold)
        else:
            worksheetSummary.write(line, i, summaryFields[i], bold)

    line = 1
    for ID in dictData:
        for j in range(len(summaryFields)):
            field = summaryFields[j]
            if field != "":
                if field in dictData[ID]:
                    toWrite = dictData[ID][field]
                else:
                    toWrite = ""

                if dictData[ID]["color"] == "red":
                    pass
                    #worksheetSummary.write(line, j, toWrite, red)
                elif dictData[ID]["color"] == "green":
                    worksheetSummary.write(line, j, toWrite, green)
                elif dictData[ID]["color"] == "blue":
                    worksheetSummary.write(line, j, toWrite, blue)
                elif dictData[ID]["color"] == "grey":
                    worksheetSummary.write(line, j, toWrite, grey)
        line = line+1


    worksheetCADD = workbook.add_worksheet("Info from CADD")

    line = 0
    for i in range(len(infoFromCADD)):
        if infoFromCADD[i] == "ID on VCF":
            worksheetCADD.write(line, i, "VCF ID", bold)
        else:
            worksheetCADD.write(line, i, infoFromCADD[i], bold)

    line = 1
    for ID in dictData:
        for j in range(len(infoFromCADD)):
            field = infoFromCADD[j]
            if field != "":
                if field in dictData[ID]:
                    toWrite = dictData[ID][field]
                else:
                    toWrite = ""

                if dictData[ID]["color"] == "red":
                    pass
                    #worksheetCADD.write(line, j, toWrite, red)
                elif dictData[ID]["color"] == "green":
                    worksheetCADD.write(line, j, toWrite, green)
                elif dictData[ID]["color"] == "blue":
                    worksheetCADD.write(line, j, toWrite, blue)
                elif dictData[ID]["color"] == "grey":
                    worksheetCADD.write(line, j, toWrite, grey)
        line = line + 1

    worksheetCLINVAR = workbook.add_worksheet("Info from CLINVAR")

    line = 0
    for i in range(len(infoFromCLINVAR)):
        if infoFromCLINVAR[i] == "ID on VCF":
            worksheetCLINVAR.write(line, i, "VCF ID", bold)
        else:
            worksheetCLINVAR.write(line, i, infoFromCLINVAR[i], bold)

    line = 1
    for ID in dictData:
        for j in range(len(infoFromCLINVAR)):
            field = infoFromCLINVAR[j]
            if field != "":
                if field in dictData[ID]:
                    toWrite = dictData[ID][field]
                else:
                    toWrite = ""

                if dictData[ID]["color"] == "red":
                    pass
                    #worksheetCLINVAR.write(line, j, toWrite, red)
                elif dictData[ID]["color"] == "green":
                    worksheetCLINVAR.write(line, j, toWrite, green)
                elif dictData[ID]["color"] == "blue":
                    worksheetCLINVAR.write(line, j, toWrite, blue)
                elif dictData[ID]["color"] == "grey":
                    worksheetCLINVAR.write(line, j, toWrite, grey)
        line = line + 1


    worksheetInSilico = workbook.add_worksheet("Info from In Silico Classifier")

    line = 0
    for i in range(len(infoFromInSilico)):
        if infoFromInSilico[i] == "ID on VCF":
            worksheetInSilico.write(line, i, "VCF ID", bold)
        else:
            worksheetInSilico.write(line, i, infoFromInSilico[i], bold)

    line = 1
    for ID in dictData:
        for j in range(len(infoFromInSilico)):
            field = infoFromInSilico[j]
            if field != "":
                if field in dictData[ID]:
                    toWrite = dictData[ID][field]
                else:
                    toWrite = ""

                if dictData[ID]["color"] == "red":
                    pass
                    #worksheetInSilico.write(line, j, toWrite, red)
                elif dictData[ID]["color"] == "green":
                    worksheetInSilico.write(line, j, toWrite, green)
                elif dictData[ID]["color"] == "blue":
                    worksheetInSilico.write(line, j, toWrite, blue)
                elif dictData[ID]["color"] == "grey":
                    worksheetInSilico.write(line, j, toWrite, grey)
        line = line + 1

    workbook.close()

def createExcelOnToSendFolder(TSV, folder, ind):
    #code adapted from  https://www.geeksforgeeks.org/convert-a-tsv-file-to-excel-using-python/

    workbook = Workbook(f"{folder}/{ind}.xlsx")
    worksheet = workbook.add_worksheet()

    readTSV = csv.reader(open(TSV, 'r', encoding='utf-8'), delimiter='\t')
    for row, data in enumerate(readTSV):
        worksheet.write_row(row, 0, data)
    workbook.close()



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Annotate VCF')
    requiredAnnotation = parser.add_argument_group("Required arguments for Annotation")
    requiredAnnotation.add_argument('-v', '--vcf', help='VCF to be annotated', required=True)
    requiredAnnotation.add_argument('-g', '--genes', help='File with gene list to be annotated. Four columns required '
                                                       '(separated by tab):Gene Name, Chromosome, Begin of the gene (bp), '
                                                       'End of the gene (bp).', required=True)
    requiredAnnotation.add_argument('-i', '--interval', help='size (in bp) of the region flanking the genes to be '
                                                             'included. Default = 1000', required=True, default = 1000)
    requiredAnnotation.add_argument('-o', '--outputName', help='Prefix of output files', required=True)
    requiredAnnotation.add_argument('-O', '--outputFolder', help='Prefix of output files', required=True)
    requiredAnnotation.add_argument('-s', '--sampleSheet', help='Sample sheet used in the iaap', required=True)


    optional = parser.add_argument_group("Optional arguments")
    optional.add_argument('-q', '--qc', help='QC parameters to VCF', required=False, default="")
    optional.add_argument('-p', '--previousSearch', help='File with MyVariant previous search', required=False,
                          default="")

    programs = parser.add_argument_group("Programs path")
    optional.add_argument('-B', '--bcftools', help='Path to bcftools (default: bcftools)', required=False,
                          default="bcftools")
    optional.add_argument('-P', '--plink2', help='Path to plink2 (default: plink2)', required=False,
                          default="plink2")

    args = parser.parse_args()

    os.system(f"mkdir {args.outputFolder}")
    logFile = open(f"{args.outputFolder}/{args.outputName}.log", "w")

    #Change 9/13 -> We do not exclude based on VCF parameters
    #VCFQC = basicVCFQC(VCF, parametersQC, args.bcftools, vcfFolder, args.outputName, logFile)
    parametersQC = readFileQC(args.qc)

    VCFRegion, dictGenes = extractGenes(args.genes, args.interval, args.vcf, args.outputFolder, args.outputName, args.bcftools,
                                        args.sampleSheet, logFile)
    anotFolder = createFolder(f"{args.outputFolder}/Anot/", logFile)
    infoQC = basicVCFQCToIncludeAsInfo(VCFRegion, parametersQC, anotFolder)

    infoQC, indQC = basicGenotypingQC(VCFRegion, args.outputFolder, args.outputName, args.plink2, infoQC, logFile)


    start = time.time()

    variantInfo = {}
    if args.previousSearch != "":
        variantInfo = readPreviousSearchFile(args.previousSearch)

    variantInfo, pred = openVCFAndSearch(VCFRegion, variantInfo, anotFolder, args.outputName, infoQC, dictGenes)
    writeFiles(VCFRegion, variantInfo, pred, anotFolder, infoQC, indQC, dictGenes, args.sampleSheet, logFile)

    end = time.time()
    diff = end - start

    logFile.write(f'Annotation using MyVariants\nTime: {diff} seconds\n\n')
    logFile.close()