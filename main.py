import myvariant
import argparse
import gzip
import os

mv = myvariant.MyVariantInfo()

#=========================================================== Utils ====================================================
def createFolder(folderName):
    line = f"mkdir {folderName}"
    execute(line)

    return folderName

def execute(commandLine, toRun = False):
    print(commandLine)
    if toRun:
        os.system(commandLine)

#=========================================================== QC ====================================================
def gtc2VCF(bcftools, bpm, egt, csv, folderGTC, allGTCFolder, genomeReference, vcfFolder, threads, outputName):
    for folder in os.listdir(folderGTC):
        line = f"cp {folderGTC}/{folder}/* {allGTCFolder}"
        execute(line)

    line = f"{bcftools} +gtc2vcf --no-version -Oz --bpm {bpm} --egt {egt} --csv {csv} --gtcs {allGTCFolder} " \
           f"--fasta-ref {genomeReference} --output {vcfFolder}/{outputName}.vcf.gz --threads {threads} " \
           f"--use-gtc-sample-names --adjust-clusters"
    execute(line)

    line = f"{bcftools} sort -Oz -o {vcfFolder}/{outputName}_Sort.vcf.gz {vcfFolder}/{outputName}.vcf.gz"
    execute(line)

    line = f"{bcftools} norm -Oz -c x -f {genomeReference} -o {vcfFolder}/{outputName}_Norm.vcf.gz {vcfFolder}/{outputName}_Sort.vcf.gz"
    execute(line)

    line = f"{bcftools} index {vcfFolder}/{outputName}_Norm.vcf.gz"
    execute(line)


    #line = f"{bcftools} plugin fill-tags -Oz -o  {vcfFolder}/{outputName}_Fill.vcf.gz {vcfFolder}/{outputName}_Norm.vcf.gz"
    #execute(line)
    line = f"{bcftools} index {vcfFolder}/{outputName}_Norm.vcf.gz"
    execute(line)

    return f"{vcfFolder}/{outputName}_Norm.vcf.gz"

def generateGTC(iaap, bpm, egt, folder, outFolder, threads):
    threadsTest = 4

    line = f"{iaap} gencall {bpm} {egt} -f {folder} --output-gtc {outFolder} --gender-estimate-call-rate-threshold -0.1 -t {threadsTest}"
    execute(line)

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


def basicVCFQC(VCF, parametersQC, bcftools, vcfFolder, vcfName):
    line = f'{bcftools}  view -i \"'

    first = True

    for parameter in parametersQC:
        if first:
            line = f'{line}INFO/{parameter}{parametersQC[parameter]}'
            first = False
        else:
            line = f'{line} && INFO/{parameter}{parametersQC[parameter]}'

    line = f'{line}\" -Oz -o {vcfFolder}/{vcfName}_QC.vcf.gz {VCF}'
    execute(line)

    return f"{vcfFolder}/{vcfName}_QC.tar.gz"


def basicGenotypingQC(VCF, folder, name, plink2):
    folderPGEN = f"{folder}/PGEN"
    execute(f"mkdir {folderPGEN}", True)

    #Convert
    execute(f"{plink2} --vcf {VCF} --make-pfile --out {folderPGEN}/{name}", True)

    #Remove missing
    execute(f"{plink2} --pfile {folderPGEN}/{name} --geno 0.05 --mind 0.05 --make-pgen --out {folderPGEN}/{name}_missing", True)

    #Remove big deviation on HWE
    execute(f"{plink2} --pfile {folderPGEN}/{name}_missing --hwe 1e-10 --make-pgen --out {folderPGEN}/{name}_HWE", True)

    #PGEN to VCF
    execute(f"{plink2} --pfile {folderPGEN}/{name}_HWE --recode vcf-iid --out {folderPGEN}/{name}_QC", True)
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
                    dictInfo[ID][header] = ""
                    dictInfo[ID][header] = split[i]
            dictInfo[ID]["check"] = True

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
    print(f"\tThe VCF has {count} varaints, which {countToSearch} is new \n")
    return searches


def dbsnpData(dictData):
    alleles = []
    allelesFreq = []

    rsID = dictData['rsid']
    if 'alleles' in dictData:
        for observation in dictData['alleles']:
            allele = observation['allele']
            if 'dbgap_popfreq' not in observation['freq']:
                freq = "NA"
            else:
                freq = observation['freq']['dbgap_popfreq']

            alleles.append(allele)
            allelesFreq.append(freq)
    else:
        alleles = ["NA", "NA"]
        allelesFreq = ["NA", "NA"]
    return alleles, allelesFreq, rsID


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
    for data in mv.getvariants(toSearch, fields=fields, returnall = True):
        if "notfound" in data:
            notFound.append(data['query'])
        else:
            if isinstance(data, dict):

                #DBSNP
                rs = ""
                alleles = ["NA", "NA"]
                freqs = ["NA", "NA"]

                #CADD
                annotype = ""
                consdetail = ""
                consequence = ""
                consscore = ""
                geneID = ""
                geneName = ""
                polyphen = ""
                sift = ""
                phred = ""

                #CLINVAR
                name = ""
                preferredName = ""
                omimList = []
                mondoList = []
                medgenList = []
                orphanetList = []
                humanPhenotypeOntologyList = []

                if 'cadd' in data:
                    annotype, consdetail, consequence, consscore, geneID, geneName, polyphen, sift, phred = caddData(data["cadd"])

                if 'dbsnp' in data:
                    alleles, freqs, rs = dbsnpData(data["dbsnp"])

                if 'clinvar' in data:
                    name, preferredName, omimList, mondoList, medgenList, orphanetList, humanPhenotypeOntologyList = getRCVData(data["clinvar"]["rcv"])

                omimOut = getDataToPrint(omimList)
                mondoOut = getDataToPrint(mondoList)
                medgenOut = getDataToPrint(medgenList)
                orphanetOut = getDataToPrint(orphanetList)
                humanPhenotypeOntologyOut = getDataToPrint(humanPhenotypeOntologyList)

                annotypeOut = getDataToPrint(annotype)
                consdetailOut = getDataToPrint(consdetail)
                consequenceOut = getDataToPrint(consequence)
                consscoreOut = getDataToPrint(consscore)

                keys = ['rs', 'Allele 1', 'Allele 2', 'Freq(A1)', 'Freq(A2)', 'CADD Phred', 'CADD Annotype',
                        'CADD Consequence Detail', 'CADD Consequence', 'CADD Consequence Score', 'CADD Gene ID',
                        'CADD Gene Name', 'Sift Cathegory', 'Polyphen Cathegory', 'CLINVAR Name', 'CLINVAR Preferred Name',
                        'CLINVAR Preferred Name', 'OMIM IDs', 'MONDO IDs', 'MEDGEN IDs', 'Orphanet IDs', 'Human Phenotype Ontology']

                query = data["query"]
                variantInfo[query] = {}
                for key in keys:
                    variantInfo[query][key] = ""

                variantInfo[query]['rs'] = rs
                variantInfo[query]['Allele 1'] = alleles[0]
                variantInfo[query]['Allele 2'] = alleles[1]

                if freqs[0] == "" or "NA":
                    freqs[0] = "-1"
                variantInfo[query]['Freq(A1)'] = float(freqs[0])

                if freqs[1] == "" or "NA":
                    freqs[1] = "-1"
                variantInfo[query]['Freq(A2)'] = float(freqs[1])
                if phred == "":
                    phred = "-1"
                variantInfo[query]['CADD Phred'] = float(phred)

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
              'dbsnp.alleles', 'dbsnp.rsid', 'cadd.gene.gene_id', 'cadd.gene.genename', 'cadd.annotype',
              'cadd.consdetail', 'cadd.consequence', 'cadd.consscore', 'cadd.sift.cat', 'cadd.polyphen.cat',
              'cadd.phred']

    toSearch = list(toSearchDict.keys())
    notFound, variantInfo = makeSearch(toSearch, variantInfo, fields)

    notFoundToSearch = {}
    for ID in notFound:
        A2 = toSearchDict[ID]["A1"]
        A1 = toSearchDict[ID]["A2"]
        prefix = toSearchDict[ID]["prefix"]

        notFoundToSearch[f"{prefix}{A1}>{A2}"] = ID

    toSearchNotFound = list(notFoundToSearch.keys())
    notFound, variantInfo = makeSearch(toSearchNotFound, variantInfo, fields)

    return variantInfo

def openVCFAndSearch(vcfName, variantInfo, folderAnot, name, correspondence):
    toSearch = readVCFs(vcfName, variantInfo)
    variantInfo = searchSNPs(variantInfo, toSearch)

    header = True
    fileReference = open(f"{folderAnot}/{name}_AnnotTable.tsv", "w")
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

    if correspondence != "":
        fileCorrespondence = open(correspondence)
        dictCorrespondence = {}
        for line in fileCorrespondence:
            split = line.strip().split()
            dictCorrespondence[split[0]] = split[1]

    print(f"Open VCF file again({vcfName}) ")

    if vcfName[-2:] == ".gz":
        file = gzip.open(vcfName)
    else:
        file = open(vcfName)

    header = True
    headerPrinted = False
    dictOut = {}
    for line in file:
        if vcfName[-2:] == ".gz":
            line = line.decode("utf-8")

        if header:
            if "#CHROM" in line:
                headerIDsTemp = line.strip().split()
                headerIDs = []

                for i in range(len(headerIDsTemp)):
                    if i < 9:
                        headerIDs.append(ID)
                    else:
                        ID = headerIDsTemp[i]
                        if ID not in dictOut:
                            dictOut[ID] = open(f"{folderAnot}/Annot_{ID}.tsv", 'w')
                            headerIDs.append(ID)
                        else:
                            count = 1
                            while f"{ID}_{count}" in headerIDs:
                                count = count + 1

                            headerIDs.append(f"{ID}_{count}")
                            dictOut[f"{ID}_{count}"] = open(f"{folderAnot}/Annot_{ID}_{count}.tsv", 'w')

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

                if variantInfo[query]['Freq(A1)'] < variantInfo[query]['Freq(A2)']:
                    maf = variantInfo[query]['Freq(A1)']
                else:
                    maf = variantInfo[query]['Freq(A2)']

                if variantInfo[query]['CADD Phred'] > 15 and maf < 0.05:
                    for i in range(9, len(split)):
                        ID = headerIDs[i]
                        if correspondence != "":
                            ID = dictCorrespondence[ID]

                        information = split[i].split(":")
                        if "1" in information[0]:
                            information[0].replace("0", A0).replace("1", A1)
                            dictOut.write(f"{query}\t{information[0]}")
                            for item in variantInfo[query]:
                                dictOut.write(f"\t{variantInfo[query][item]}")
                            dictOut.write(f"\n")
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



    folder = createFolder(args.outputFolder)
    folderGTC = createFolder(folder+"/GTC/")
    #generateGTC(args.iaap, args.bpm, args.egt, args.folder, folderGTC, args.threads)
    #generatePED(args.iaap, args.bpm, args.egt, args.folder, folderGTC, args.threads)

    allGTCFolder = createFolder(folder+"/GTCAll/")
    vcfFolder = createFolder(folder+"/VCFs/")
    VCF = gtc2VCF(args.bcftools, args.bpm, args.egt, args.csv, folderGTC, allGTCFolder, args.genomeReference, vcfFolder, args.threads, args.outputName)
    parametersQC = readFileQC(args.qc)

    VCFQC = basicVCFQC(VCF, parametersQC, args.bcftools, vcfFolder, args.outputName)

    VCFFinal = basicGenotypingQC(VCFQC, folder, args.outputName, args.plink2)

    anotFolder = createFolder(folder + "/Anot/")
    variantInfo = {}
    if args.previousSearch != "":
        variantInfo = readPreviousSearchFile(args.previousSearch)

    variantInfo = openVCFAndSearch(VCFFinal, variantInfo, anotFolder, args.outputName, args.correspondenceList)


