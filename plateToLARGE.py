file = open("PlateToIlluminaCode.txt")

dictPlateIllumina = {}

for line in file:
    split = line.strip().split(",")
    dictPlateIllumina[split[0]] = split[1]

file.close()
file = open("IlluminaCodeToLARGE.txt")
fileOut = open("Correspondence.txt", 'w')

for line in file:
    split = line.strip().split("\t")
    fileOut.write(f"{dictPlateIllumina[split[0]]}\t{split[1]}\n")

fileOut.close()