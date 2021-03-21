import csv
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", help="Input File")

args = parser.parse_args()

print("Input {}".format(args.input))
name = args.input
#name = "WD 2TB/FiguresPapers/FabiansPaper_smiFISH/Fabian_smFISH/2019-05-16/Outlines_DAPI_Max/467_TPI1_FBA1_PGK1_GPM1/_FQ_outline/467_TPI1-488_FBA1-546_PGK1-590_GPM1-647_Series010_1_cmle_ch00_outline.txt"
outDet = name[:-4] + "_Cell_OUTLINES_LONG.csv"
nucoutDet = name[:-4] + "_Nuc_OUTLINES_LONG.csv"
cellStarts = []
cellEnds = []
nucStarts = []
nucEnds = []
xpos = []
ypos = []
print("Generating outlines from" + " "+name)
#FISHQuants =
with open(name) as tsv:
    reader = csv.reader(tsv, delimiter="\t")
    FISH = list(reader)
#Get cell starts, spot starts and spot ends. These can be used to
# 1) Attribute spots to specific cells
# 2) Determine the number of spots for each cell
# 3) Determine the total number of spots in the image and therefore the average spots per cell
i = 0
while i < len(FISH):
    if "CELL_START" in FISH[i]:
        cellStarts.append(i)
    elif "Nucleus_START" in FISH[i]:
        nucStarts.append(i)
    elif "Nucleus_END" in FISH[i]:
        nucEnds.append(i)
    elif "CELL_END" in FISH[i]:
        cellEnds.append(i)
    elif "X_POS" in FISH[i]:
        xpos.append(i)
    elif "Y_POS" in FISH[i]:
        ypos.append(i)
    i += 1
i = 0
#if len(cellStarts) != len(spotStarts):
#    print("Number of Cells does not equal number of spots cell identifiers, check the original .txt for cells without any spots")
#    exit(code = 1)
#else:
print("Number of Cells =" + " " + str(len(cellStarts)))
cellHeader = "Cell"
countHeader = "Spots"
xHeader = "X"
yHeader = "Y"
detHeader = ['X_POS', 'Y_POS']
countRows = [[cellHeader] + [countHeader]]
detRows = [[cellHeader] + detHeader[0:]]
while i < len(cellStarts):
    cellNum = i+1
    x = xpos[i]
    y = ypos[i]
    xp = FISH[x][1:]
    yp = FISH[y][1:]
    z = 0
    while z < len(xp)-1:
        current = [str(cellNum), str(xp[z]), str(yp[z])]
        detRows.append(current)
        z = z + 1
    i = i + 1
with open(outDet, "w") as d:
    writer = csv.writer(d)
    writer.writerows(detRows)
    d.close()
print("Cell outlines Completed")
detRows = [[cellHeader] + detHeader[0:]]
i = 0
while i < len(nucStarts):
    cellNum = i+1
    x = xpos[i]
    y = ypos[i]
    xp = FISH[x][1:]
    yp = FISH[y][1:]
    z = 0
    while z < len(xp)-1:
        current = [str(cellNum), str(xp[z]), str(yp[z])]
        detRows.append(current)
        z = z + 1
    i = i + 1
with open(nucoutDet, "w") as d:
    writer = csv.writer(d)
    writer.writerows(detRows)
    d.close()
print("Nucleus outlines Completed")
