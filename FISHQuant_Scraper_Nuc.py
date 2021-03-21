import csv
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", help = "Input File")

args=parser.parse_args()

print("Input {}".format(args.input))
name = args.input
out = name[:-4] +"_SCRAPED.csv"
outDet = name[:-4] +"_SCRAPED_DETAILED.csv"
cellStarts = []
spotStarts =[]
spotEnds = []
print("Scraping data from" + " "+name)
#FISHQuants =
with open(name) as tsv:
    reader =  csv.reader(tsv, delimiter = "\t")
    FISH = list(reader)
#Get cell starts, spot starts and spot ends. These can be used to
# 1) Attribute spots to specific cells
# 2) Determine the number of spots for each cell
# 3) Determine the total number of spots in the image and therefore the average spots per cell
i=0
while i <len(FISH):
    if "CELL_START" in FISH[i]:
        if "SPOTS_START" in FISH[i+10]:
            cellStarts.append(i)
    elif "SPOTS_START" in FISH[i]:
        spotStarts.append(i)
    elif "SPOTS_END" in FISH[i]:
        spotEnds.append(i)
    i += 1
i=0
if len(cellStarts) != len(spotStarts):
    print("Number of Cells does not equal number of spots cell identifiers, check the original .txt for cells without any spots")
    exit(code = 1)
else:
    print("Number of Cells =" + " " + str(len(cellStarts)))
    cellHeader = "Cell"
    countHeader = "Spots"
    detHeader = FISH[spotStarts[0]+1]
    countRows = [[cellHeader] + [countHeader]]
    detRows = [[cellHeader] + detHeader[0:]]
    while i < len(cellStarts):
        cellNum = i+1
        spotCount = (spotEnds[i]-1)-(spotStarts[i]+2)
        detStart = spotStarts[i]+2
        detEnd = spotEnds[i]-1
        x = detStart
        while x < detEnd:
            detailed = FISH[x]
            spotsDetailed = [str(cellNum)] + detailed
            detRows.append(spotsDetailed)
            x = x+1
        countRows.append([str(cellNum)] + [str(spotCount)])
        i = i+1
    with open(out, "w") as f:
        writer = csv.writer(f)
        writer.writerows(countRows)
        f.close()
    with open(outDet, "w") as d:
        writer = csv.writer(d)
        writer.writerows(detRows)
        d.close()
    print("Scraping Completed")
