import sys
import pandas as pd

snapGenePath = sys.argv[1]
gffPath = sys.argv[2]

def snapToGff(snapGenePath, gffPath):
    snapGene = pd.read_table(snapGenePath, delimiter="\t", names=["name", "coords", "length", "direction", "type"])
    snapGene['coords'] = snapGene['coords'].str.replace(',', '')
    snapGene[['start', 'end']] = snapGene['coords'].str.split("\..", expand=True)
    snapGene.drop(columns=['coords'], inplace=True)
    snapGene['seqid'] = snapGenePath.split('Features from ')[1].split('.txt')[0]
    snapGene['source'] = 'snapGene'
    snapGene['strand'] = snapGene['direction'].apply(lambda x: "+" if x == "=>" else "-" if x == "<=" else ".")
    snapGene['type'] = snapGene['type'].str.strip()
    snapGene['phase'] = snapGene["type"].apply(lambda x: "0" if x == "CDS" else ".")
    snapGene.drop(columns=['direction'], inplace=True)
    snapGene['score'] = "."
    snapGene['name'] = snapGene['name'].str.strip()
    snapGene['name'] = snapGene['name'].str.replace(' ', '_')
    snapGene['name'] = snapGene['name'].apply(lambda x: "ID=" + x)

    #Go through and check if any entries in 'name' are duplicated, if so add a number (which iterates up for each example of the name) to the end of the name
    nameCounts = snapGene['name'].value_counts()
    for name, count in nameCounts.items():
        if count > 1:
            nameIndices = snapGene[snapGene['name'] == name].index
            for i, index in enumerate(nameIndices):
                snapGene.at[index, 'name'] = snapGene.at[index, 'name'] + "_" + str(i + 1)

    snapGene = snapGene[['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'name']]
    snapGene.to_csv(gffPath, sep="\t", header=False, index=False)
    return(snapGene)

snapToGff(snapGenePath, gffPath)
