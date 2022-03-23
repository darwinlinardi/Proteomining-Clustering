import os
import sys
import pandas as pd

if len(sys.argv) == 1:
    print("""The script identifies BGC by the location of the ORF in the genome.
GeneMark can differentiate the nodes where the proteins are from. As such this is suitable for partial genomes.
Usage: python genome-script.py input proteindb genegap min-members
Arguments:
Input = excel or csv file with columns ORF No., log2(fold-change), -log10(p-value)
proteindb = fasta or faa files.
Optional:
genegap: Maximum allowed gaps between ORFs, Default = 5
min-members: Minimum nodes in a cluster to be considered, Default = 8
""")
    exit()

genegap = 5
minmember = 8

print("Save file as:")
out_name = input('> ')
dirName = "./" + out_name + "-network"
# Create target Directory if don't exist
if not os.path.exists(dirName):
    os.mkdir(dirName)
else:
    pass

input_file = sys.argv[1]
genemarkfile = sys.argv[2]
if len(sys.argv) == 4:
    genegap = int(sys.argv[3])
elif len(sys.argv) ==  5:
    genegap = int(sys.argv[3])
    minmember = int(sys.argv[4])
genemark_temp = open(genemarkfile, "r")
genemark_file = genemark_temp.readlines()
genemark_temp.close()
genemark = []
for temp_line in genemark_file:
    if ">" in temp_line:
        genemark.append(temp_line)
#Contains the start and end location of genes(key) and in which nodes
nodedict = {}
genestartdict = {}
geneenddict = {}
for i in genemark:
    print(i)
    if "DECOY" in i or "cont" in i:
        continue
    elif 'gene_' in i:
        genenum = int(i.split('gene_')[1].split('|')[0])
        if "NODE_" in i:
            nodedict[genenum] = int(i.split('NODE_')[1].split('_')[0])
            genestartdict[genenum] = int(i.split('|')[4])
            geneenddict[genenum] = int(i.split('|')[5].split('\t')[0])
        else:
            nodedict[genenum] = 1
            genestartdict[genenum] = 0
            geneenddict[genenum] = 0

    else:
        if "." in i:
            genenum = int(i.split('.')[-1])
        elif "_" in i:
            genenum = int(i.split('.')[-1])
        else:
            print("Use '.' or '_' as delimiter before the ORF No.")
        nodedict[genenum] = 1
        genestartdict[genenum] = 0
        geneenddict[genenum] = 0



networkfile = pd.read_csv(input_file, header=0, names=['prot_name', 'log2(FC)', 'NEG(log10(pval))'], index_col=None)
networkfile = networkfile.sort_values(networkfile.columns[0], ignore_index=True)
temp_category_list = []
temp_nodenum = []
temp_genestart = []
temp_geneend = []
#Assigns scores to the proteins
for i in range(len(networkfile)):
    prot_num = networkfile.loc[i,"prot_name"]
    temp_nodenum.append(nodedict[prot_num])
    temp_genestart.append(genestartdict[prot_num])
    temp_geneend.append(geneenddict[prot_num])
    if networkfile.loc[i,"log2(FC)"] >= 0.584963 and networkfile.loc[i,"NEG(log10(pval))"] >= 1.30103 :
        if networkfile.loc[i, "log2(FC)"] >= 1 and networkfile.loc[i, "NEG(log10(pval))"] >= 2:
            temp_category_list.append(10)
        else:
            temp_category_list.append(5)
    elif networkfile.loc[i,"log2(FC)"] <= -0.584963 and networkfile.loc[i,"NEG(log10(pval))"] >= 1.30103:
        if networkfile.loc[i, "log2(FC)"] <= -1 and networkfile.loc[i, "NEG(log10(pval))"] >= 2:
            temp_category_list.append(-10)
        else:
            temp_category_list.append(-5)
    else:
        temp_category_list.append(0)
networkfile['Category'] = temp_category_list
networkfile['Node'] = temp_nodenum
networkfile['Gene Start'] = temp_genestart
networkfile['Gene End'] = temp_geneend
networkfile = networkfile[networkfile.Category != 0]
networkfile = networkfile.reset_index(drop=True)
#This generates the edges of the nodes
edgelist = []
flaggednodes = []
for i in range(len(networkfile)):
    j = i + 1
    while j < len(networkfile):
        dist_temp = networkfile.loc[j,"prot_name"] - networkfile.loc[i,"prot_name"]
        nodei = networkfile.loc[i,"Node"]
        nodej = networkfile.loc[j,"Node"]
        #Checks if the two genes are within the threshold limit and lies within the same node
        if (nodei == nodej) and (dist_temp <= genegap):
            temp_edge = "%s,%s" %(networkfile.loc[i,"prot_name"], networkfile.loc[j,"prot_name"])
            if dist_temp == genegap:
                flaggednodes.append(networkfile.loc[i,"prot_name"])
            edgelist.append(temp_edge)
            if networkfile.loc[j,"prot_name"] == 4914:
                print(networkfile.loc[i,"prot_name"])
                print(nodei)
                print(nodej)
        else:
            break
        j+= 1
#Determines number of clusters
clustdict = {} # A dictionary with cluster number as key to access the list of nodes that are in a cluster
fulledgedict = {} # A dictionary with cluster number as key to access the edges that are associated with the cluster
for i in range(len(edgelist)):
    entry1 = edgelist[i].split(',')[0]
    entry2 = edgelist[i].split(',')[1]
    j = 1
    k = len(clustdict.keys()) + 1
    if i == 0:
        clust1 = []
        clust1.append(entry1)
        clust1.append(entry2)
        clustdict[1] = clust1
        fulledgedict[1] = '{"source": "%s", "target": "%s", "value": 10},' % (entry1,entry2)

    else:
        #Scans through dictionary to check if the new edges are connected to them.
        while j < k:
            temp_clust = clustdict[j]
            #if any of the proteins are in an existing network, it will append to the network and deduplicate it
            if entry1 in temp_clust or entry2 in temp_clust:
                temp_clust.append(entry2)
                temp_clust.append(entry1)
                temp_clust2 = list(set(temp_clust))
                clustdict[j] = temp_clust2
                temp_edge = fulledgedict[j]
                temp_edge += '{"source": "%s", "target": "%s", "value": 10},' % (entry1,entry2)
                fulledgedict[j] = temp_edge
                j = 1
                break
            else:
                j += 1
        else:
            clustnew = []
            clustnew.append(entry1)
            clustnew.append(entry2)
            clustdict[k] = clustnew
            fulledgedict[k] = '{"source": "%s", "target": "%s", "value": 10},' % (entry1,entry2)

clusternumber = []
for i in range(len(networkfile)):
    prot_name = networkfile.loc[i,"prot_name"]
    for j in clustdict.keys():
        if len(clustdict[j]) < minmember:
            continue
        else:
            for k in clustdict[j]:
                if prot_name == int(k):
                    clusternumber.append(j)
                    break
                else:
                    continue
    if i+1 != len(clusternumber):
        clusternumber.append(0)

networkfile["Cluster"] = clusternumber
#function that generates the html of the network files
def generatehtml(nodelist, edgedict, outname, size):
    #nodelist is a dataframe containing the protname and score that indicate the color
    #edgelist is a list of clusters that is to be included in the html
    #edgedict is a dictionary with cluster number as key to access the edges that are associated with the cluster
    #size is the size of the html file used for the generation of the cluster graph
    outputjson = """{"nodes":["""
    temp_clusterlist = []
    for i in list(nodelist.index.values):
        protid = nodelist.loc[i,"prot_name"]
        protscore = nodelist.loc[i,"Category"]
        temp_clusterlist.append(nodelist.loc[i,"Cluster"])
        if protscore == 10:
            protgroup = 1
            groupcolor = "#00ff00"
        elif protscore == 5:
            protgroup = 2
            groupcolor = "#007200"
        elif protscore == -5:
            protgroup = 3
            groupcolor = "#720000"
        elif protscore == -10:
            protgroup = 4
            groupcolor = "#ff0000"
        else:
            print("Something is amiss!")
        outputjson += '{"id": "%s", "group": %d, "color": "%s"},' % (protid, protgroup, groupcolor)
    outputjson = outputjson[:-1]
    outputjson += """],
    "links": ["""
    temp_clusterlist = set(temp_clusterlist)
    for i in temp_clusterlist:
        if i == 0:
            continue
        else:
            outputjson += edgedict[i]

    outputjson = outputjson[:-1]
    outputjson += """]
    }
    """
    temp_file = open(dirName + "/" + outname + "-orfnetwork.js", "w")
    temp_file.write(outputjson)
    temp_file.close()

    if size == "Large":
        temp_loc = os.path.dirname(os.path.realpath(__file__))
        temp_file = open("""%s\cluster-network.html""" % temp_loc , "r")
        htmlfile = temp_file.readlines()
        temp_file.close()
        for i in range(len(htmlfile)):
            if "bananaslammajamma" in htmlfile[i]:
                htmlfile[i] = '''d3.json("%s", function(error, graph) {''' % (outname + "-orfnetwork.js")
                break
        htmlout = ('').join(htmlfile)
        temp_file = open(dirName + "/" + outname + ".html", "w")
        temp_file.write(htmlout)
        temp_file.close()
    elif size == "Small":
        temp_loc = os.path.dirname(os.path.realpath(__file__))
        temp_file = open("""%s\cluster-network-small.html""" % temp_loc , "r")
        htmlfile = temp_file.readlines()
        temp_file.close()
        for i in range(len(htmlfile)):
            if "bananaslammajamma" in htmlfile[i]:
                htmlfile[i] = '''d3.json("%s", function(error, graph) {''' % (outname + "-orfnetwork.js")
        htmlout = ('').join(htmlfile)
        temp_file = open(dirName + "/" + outname + ".html", "w")
        temp_file.write(htmlout)
        temp_file.close()


#Generate output with all of the data
generatehtml(networkfile, fulledgedict, out_name+"-raw-network" , "Large")

temp_rawscore = []
temp_score = []
for i in list(networkfile.index.values):
    if networkfile.loc[i,"Cluster"] == 0:
        temp_rawscore.append(0)
        temp_score.append(0)
    else:
        rawscore = networkfile[networkfile.Cluster == networkfile.loc[i,"Cluster"]]['Category'].sum()
        temp_col = networkfile[networkfile.Cluster == networkfile.loc[i,"Cluster"]]["prot_name"]
        temp_range = temp_col.max() - temp_col.min()
        finalscore = abs(rawscore/temp_range)
        temp_rawscore.append(rawscore)
        temp_score.append(finalscore)
networkfile["Raw_Score"] = temp_rawscore
networkfile["Score"] = temp_score

clean_network = networkfile[networkfile.Cluster != 0]
clean_network = networkfile[networkfile.Score.between(1.5,3)]
clean_network_highscore = networkfile[networkfile.Score >= 3.0]
finalclusterlist = clean_network["Cluster"].tolist()
finalclusterlist = set(finalclusterlist)
generatehtml(clean_network, fulledgedict, out_name+"-fringe-network" , "Large")
generatehtml(clean_network_highscore, fulledgedict, out_name+"-hiscore-network" , "Large")
for i in finalclusterlist:
    temp_network = clean_network[clean_network.Cluster == i]
    temp_name = out_name+"-cluster"+str(i)
    generatehtml(temp_network, fulledgedict, temp_name , "Small")
temp_outname = out_name + "-output.xlsx"
networkfile.to_excel(temp_outname)
