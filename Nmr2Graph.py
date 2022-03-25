import pybel
import pandas as pd
from pybel.dsl import Protein
from pybel.dsl import Abundance
from pybel.dsl import Pathology
from pybel.dsl import BiologicalProcess
from pybel.dsl import Population
from pybel.dsl import Gene
from pybel.dsl import MicroRna
from pybel.dsl import Fragment
import chembl_webresource_client
import openpyxl
import networkx as nx
from pybel.io.jupyter import to_jupyter
import matplotlib.pyplot as plt
import chembl_webresource_client
from chembl_webresource_client.new_client import new_client
import pubchempy
import pickle
import re

#function to retrieve mechanisms from ChEMBL

#nmrGraph = pybel.BELGraph(name='NMRgraph')
nmrData = pd.read_csv('NMR_data_new.csv')
nmr2chembl = pd.read_excel('Protein-Screening-Summary_CHEMBL30_Actives_teams_max08tan.xlsx')
#itmp_chem = pd.read_csv('C:\\Users\\reagon.karki\\Documents\\ITMP\\EOSC Future\\ITMP ChEMBL data\\CHEMBL4495564.csv',sep=';',usecols=['ChEMBL ID'])
nmr2pchem = pd.read_excel('PubChem_Actives_Literature_NCATS.xlsx',sheet_name='Similar_065')


def RetMech(chemblIds):
    getMech = new_client.mechanism
    mechList = []
    for i in range(len(chemblIds)):
        mechs = getMech.filter(molecule_chembl_id=chemblIds[i]).only(['mechanism_of_action','target_chembl_id'])
        #mechs = getMech.filter(molecule_chembl_id=chemblIds[i])
        print(mechs)
        mechList.append(list(mechs))
    named_mechList = dict(zip(chemblIds,mechList))
    named_mechList = {k: v for k, v in named_mechList.items() if v}
    return(named_mechList)

def RetDrugInd(chemblIDs):
    getDrugInd = new_client.drug_indication
    drugIndList = []
    for i in range(len(chemblIDs)):
        drugInd = getDrugInd.filter(molecule_chembl_id=chemblIDs[i]).only('mesh_heading')
        #drugInd = getDrugInd.filter(molecule_chembl_id=chemblIDs[i])
        print(drugInd)
        drugIndList.append(list(drugInd))
    named_drugIndList = dict(zip(chemblIDs,drugIndList))
    named_drugIndList = {k: v for k, v in named_drugIndList.items() if v}
    return(named_drugIndList)

def RetAct(chemblIds,j):
    GetAct = new_client.activity
    ActList = []
    #for i in range(len(chemblIds)):
    #print(chemblIds[0])
    for i in range(len(itmp_chem)):
        filtered_list=['assay_chembl_id','assay_type','pchembl_value','target_chembl_id','target_organism','bao_label']
        acts = GetAct.filter(molecule_chembl_id=chemblIds[i],pchembl_value__isnull=False,assay_type_iregex='(B|F)',target_organism='Homo sapiens').only(filtered_list)
        #acts = GetAct.filter(molecule_chembl_id=chemblIds[i],pchembl_value__isnull=False,target_organism='Homo sapiens').only(filtered_list)
        #acts = GetAct.filter(molecule_chembl_id=chemblIds[i],pchembl_value__isnull=False)
        j=j+1

        acts = [d for d in acts if d.get('target_organism') == 'Homo sapiens']
        acts = [d for d in acts if d.get('bao_label') == 'single protein format']
        acts = [d for d in acts if d.get('type') in ['Ki', 'IC50']]
        acts = [d for d in acts if float(d.get('pchembl_value')) >= 6]
        print(j)
        #print(len(acts))
        acts = acts[:5]
        print(acts)
        ActList.append(list(acts))
    #print(ActList)
    named_ActList = dict(zip(chemblIds,ActList))
    named_ActList = {k: v for k, v in named_ActList.items() if v}
    return(named_ActList)

def filter_graph(mainGraph, vprotList):
    nsp_list = []
    for u, v, data in mainGraph.edges(data=True):
        if u.name in vprotList or v.name in vprotList:
            nsp_list.append(u)
            nsp_list.append(v)

    for u, v, data in mainGraph.edges(data=True):
        if 'Similarity' not in data:
            continue
        if data['Similarity'] >= .65 and (u in nsp_list or v in nsp_list):
            # yield u, v, k, data
            nsp_list.append(u)
            nsp_list.append(v)

    # print(nsp_list)

    nsp_graph = mainGraph.subgraph(nsp_list)
    # return(to_jupyter(nsp_graph))
    final_nsp = [node for node in nsp_graph.nodes() if isinstance(node, pybel.dsl.Protein) and node.name in vprotList]
    # print(final_nsp)
    for u, v, data in nsp_graph.edges(data=True):
        if u.namespace in ['Enamine', 'ChEMBL', 'CID'] and v.namespace in ['Enamine', 'ChEMBL', 'CID']:
            final_nsp.append(u)
            final_nsp.append(v)

    # return(nsp_graph.subgraph(final_nsp))
    nsp_graph = nsp_graph.subgraph(final_nsp)
    # chembl_act = []
    nsp_nodes = [node for node in nsp_graph.nodes()]
    for u in nsp_graph.nodes():
        if u.namespace == 'ChEMBL':
            # test = nmrGraph.neighbors(u)
            chem = [n for n in mainGraph.neighbors(u)]
            nsp_nodes.extend(chem)

    return (mainGraph.subgraph(nsp_nodes))

#MechList to graph
for i in mechList:
    #print(i)
    #break
    #print(named_mechList[i])
    #print(len(named_mechList[i]))
    #break
    for j in range(len(mechList[i])):
        #print(named_mechList[i][j]['mechanism_of_action'])
        #print(named_mechList[i][j]['target_chembl_id'])
        #print(i)
        #break
        nmrGraph.add_association(Pathology(namespace='ChEMBL',name=i),BiologicalProcess(namespace='MOA',name=mechList[i][j]['mechanism_of_action']),citation='ChEMBL database',evidence='ChEMBL query')
        if not mechList[i][j]['target_chembl_id'] == None:
            nmrGraph.add_association(Pathology(namespace='ChEMBL',name=i),MicroRna(namespace='HP',name=mechList[i][j]['target_chembl_id']),citation='ChEMBL database',evidence='ChEMBL query')


#DrugIndList to graph

for i in drugIndList:
    # print(i)
    # break
    # print(named_drugIndList[i])
    # print(len(named_drugIndList[i]))
    # break
    for j in range(len(drugIndList[i])):
        # print(named_drugIndList[i][j]['mesh_heading'])
        # print(i)
        # break
        nmrGraph.add_association(Pathology(namespace='ChEMBL', name=i),
                                 Population(namespace='Disease', name=drugIndList[i][j]['mesh_heading']),
                                 citation='ChEMBL database', evidence='ChEMBL query')


#ActList to graph
#Only for NMR
for i in ActList:
    for j in range(len(ActList[i])):
        #print(named_mechList[i][j]['mechanism_of_action'])
        #print(named_mechList[i][j]['target_chembl_id'])
        #print(i)
        #break
        #nmrGraph.add_association(Pathology(namespace='ChEMBL',name=i),BiologicalProcess(namespace='MOA',name=mechList[i][j]['mechanism_of_action']),citation='ChEMBL database',evidence='ChEMBL query')
        if not ActList[i][j]['target_chembl_id'] == None:
            nmrGraph.add_association(Pathology(namespace='ChEMBL',name=i),MicroRna(namespace='HP',name=ActList[i][j]['target_pref_name']),citation='ChEMBL database',evidence='ChEMBL query')


#a = ['CHEMBL455','CHEMBL960']
#test= RetMech(a)
#print(test)

#itmp_moa = RetMech(list(itmp_chem['ChEMBL ID']),0)
print(len(itmp_moa))

#itmp_chem_moa2pickle = open('itmp_chem_MechList2Pickle','wb')
pickle.dump(itmp_moa,itmp_chem_moa2pickle)
itmp_chem_moa2pickle.close()

# infile = open('itmp_chem_moa2pickle','rb')
# itmp_chem_moa = pickle.load(infile)
# infile.close()

def RetDrugInd(chemblIDs):
    getDrugInd = new_client.drug_indication
    drugIndList = []
    for i in range(len(chemblIDs)):
        drugInd = getDrugInd.filter(molecule_chembl_id=chemblIDs[i]).only('mesh_heading')
        print(drugInd)
        drugIndList.append(list(drugInd))
    return(drugIndList)