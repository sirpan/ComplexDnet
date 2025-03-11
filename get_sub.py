import os
from xml.etree import ElementTree as ET
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from rdkit.Chem import DataStructs
from rdkit.Chem import Draw


def makedir(dir_name):
    if os.path.exists('./' + dir_name) is not True:
        os.makedirs('./' + dir_name)


def set_default(desc_path='/home/fpan/PaDEL/descriptors.xml'):
    tree = ET.parse(desc_path)
    root = tree.getroot()
    for desc in root[2]:
        desc.set('value', 'false')
    root[2][5].set('value', 'true')
    tree.write(desc_path, encoding='utf-8', xml_declaration=True)
    print('descriptors.xml have been reseted \n\n')


def select_finger(num, desc_path='/home/fpan/PaDEL/descriptors.xml'):
    tree = ET.parse(desc_path)
    root = tree.getroot()
    root[2][5].set('value', 'false')
    root[2][num].set('value', 'true')
    tree.write(desc_path, encoding='utf-8', xml_declaration=True)
    print('descriptors.xml setted to all calc  \n\n')


def df_reciever(finger, num):
    padel_path = '/home/fpan/PaDEL/PaDEL-Descriptor.jar'
    desc_path = '/home/fpan/PaDEL/descriptors.xml'
    input_file = 'huang.smi'
    # input_file = '/home/zhyu/Endocrine_Disruption/Origin_data/1.smi'
    output_file = './fingerprint/' + finger + '.out'
    # output_file = '/home/zhyu/Endocrine_Disruption/Origin_data/1.out'

    cmd = f'java -Djava.awt.headless=true -jar {padel_path} -maxruntime -1 -waitingjobs -1 -descriptortypes {desc_path} -dir {input_file} -file {output_file} -fingerprints -detectaromaticity -removesalt -standardizenitro -retainorder -log'

    select_finger(num)
    os.system(cmd)
    set_default()
    os.remove('./fingerprint/' + finger + '.out.log')
    # os.remove('/home/zhyu/Endocrine_Disruption/Origin_data/1.out.log')


def handle_finger(finger):
    infile = open('./fingerprint/' + finger + '.out', 'r', encoding='utf-8')
    # infile = open('/home/zhyu/Endocrine_Disruption/Origin_data/1.out', 'r', encoding='utf-8')
    outfile = open('./fingerprint/' + finger + '_True.txt', 'w', encoding='utf-8')
    # outfile = open('/home/zhyu/Endocrine_Disruption/Origin_data/1_True.txt', 'w', encoding='utf-8')

    dir_FP = {}
    dir_num = {}
    for index, content in enumerate(infile):
        content = content.strip().split(',')
        if index == 0:
            for num, FP in enumerate(content):
                dir_num[num] = FP
        else:
            for num, FP in enumerate(content):
                title = dir_num.get(num)
                if FP == '1':
                    if content[0].strip('"') not in dir_FP:
                        dir_FP[content[0].strip('"')] = []
                    dir_FP[content[0].strip('"')].append(title)
    for key, values in dir_FP.items():
        for sub in values:
            outfile.write('DRUG' + '\t' + key + '\t' + 'SUB' + '\t' + sub + '\t' + '1' + '\n')
    # os.remove('/home/zhyu/Endocrine_Disruption/Endocrine_disruption/fingerprint/' + finger + '.out')
    outfile.close()

def SMILES_to_Morgan_new(filename, radius, useFeatures=False, ofile='', nBits=1024):
    infile = open(filename, 'r')
    if not useFeatures:
        outfile = open(ofile + "ECFP" + str(radius * 2) + "_True.out", 'w')
        outfile1 = open(ofile + "ECFP" + str(radius * 2) + "_True.txt", 'w')
    else:
        outfile = open(ofile + "FCFP" + str(radius * 2) + "_True.out", 'w')
        outfile1 = open(ofile + "FCFP" + str(radius * 2) + "_True.txt", 'w')
    
    title = []
    for x in range(1, nBits + 1):
        title.append('Morgan' + str(x))
    outfile.write('Name' + ',' + ','.join(title) + '\n')
    
    for id, line in enumerate(infile):
        smiles = line.rstrip('\n')
        try:
            mol = Chem.MolFromSmiles(smiles)
            fp1_morgan = AllChem.GetMorganFingerprint(mol, radius, useFeatures=useFeatures)
            fp1_morgan_hashed = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits, useFeatures=useFeatures)
            
            list_bits = fp1_morgan_hashed.ToBitString()
            bit_list = [list_bits[i] for i in range(len(list_bits))]
            
            outfile.write('"AUTOGEN_Drug_' + str(id + 1) + '"' + ',' + ','.join(bit_list) + '\n')
            
            for sub in fp1_morgan.GetNonzeroElements().keys():
                outfile1.write('DRUG' + '\t' + 'AUTOGEN_Drug_' + str(id + 1) + '\t' + 'SUB' + '\t' + str(sub) + '\t' + '1' + '\n')
        except Exception as e:
            print(f"Error processing {smiles} in file {filename}: {e}")
    
    infile.close()
    outfile.close()
    outfile1.close()

def SMILES_to_Morgan(filename, radius, useFeatures = False, ofile = '', nBits = 1024):
    #
    infile = open(filename, 'r')
    if not useFeatures:
        outfile = open(ofile + "ECFP" + str(radius * 2) + "_True.out", 'w')
        outfile1 = open(ofile + "ECFP" + str(radius * 2) + "_True.txt", 'w')
    else:
        outfile = open(ofile + "FCFP" + str(radius * 2) + "_True.out", 'w')
        outfile1 = open(ofile + "FCFP" + str(radius * 2) + "_True.txt", 'w')
    # outfile = open("Morgan_False_" + str(radius) + "_" + str(nBits) + ".txt", 'w')
    title = []
    for x in range(1, nBits + 1):
        title.append('Morgan' + str(x))
    outfile.write('Name' + ',' + ','.join(title) + '\n')
    for id, line in enumerate(infile):
        smiles = line.rstrip('\n')
        try:
            mol = Chem.MolFromSmiles(smiles)
            fp1_morgan = AllChem.GetMorganFingerprint(mol, radius, useFeatures=useFeatures)
            fp1_morgan_hashed = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits, useFeatures=useFeatures)
        except:
            print(filename, smiles)
        list = []
        for i in range(0, len(fp1_morgan_hashed.ToBitString())):
            list.append(fp1_morgan_hashed.ToBitString()[i])
        outfile.write('"AUTOGEN_Drug_' + str(id + 1) + '"' + ',' + ','.join(list) + '\n')
        for sub in fp1_morgan.GetNonzeroElements().keys():
            outfile1.write('DRUG' + '\t' + 'AUTOGEN_Drug_' + str(id + 1) + '\t' + 'SUB' + '\t' + str(sub) + '\t' + '1' + '\n')
def GetRdkitDescriptors(smile_list,csv_file:str):
    mol_list=[]
    for smi in smile_list:
        mol_list.append(Chem.MolFromSmiles(smi))
    mols = [mol for mol in mol_list]
    descs = [desc_name[0] for desc_name in Descriptors._descList]
    desc_calc = MoleculeDescriptors.MolecularDescriptorCalculator(descs)
    descriptors = pd.DataFrame([desc_calc.CalcDescriptors(mol) for mol in mols])
    descriptors.columns = descs
    descriptors.index = smile_list
    index_list = list(map(str,list(range(len(mols)))))
    y = pd.DataFrame(index_list)
    y.index = smile_list
    y.columns = ["index"]
    dataset = pd.concat([y, descriptors], axis=1)
    dataset.to_csv(csv_file)


if __name__ == "__main__":
    makedir('./fingerprint')
    # for i in range(2,3):
    #     SMILES_to_Morgan('all_Can_SMILES.smi', radius=i, useFeatures=False, ofile = './bSDTNBI_FP/fingerprint/')
    #     SMILES_to_Morgan('all_Can_SMILES.smi', radius=i, useFeatures=True, ofile = './bSDTNBI_FP/fingerprint/')
    for num, finger in enumerate(
            ['CDK', 'CDKExt', 'EState', 'GraphFP', 'MACCS', 'PubChemFP', 'SubFP', 'KRFP',
             'AP2D','Morgan' ]):
        if finger not in ['Morgan']:
            df_reciever(finger, num)
            handle_finger(finger)

        else:
            for i in range(4):
                SMILES_to_Morgan('targetmol_compound.smi', radius=i, useFeatures=False)
                #                 # SMILES_to_Morgan('./Test/Extra_drug.smi', radius=i, useFeatures=False, ofile = './Test/')
                SMILES_to_Morgan('targetmol_compound.smi', radius=i, useFeatures=True)
                # SMILES_to_Morgan('./Test/Extra_drug.smi', radius=i,  useFeatures=True, o



