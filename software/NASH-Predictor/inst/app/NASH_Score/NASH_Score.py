import numpy as np
import pandas as pd

#匹配FINDER
def cal_nash_score(clock_gene, KEY_PATHWAY):
    Finder_node=pd.read_table('./FINDER/code/results/StepRatio_0.0100/Finder_PPI.txt',header=0)#Finder分数
    Finder_node.columns=['node_renum']
    Finder_score=pd.read_table('./FINDER/code/results/StepRatio_0.0100/MaxCCList_Strategy_Finder_PPI.txt',header=None)#Finder分数
    Finder_score.columns=['Finder_score']
    Finder_score=pd.concat([Finder_node,Finder_score],axis=1)
    print(Finder_score)
    Finder_node_map=pd.read_table('./ppi_network/Finder_node_mapping.txt',sep='\t',header=None)#Finder节点匹配信息
   
    Finder_node_map.columns=['node_renum','GeneID']
    
    big_block_INFOR=pd.read_csv(r'./data/Human.GRCh38.p13.csv')#基因注释

    big_block_key=pd.merge(clock_gene,big_block_INFOR,on='GeneID')
    Finder_score=pd.merge(Finder_node_map,Finder_score,on="node_renum")#finder得分信息
    big_block_key=pd.merge(Finder_score,big_block_INFOR,on='GeneID')
    
    col_list = ['GeneID', 'Symbol', 'Finder_score','EnsemblGeneID']
    big_block_key=big_block_key.loc[:,col_list]
    
    #通路计算
    pathway_df = pd.read_table(r'./data/ontology_network.tsv', sep='\t')
    pathway_df_KEGG = pathway_df[pathway_df['DB'] == 'KEGG']

    KEY_gene_PATHWAY = pathway_df_KEGG[pathway_df_KEGG['term_id'].isin(list(KEY_PATHWAY['Term']))]
    sim_total = []
    for i in big_block_key['EnsemblGeneID']:
        pathway_num = pathway_df_KEGG[pathway_df_KEGG['gene_id'] == i]
        if pathway_num.shape[0] == 0:
            sim = 0
        else:
            sim = KEY_gene_PATHWAY[KEY_gene_PATHWAY['gene_id'] == i].shape[0] / pathway_num.shape[0]

        sim_total.append(sim)
    big_block_key['sim']=sim_total
    GAPR_results=pd.read_csv(r'GAPR_Results.csv')
    big_block_key_all=pd.merge(big_block_key,GAPR_results,left_on='GeneID',right_on='AP')
    big_block_key_all['Score_total'] = big_block_key_all['Finder_score'] * big_block_key_all['LCC_Changes'] * big_block_key_all[
        'sim']
    big_block_key_all=big_block_key_all.sort_values(by='Score_total',ascending=False)
    score_min = big_block_key_all['Score_total'].min()
    score_max = big_block_key_all['Score_total'].max()
    big_block_key_all['Score_total'] = (big_block_key_all['Score_total'] - score_min) / (score_max - score_min)
    big_block_key_all.to_csv(r'./NASH_Score/NASH_Score.csv',index=False)
    return big_block_key_all
