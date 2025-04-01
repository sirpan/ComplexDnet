# Process_ppi.py
import pandas as pd
import os

def process_data(gene_path, gene_info_path, ppi_path, output_dir):
    # 读取数据
    gene = pd.read_csv(gene_path)
    gene_information = pd.read_csv(gene_info_path, sep='\t')
    PPI_2022 = pd.read_csv(ppi_path, sep='\t')
    
    print(f"初始基因数量: {len(gene)}")
    print(f"原始PPI网络边数: {len(PPI_2022)}")

    # 基因匹配（不区分大小写）
    gene_sc = gene_information[gene_information['Symbol'].str.upper().isin(gene['ensembl_id'].str.upper())]
    gene_sc = gene_sc.drop_duplicates(subset='GeneID')

    # 筛选PPI网络
    protein_a = PPI_2022['#Protein A'].isin(gene_sc['GeneID'])
    protein_b = PPI_2022['Protein B'].isin(gene_sc['GeneID'])
    filtered_ppi = PPI_2022[protein_a | protein_b]
    filtered_ppi = filtered_ppi.iloc[:, [0, 1]]  # 保留前两列

    # 新增步骤：在生成映射前删除自环边（基于原始ID）
    print(f"过滤自环边前边数: {len(filtered_ppi)}")
    filtered_ppi = filtered_ppi[filtered_ppi['#Protein A'] != filtered_ppi['Protein B']]
    print(f"过滤自环边后边数: {len(filtered_ppi)}")

    # 生成GAPR文件（保留原始ID）
    gapr_path = os.path.join(output_dir, "GAPR_PPI.csv")
    filtered_ppi.to_csv(gapr_path, index=False)

    # 生成新的节点映射（基于过滤后的数据）
    nodes = set(filtered_ppi['#Protein A']).union(set(filtered_ppi['Protein B']))
    node_mapping = {node: idx for idx, node in enumerate(sorted(nodes))}  # 排序保证可重复性

    # 写入映射文件
    mapping_path = os.path.join(output_dir, "Finder_node_mapping.txt")
    with open(mapping_path, 'w') as f:
        for node, idx in node_mapping.items():
            f.write(f"{idx}\t{node}\n")

    # 写入Finder网络文件（此时不需要再检查自环边）
    finder_path = os.path.join(output_dir, "Finder_PPI.txt")
    with open(finder_path, 'w') as f:
        for _, row in filtered_ppi.iterrows():
            node1 = node_mapping[row['#Protein A']]
            node2 = node_mapping[row['Protein B']]
            f.write(f"{node1} {node2} {{}}\n")  # 格式化为空格分隔

    # 统计信息
    print(f"最终节点数量: {len(nodes)}")
    print(f"有效边数量: {len(filtered_ppi)}")
    print("PPI处理完成！")
