import pandas as pd
import networkx as nx
from itertools import chain

def gapr_analysis(ppi_path, output_dir):
    """执行GAPR核心分析流程"""
    # 读取PPI网络
    df = pd.read_csv(ppi_path, sep='\t')
    G = nx.from_pandas_edgelist(df, '#Protein A', 'Protein B')
    
    # 执行割点分析
    AP_list, components_changes, lcc_changes, rounds = [], [], [], []
    current_round = 1
    
    while True:
        articulation_points = list(nx.articulation_points(G))
        if not articulation_points:
            break
            
        # 计算网络变化
        comp_changes, lcc_chgs = calculate_network_changes(G, articulation_points)
        
        # 记录结果
        AP_list.extend(articulation_points)
        components_changes.extend(comp_changes)
        lcc_changes.extend(lcc_chgs)
        rounds.extend([current_round] * len(articulation_points))
        
        # 更新网络
        G.remove_nodes_from(articulation_points)
        current_round += 1
    
    # 生成结果数据框
    out_results= pd.DataFrame({
        'AP': AP_list,
        'Layer': rounds,
        'Component_Changes': components_changes,
        'LCC_Changes': lcc_changes
    })
    out_results.to_csv(output_dir, index=False)

def calculate_network_changes(G, points):
    """计算网络指标变化"""
    comp_changes, lcc_chgs = [], []
    original_lcc = len(max(nx.connected_components(G), key=len, default=[]))
    
    for point in points:
        # 复制网络以避免修改原图
        G_copy = G.copy()
        G_copy.remove_node(point)
        
        # 计算连通分量变化
        new_comp = nx.number_connected_components(G_copy)
        comp_changes.append(new_comp - nx.number_connected_components(G))
        
        # 计算LCC变化
        new_lcc = len(max(nx.connected_components(G_copy), key=len, default=[]))
        lcc_chgs.append(1 - new_lcc/original_lcc if original_lcc else 0)
    
    return comp_changes, lcc_chgs
