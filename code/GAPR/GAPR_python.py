import pandas as pd
import networkx as nx
import os
from itertools import chain

def save_network_to_file(G, round_number, base_path=r"D:\metabolic network\PPI\GAPR"):
    """
    将网络的边列表保存到文本文件中
    """
    file_name = f"network_round_{round_number}.txt"
    file_path = os.path.join(base_path, file_name)
    with open(file_path, 'w') as file:
        for edge in G.edges():
            file.write(f"{edge[0]} {edge[1]}\n")
    return file_path

def find_remove_articulation_points_in_AP_list(G):
    """
    按轮次查找并移除割点（AP节点）及其相关的边
    返回值：按轮次记录的移除割点的列表
    """
    AP_list = []
    list_components_change = []
    list_lcc_change = []

    saved_files = []
    round_number = 1
    round_list = []

    while True:
        # 使用NetworkX库中的函数来找到所有的割点
        articulation_points = list(nx.articulation_points(G))

        print(len(articulation_points))

        round_list.append(len(articulation_points)*[round_number])

        # 如果没有割点，结束循环
        if not articulation_points:
            break

        # 记录这一轮的割点
        AP_list.append(articulation_points)
        temp_components_change,temp_lcc_change = find_articulation_points_and_calculate_components_change(G,articulation_points)
        list_components_change.append(temp_components_change)
        list_lcc_change.append(temp_lcc_change)

        # 移除这些割点及其相关的边
        G.remove_nodes_from(articulation_points)
        saved_file_path = save_network_to_file(G, round_number)
        saved_files.append(saved_file_path)
        round_number += 1

    return AP_list,list_components_change,list_lcc_change,round_list


def calculate_connected_components(G):
    """
    计算给定网络的连通分量数量
    """
    return nx.number_connected_components(G)


def largest_connected_component_size(G):
    """
    计算给定网络的最大连通分量（LCC）的大小
    """
    if nx.is_empty(G):
        return 0
    largest_cc = max(nx.connected_components(G), key=len)
    return len(largest_cc)


# 在移除割点之前和之后计算LCC的大小，并计算（1-LCCnew/LCC)
def find_articulation_points_and_calculate_components_change(G, articulation_points):
    """
    查找割点并计算移除每个割点前后连通分量数量的变化及 1-LCCnew/LCC
    返回值：割点列表及其对应的连通分量变化和 1-LCCnew/LCC
    """
    list_components_change = []
    list_lcc_change = []

    list_test = []
    for point in articulation_points:
        original_components = calculate_connected_components(G)
        original_lcc_size = largest_connected_component_size(G)

        # 需要提前建立网络的copy
        G_copy = G.copy()
        # 移除割点
        G_copy.remove_node(point)
        new_components = calculate_connected_components(G_copy)
        new_lcc_size = largest_connected_component_size(G_copy)

        # 计算连通分量的变化和 1-LCCnew/LCC值
        components_change = new_components - original_components
        lcc_change = (1 - new_lcc_size / original_lcc_size) if original_lcc_size else 0

        # 将结果添加到列表中
        list_components_change.append(components_change)
        list_lcc_change.append(lcc_change)

        # 重新添加割点以保持网络状态;别扯淡了
        # G.add_node(point)

        # list_test.append((point,original_components,new_components))
        # list_test.append((point,original_lcc_size,new_lcc_size))

    # print(list_test)
    return list_components_change,list_lcc_change



# df = pd.read_csv('PPI_AD_324.txt', sep='\t')
df = pd.read_csv(r'D:\metabolic network\PPI\GAPR\NASH_PPI.txt', sep='\t')

# Create a graph from the DataFrame
G = nx.from_pandas_edgelist(df, 'Node1', 'Node2')

AP_list,list_components_change,list_lcc_change,round_list = find_remove_articulation_points_in_AP_list(G)


# print(round_list)
# 输出每轮移除的割点
# print(AP_list_of_articulation_points)

flattened_AP_list = list(chain.from_iterable(AP_list))
flattened_list_connected_component = list(chain.from_iterable(list_components_change))
flattened_list_lcc_change= list(chain.from_iterable(list_lcc_change))
flattened_round_list = list(chain.from_iterable(round_list))

df_final = pd.DataFrame({'AP':flattened_AP_list,'layer':flattened_round_list,'connected component changes':flattened_list_connected_component,'LCC changes':flattened_list_lcc_change})

# df_final.to_excel('test_AD_AP_GAPR_python.xlsx',index=0)
df_final.to_excel(r'D:\metabolic network\PPI\GAPR\test_2022NASH_AP.xlsx',index=0)

