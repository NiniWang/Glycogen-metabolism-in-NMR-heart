import os
import pandas as pd

# 获取文件夹的名称
folder_names = ['hypo_mouse_B_1', 'hypo_mouse_B_2', 'hypo_mouse_B_3', 'hypo_mouse_H_1', 'hypo_mouse_H_2', 'hypo_mouse_H_3', 'hypo_mouse_K_1', 'hypo_mouse_K_2', 'hypo_mouse_K_3', 'hypo_mouse_L_1', 'hypo_mouse_L_2', 'hypo_mouse_L_3', 'norm_mouse_B_1', 'norm_mouse_B_2', 'norm_mouse_B_3', 'norm_mouse_H_1', 'norm_mouse_H_2', 'norm_mouse_H_3', 'norm_mouse_K_1', 'norm_mouse_K_2', 'norm_mouse_K_3', 'norm_mouse_L_1', 'norm_mouse_L_2', 'norm_mouse_L_3']

# 读取每个文件的 quant.sf 文件，只选择 "NumReads" 列，并为 DataFrame 指定列名
dataframes = []
for folder_name in folder_names:
    file_path = os.path.join(folder_name, 'quant.sf')
    df = pd.read_csv(file_path, sep='\t', usecols=['Name', 'NumReads'], index_col='Name')
    dataframes.append(df)

# 合并文件，并为每个 DataFrame 指定列名
merged_data = pd.concat(dataframes, axis=1, keys=folder_names)

# 将合并的数据保存到新文件
merged_data.to_csv('merged_mouse_NumReads.csv', sep='\t')

