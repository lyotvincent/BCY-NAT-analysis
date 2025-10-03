import matplotlib.pyplot as plt
import numpy as np

# 数据
data = [
    {
        'category': 'Epithelial',
        'go_terms': ['intermediate filament', 'epithelial cell differentiation', 'cell-cell junction'],
        'values': [20.88, 15.94, 8.74]
    },
    {
        'category': 'Endothelial',
        'go_terms': ['detection of chemical stimulus involved in sensory perception', 'vasculature development', 'miRNA-mediated post-transcriptional gene silencing'],
        'values': [47.30, 11.83, 9.06]
    },
    {
        'category': 'T cell',
        'go_terms': ['leukocyte activation', 'regulation of leukocyte activation', 'external side of plasma membrane'],
        'values': [55.97, 42.50, 30.38]
    },
    {
        'category': 'B cell',
        'go_terms': ['regulation of leukocyte activation', 'cell activation', 'adaptive immune response'],
        'values': [29.54, 25.84, 23.90]
    },
    {
        'category': 'Myeloid',
        'go_terms': ['immune response-regulating signaling pathway', 'inflammatory response', 'regulation of cell activation'],
        'values': [43.84, 43.82, 31.90]
    },
    {
        'category': 'Mast',
        'go_terms': ['detection of chemical stimulus involved in sensory perception of smell', 'regulation of inflammatory response', 'sialic acid binding'],
        'values': [21.41, 10.29, 7.23]
    },
    {
        'category': 'Plasma',
        'go_terms': ['cell activation', 'cytokine-mediated signaling pathway', 'regulation of cell activation'],
        'values': [14.44, 10.70, 10.03]
    },
    {
        'category': 'Fibroblast',
        'go_terms': ['extracellular matrix', 'miRNA-mediated post-transcriptional gene silencing', 'collagen trimer'],
        'values': [35.08, 28.49, 13.81]
    },
    {
        'category': 'PeriVascular',
        'go_terms': ['detection of chemical stimulus involved in sensory perception', 'regulation of system process', 'regulatory ncRNA-mediated post-transcriptional gene silencing'],
        'values': [27.39, 12.50, 9.92]
    }
]

# 设置条形图的宽度
bar_width = 0.02

# 创建一个图形和多个子图
fig, axes = plt.subplots(nrows=len(data), ncols=1, figsize=(10, 1 * len(data)), sharex=True)

# 遍历每个子图和数据集
for i, (ax, dataset) in enumerate(zip(axes, data)):
    category = dataset['category']
    go_terms = dataset['go_terms']
    values = dataset['values']
    
    # 按 values 降序排序
    sorted_terms = [term for _, term in sorted(zip(values, go_terms), reverse=True)]
    sorted_values = sorted(values, reverse=True)
    
    # 计算Y轴位置，缩短间距，并反转顺序
    y_positions = np.arange(len(sorted_terms))[::-1] * 0.03  # 反转顺序
    
    # 绘制条形图
    ax.barh(y_positions, sorted_values, height=bar_width, color='darkred')
    
    # 设置Y轴标签
    ax.set_yticks(y_positions)
    ax.set_yticklabels(sorted_terms)
    
    # 放大Y轴标签字号
    ax.tick_params(axis='y', labelsize=12)
    
    # 只保留第一个条形图的标题
    if i == 0:
        ax.set_title('Top 3 enriched GO terms')
    else:
        ax.set_title('')
    
    # 只保留最后一个条形图的X轴标签
    if i == len(data) - 1:
        ax.set_xlabel('-log10(P)')
    else:
        ax.set_xlabel('')
    
    # 显示X轴刻度
    ax.xaxis.set_visible(True)

# 调整子图间的间距
plt.tight_layout()

# 显示图表
plt.savefig('go_enrichment.png')