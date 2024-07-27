import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Load data from csv file
GNCommunityData = pd.read_csv('GNCommunity.csv').sort_values(by='Gate num')
GreedyData = pd.read_csv('GREEDY.csv').sort_values(by='Gate num')
KaHyParData = pd.read_csv('KaHyPar.csv').sort_values(by='Gate num')
PartitionScheme1Data = pd.read_csv('PartitionScheme1.csv').sort_values(by='Gate num')
PartitionScheme2Data = pd.read_csv('PartitionScheme2.csv').sort_values(by='Gate num')

df_dic = {
    'GNCommunity': GNCommunityData,
    'Greedy': GreedyData,
    'KaHyPar': KaHyParData,
    'PartitionScheme1': PartitionScheme1Data,
    'PartitionScheme2': PartitionScheme2Data
}

X = GNCommunityData['Gate num']

# output time column

for key, value in df_dic.items():
    plt.plot(X, value['Time'], label=key)

plt.xlabel('Gate num')
plt.ylabel('Time (ms)')
plt.title('Contraction Time')
plt.legend()
plt.savefig('ContractionTime.png')
plt.clf()

# output contraction cost column

for key, value in df_dic.items():
    plt.plot(X, value['Contraction cost (log)'], label=key)

plt.xlabel('Gate num')
plt.ylabel('Contraction cost (log)')
plt.title('Contraction Cost (for classical tensor contraction)')
plt.legend()
plt.savefig('ContractionCost.png')
plt.clf()

# output node num final column

for key, value in df_dic.items():
    plt.plot(X, value['Node num final'], label=key)

plt.xlabel('Gate num')
plt.ylabel('Node num final')
plt.title('Node num final for TDD')
plt.legend()
plt.ylim(0, 10000)
plt.savefig('NodeNumFinal.png')
plt.clf()

# output node num max column

for key, value in df_dic.items():
    plt.plot(X, value['Node num max'], label=key)

plt.xlabel('Gate num')
plt.ylabel('Node num max')
plt.title('Node num max for TDD')
plt.legend()
plt.ylim(0, 10000)
plt.savefig('NodeNumMax.png')
plt.clf()

# plot together

fig, axs = plt.subplots(2, 2, figsize=(15, 10))

# Plot Time
for key, value in df_dic.items():
    axs[0, 0].plot(X, value['Time'], label=key)
axs[0, 0].set_xlabel('Gate num')
axs[0, 0].set_ylabel('Time (ms)')
axs[0, 0].set_title('Contraction Time')

# Plot Contraction Cost
for key, value in df_dic.items():
    axs[0, 1].plot(X, value['Contraction cost (log)'], label=key)
axs[0, 1].set_xlabel('Gate num')
axs[0, 1].set_ylabel('Contraction cost (log)')
axs[0, 1].set_title('Contraction Cost (for classical tensor contraction)')

# Plot Node num final
for key, value in df_dic.items():
    axs[1, 0].plot(X, value['Node num final'], label=key)
axs[1, 0].set_xlabel('Gate num')
axs[1, 0].set_ylabel('Node num final')
axs[1, 0].set_title('Node num final for TDD')
axs[1, 0].set_ylim(0, 10000)  # Set y-axis range

# Plot Node num max
for key, value in df_dic.items():
    axs[1, 1].plot(X, value['Node num max'], label=key)
axs[1, 1].set_xlabel('Gate num')
axs[1, 1].set_ylabel('Node num max')
axs[1, 1].set_title('Node num max for TDD')
axs[1, 1].set_ylim(0, 10000)  # Set y-axis range

# Merge legend
handles, labels = axs[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', ncol=len(df_dic))

plt.tight_layout(rect=[0, 0.05, 1, 1])
plt.savefig('CombinedPlot.png')
plt.show()