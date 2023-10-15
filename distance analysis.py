#!/usr/bin/env python
# coding: utf-8

# In[1]:


from Bio import AlignIO
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


# ## 60seq:

# In[2]:


msafile = "./data/60seq.aln"
alignment = AlignIO.read(msafile, "clustal") 

columns_of_interest = [408, 450, 452, 456, 566, 687]
def get_nongapped_distances(alignment, columns_of_interest):
    distances = []
    seq_ids = []
    for record in alignment:
        seq_id = record.id
        seq = record.seq
        row_distances = []
        total_length = 0
        for i in range(len(columns_of_interest) - 1):
            row_distance = 0
            for j in range(columns_of_interest[i], columns_of_interest[i + 1]):
                if "-" not in seq[j]:
                    row_distance += 1
            row_distances.append(row_distance)
            total_length += row_distance
        distances.append(tuple(row_distances + [total_length]))
        seq_ids.append(seq_id.split("|")[1])
    return np.array(distances), np.array(seq_ids)
# nongapped_distances, seq_ids = get_nongapped_distances(alignment, columns_of_interest)
# for i in range(len(seq_ids)):
#     print( "408-450, 450-452, 452-456, 456-566, 566-687, total lenght--->", f"{seq_ids[i]}: {nongapped_distances[i]}")
# print("ACC: [ Q-H1 H1-H2 H2-E E-H3 H3-T Total]")
# for i in range(len(seq_ids)):
#     print(f"{seq_ids[i]}: {nongapped_distances[i]}")
distances, seq_ids = get_nongapped_distances(alignment, columns_of_interest)
data = {"ACC": seq_ids}
labels = ["Q-H1", "H1-H2", "H2-E", "E-H3", "H3-T"]
for i in range(len(columns_of_interest) - 1):
    data[labels[i]] = distances[:, i]
data["Total"] = distances[:, -1]
df1 = pd.DataFrame(data)
df1 = df1.drop(columns=["H1-H2"])
print(df1)

html = df1.to_html()
with open("df1.html", "w") as f:
    f.write(html)


# In[3]:


df1.boxplot()
plt.title('60 seq')
plt.xlabel('Columns')
plt.ylabel('Values')
# plt.show()


# ## DCAb short 50

# In[3]:


msafile = "./data/DCAb short 50.fasta"
alignment = AlignIO.read(msafile, "fasta")

columns_of_interest = [109,140,142,146,210,296]
def get_nongapped_distances(alignment, columns_of_interest):
    distances = []
    seq_ids = []
    for record in alignment:
        seq_id = record.id
        seq = record.seq
        row_distances = []
        total_length = 0
        for i in range(len(columns_of_interest) - 1):
            row_distance = 0
            for j in range(columns_of_interest[i], columns_of_interest[i + 1]):
                if "-" not in seq[j]:
                    row_distance += 1
            row_distances.append(row_distance)
            total_length += row_distance
        distances.append(tuple(row_distances + [total_length]))
        seq_ids.append(seq_id.split("|")[1])
    return np.array(distances), np.array(seq_ids)

# nongapped_distances, seq_ids = get_nongapped_distances(alignment, columns_of_interest)
# for i in range(len(seq_ids)):
#     print( "109-140,140-142,142-146,146-210,210-296, total lenght--->", f"{seq_ids[i]}: {nongapped_distances[i]}")
    
distances, seq_ids = get_nongapped_distances(alignment, columns_of_interest)
data = {"ACC": seq_ids}
labels = ["Q-H1", "H1-H2", "H2-E", "E-H3", "H3-T"]
for i in range(len(columns_of_interest) - 1):
    data[labels[i]] = distances[:, i]
data["Total"] = distances[:, -1]
df2 = pd.DataFrame(data)
df2 = df2.drop(columns=["H1-H2"])
print(df2)

html = df2.to_html()
with open("df2.html", "w") as f:
    f.write(html)


# In[5]:


df2.boxplot()
plt.title('50 seq')
plt.xlabel('Columns')
plt.ylabel('Values')
# plt.show()


# ## DCAb long 96

# In[12]:


msafile = "./data/DCAb long 96.fasta"
alignment = AlignIO.read(msafile, "fasta")

columns_of_interest = [109,140,142,146,210,296]
def get_nongapped_distances(alignment, columns_of_interest):
    distances = []
    seq_ids = []
    for record in alignment:
        seq_id = record.id
        seq = record.seq
        row_distances = []
        total_length = 0
        for i in range(len(columns_of_interest) - 1):
            row_distance = 0
            for j in range(columns_of_interest[i], columns_of_interest[i + 1]):
                if "-" not in seq[j]:
                    row_distance += 1
            row_distances.append(row_distance)
            total_length += row_distance
        distances.append(tuple(row_distances + [total_length]))
        seq_ids.append(seq_id.split("|")[1])
    return np.array(distances), np.array(seq_ids)

# nongapped_distances, seq_ids = get_nongapped_distances(alignment, columns_of_interest)
# for i in range(len(seq_ids)):
#     print("109-140,140-142,142-146,146-210,210-296, total lenght--->", f"{seq_ids[i]}: {nongapped_distances[i]}")
    
distances, seq_ids = get_nongapped_distances(alignment, columns_of_interest)
data = {"ACC": seq_ids}
labels = ["Q-H1", "H1-H2", "H2-E", "E-H3", "H3-T"]
for i in range(len(columns_of_interest) - 1):
    data[labels[i]] = distances[:, i]
data["Total"] = distances[:, -1]
df3 = pd.DataFrame(data)
df3 = df3.drop(columns=["H1-H2"])
# df3

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
print(df3)


# In[7]:


df3.boxplot()
plt.title('96 seq')
plt.xlabel('Columns')
plt.ylabel('Values')
# plt.show()


# ## DCAb 146 iter

# In[5]:


msafile = "./data/DCAb 146 iter.aln"
alignment = AlignIO.read(msafile, "clustal")

columns_of_interest = [109,140,142,146,210,296]
def get_nongapped_distances(alignment, columns_of_interest):
    distances = []
    seq_ids = []
    for record in alignment:
        seq_id = record.id
        seq = record.seq
        row_distances = []
        total_length = 0
        for i in range(len(columns_of_interest) - 1):
            row_distance = 0
            for j in range(columns_of_interest[i], columns_of_interest[i + 1]):
                if "-" not in seq[j]:
                    row_distance += 1
            row_distances.append(row_distance)
            total_length += row_distance
        distances.append(tuple(row_distances + [total_length]))
        seq_ids.append(seq_id.split("|")[1])
    return np.array(distances), np.array(seq_ids)

# nongapped_distances, seq_ids = get_nongapped_distances(alignment, columns_of_interest)
# for i in range(len(seq_ids)):
#     print("109-140,140-142,142-146,146-210,210-296, total lenght--->", f"{seq_ids[i]}: {nongapped_distances[i]}")
    
distances, seq_ids = get_nongapped_distances(alignment, columns_of_interest)
data = {"ACC": seq_ids}
labels = ["Q-H1", "H1-H2", "H2-E", "E-H3", "H3-T"]
for i in range(len(columns_of_interest) - 1):
    data[labels[i]] = distances[:, i]
data["Total"] = distances[:, -1]
df4 = pd.DataFrame(data)
df4 = df4.drop(columns=["H1-H2"])
print(df4)

html = df4.to_html()
with open("df4.html", "w") as f:
    f.write(html)


# In[9]:


df4.boxplot()
plt.title('DCAb 146 iter')
plt.xlabel('Columns')
plt.ylabel('Values')
# plt.show()


# ## bact ACA 494

# In[6]:


msafile = "./data/bact ACA Hosein 494.fasta"
alignment = AlignIO.read(msafile, "fasta")

columns_of_interest = [472,585,587,591,605,747]
def get_nongapped_distances(alignment, columns_of_interest):
    distances = []
    seq_ids = []
    for record in alignment:
        seq_id = record.id
        seq = record.seq
        row_distances = []
        total_length = 0
        for i in range(len(columns_of_interest) - 1):
            row_distance = 0
            for j in range(columns_of_interest[i], columns_of_interest[i + 1]):
                if "-" not in seq[j]:
                    row_distance += 1
            row_distances.append(row_distance)
            total_length += row_distance
        distances.append(tuple(row_distances + [total_length]))
        seq_ids.append(seq_id.split("|")[1])
    return np.array(distances), np.array(seq_ids)

# nongapped_distances, seq_ids = get_nongapped_distances(alignment, columns_of_interest)
# for i in range(len(seq_ids)):
#     print("472,585,587,591,605,747, total lenght--->", f"{seq_ids[i]}: {nongapped_distances[i]}")
    
distances, seq_ids = get_nongapped_distances(alignment, columns_of_interest)
data = {"ACC": seq_ids}
labels = ["Q-H1", "H1-H2", "H2-E", "E-H3", "H3-T"]
for i in range(len(columns_of_interest) - 1):
    data[labels[i]] = distances[:, i]
data["Total"] = distances[:, -1]
df5 = pd.DataFrame(data)
df5 = df5.drop(columns=["H1-H2"])
print(df5)

html = df5.to_html()
with open("df5.html", "w") as f:
    f.write(html)


# In[11]:


df5.boxplot()
plt.title('bact ACA 494')
plt.xlabel('Columns')
plt.ylabel('Values')
# plt.show()


# In[14]:


msafile = "./data/hum_mus_zfish_droso aligned.fasta" # metazoan data
alignment = AlignIO.read(msafile, "fasta")

columns_of_interest = [189,272,274,286,299,399]
def get_nongapped_distances(alignment, columns_of_interest):
    distances = []
    seq_ids = []
    for record in alignment:
        seq_id = record.id
        seq = record.seq
        row_distances = []
        total_length = 0
        for i in range(len(columns_of_interest) - 1):
            row_distance = 0
            for j in range(columns_of_interest[i], columns_of_interest[i + 1]):
                if "-" not in seq[j]:
                    row_distance += 1
            row_distances.append(row_distance)
            total_length += row_distance
        distances.append(tuple(row_distances + [total_length]))
        seq_ids.append(seq_id.split("|")[1])
    return np.array(distances), np.array(seq_ids)
# nongapped_distances, seq_ids = get_nongapped_distances(alignment, columns_of_interest)
# for i in range(len(seq_ids)):
#     print( "408-450, 450-452, 452-456, 456-566, 566-687, total lenght--->", f"{seq_ids[i]}: {nongapped_distances[i]}")
# print("ACC: [ Q-H1 H1-H2 H2-E E-H3 H3-T Total]")
# for i in range(len(seq_ids)):
#     print(f"{seq_ids[i]}: {nongapped_distances[i]}")
distances, seq_ids = get_nongapped_distances(alignment, columns_of_interest)
data = {"ACC": seq_ids}
labels = ["Q-H1", "H1-H2", "H2-E", "E-H3", "H3-T"]
for i in range(len(columns_of_interest) - 1):
    data[labels[i]] = distances[:, i]
data["Total"] = distances[:, -1]
df6 = pd.DataFrame(data)
df6 = df6.drop(columns=["H1-H2"])
print(df6)
df6.boxplot()
plt.title('metazoan data - 58 seq')
plt.xlabel('Columns')
plt.ylabel('Values')
# plt.show()

# html = df6.to_html()
# with open("df6.html", "w") as f:
#     f.write(html)


# In[13]:


msafile = "./data/UP ACA 994 renamed.fasta" 
alignment = AlignIO.read(msafile, "fasta")

columns_of_interest = [814,941,943,947,962,1080]
def get_nongapped_distances(alignment, columns_of_interest):
    distances = []
    seq_ids = []
    for record in alignment:
        seq_id = record.id
        seq = record.seq
        row_distances = []
        total_length = 0
        for i in range(len(columns_of_interest) - 1):
            row_distance = 0
            for j in range(columns_of_interest[i], columns_of_interest[i + 1]):
                if "-" not in seq[j]:
                    row_distance += 1
            row_distances.append(row_distance)
            total_length += row_distance
        distances.append(tuple(row_distances + [total_length]))
        seq_ids.append(seq_id[1:])
    return np.array(distances), np.array(seq_ids)
# nongapped_distances, seq_ids = get_nongapped_distances(alignment, columns_of_interest)
# for i in range(len(seq_ids)):
#     print( "408-450, 450-452, 452-456, 456-566, 566-687, total lenght--->", f"{seq_ids[i]}: {nongapped_distances[i]}")
# print("ACC: [ Q-H1 H1-H2 H2-E E-H3 H3-T Total]")
# for i in range(len(seq_ids)):
#     print(f"{seq_ids[i]}: {nongapped_distances[i]}")
distances, seq_ids = get_nongapped_distances(alignment, columns_of_interest)
data = {"ACC": seq_ids}
labels = ["Q-H1", "H1-H2", "H2-E", "E-H3", "H3-T"]
for i in range(len(columns_of_interest) - 1):
    data[labels[i]] = distances[:, i]
data["Total"] = distances[:, -1]
df7 = pd.DataFrame(data)
df7 = df7.drop(columns=["H1-H2"])
# print(df11)
df7.boxplot()
plt.title('994 seq')
plt.xlabel('Columns')
plt.ylabel('Values')
# plt.show()

# html = df7.to_html()
# with open("df7.html", "w") as f:
#     f.write(html)

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
print(df7)


# In[14]:


# fig, axs = plt.subplots(1, 5, figsize=(20, 5))

# df1.boxplot(ax=axs[0])
# axs[0].set_title('60 seq')
# axs[0].set_xlabel('Columns')
# axs[0].set_ylabel('Values')

# df2.boxplot(ax=axs[1])
# axs[1].set_title('50 seq')
# axs[1].set_xlabel('Columns')
# axs[1].set_ylabel('Values')

# df3.boxplot(ax=axs[2])
# axs[2].set_title('96 seq')
# axs[2].set_xlabel('Columns')
# axs[2].set_ylabel('Values')

# df4.boxplot(ax=axs[3])
# axs[3].set_title('146 seq')
# axs[3].set_xlabel('Columns')
# axs[3].set_ylabel('Values')

# df5.boxplot(ax=axs[4])
# axs[4].set_title('494 seq')
# axs[4].set_xlabel('Columns')
# axs[4].set_ylabel('Values')
# plt.show()


# In[15]:


# fig, axs = plt.subplots(1, 5, figsize=(15, 5), sharex=True, sharey=True)
# for df, ax in zip([df1, df2, df3, df5, df6], axs):
#     df.boxplot(ax=ax)
#     ax.set_title(f"{len(df)} seq")
#     ax.set_xlabel('Columns')
#     ax.set_ylabel('Values')
# fig.suptitle('Distances between key residues in four groups of sequences')
# plt.show()


# In[16]:


# fig, axs = plt.subplots(1, 5, figsize=(12, 5), sharex=True, sharey=True)
# for df, ax in zip([df2, df3, df1, df7,df6], axs):
#     sns.boxplot(data=df, ax=ax, showfliers=False)
#     ax.set_title(f"{len(df)} seq")
#     ax.set_xlabel('Columns')
#     ax.set_ylabel('Values')
# fig.suptitle('Distances between key residues in four groups of sequences')
# plt.show()


# In[17]:


fig, axs = plt.subplots(1, 5, figsize=(12, 5), sharex=True, sharey=True)
for df, ax in zip([df2, df3, df1, df7,df6], axs):
    sns.boxplot(data=df, ax=ax, showfliers=False)
    ax.set_title(f"{len(df)} seq")
    ax.set_xlabel('Residues')
    ax.set_ylabel('Distance')
    ax.grid(True, linestyle='--', linewidth=0.5, color='gray')  # major gridlines
    # Add minor gridlines
    ax.minorticks_on()
    ax.grid(True, which='minor', linestyle=':', linewidth=0.5, color='lightgray')

fig.suptitle('Distances between key residues in four groups of sequences')
plt.show()


# In[18]:


fig, axs = plt.subplots(1, 5, figsize=(12, 5), sharex=True, sharey=True)
titles = ["Bacterial DCA short", "Bacterial DCA long", "Eukaryotic DCA", "Bacterial ACA", "Metazoan ACA"]
for df, ax, title in zip([df2, df3, df1, df7, df6], axs, titles):
    sns.boxplot(data=df, ax=ax, showfliers=False)
    ax.set_title(title)
    ax.set_xlabel('Residues')
    ax.set_ylabel('Distance')
    ax.grid(True, linestyle='--', linewidth=0.5, color='gray')  # major gridlines
    # Add minor gridlines
    ax.minorticks_on()
    ax.grid(True, which='minor', linestyle=':', linewidth=0.5, color='lightgray')

fig.suptitle('Distances between key residues in five groups of sequences')
plt.show()


# In[ ]:




