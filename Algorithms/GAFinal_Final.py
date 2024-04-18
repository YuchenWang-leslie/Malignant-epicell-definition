#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import random
import pandas as pd
import math
import anndata
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings("ignore")
from matplotlib import rcParams
from sklearn.mixture import GaussianMixture
from sklearn.metrics import silhouette_score


# In[ ]:


#improt a anndata object and a favorabel Differencial Expression Genes list
DEG = pd.read_csv('./DEG_new.csv',sep=',')
adata = anndata.read_h5ad("/home/liujingyue/wyc_20231021/single_cell_SEQ/project/chemotherapy/epicell.h5ad")
lenDEG=len(DEG)


# In[ ]:


length=100
#选取基因数量
size=20
toprange=500
crossr=0.8
mutater=0.1
generation=5
elitenumber=2


# In[ ]:


def getweights(lenDEG):
    #生成0到lenDEG的一共lenDEG个数
    numbers = np.arange(lenDEG)
    
    #计算每个数的权重，使用正态分布
    mean = round(lenDEG/2) #正态分布的均值
    std_dev = 200 #正态分布的标准差
    
    #计算每个数对应的概率，使中间的数具有更高的概率
    probabilities = np.exp(-((numbers - mean) / std_dev)**2)

    weight1=probabilities[mean:]
    weight2=probabilities[:mean]
    return weight1, weight2


# In[ ]:


def popinit(length,size,lenDEG,toprange):
    halflen=round(length/2)
    rownum=lenDEG
    weight1, weight2=getweights(lenDEG)
    #归一化，要求概率总和为1
    weight1=weight1[:toprange]/np.sum(weight1[:toprange])
    weight2=weight2[-toprange:]/np.sum(weight2[-toprange:])
    totallist=np.arange(lenDEG)
    
    popu=[]
    for i in range(size):
        chrome=[]
        chrome1=np.random.choice(totallist[:toprange], size=halflen, replace=False, p=weight1)
        chrome2=np.random.choice(totallist[-toprange:], size=halflen, replace=False, p=weight2)
        chrome=chrome1.tolist()+chrome2.tolist()
        popu.append(chrome)
    
    return popu  #各自从toprange范围内取前500，取值范围内随机抽取50个，两个染色体相加作为一个个体，重复形成种群


# In[ ]:


def getfitness(adata,popu,length,size,DEG):
    fitness=[]
    halflen=round(length/2)
    for i in range(len(popu)):
        chrome=popu[i]
        index1=chrome[:halflen-1]
        index2=chrome[halflen:]
        top=DEG.iloc[index1,[0]]
        top.columns=['genename']
        tail=DEG.iloc[index2,[0]]
        tail.columns=['genename']
        top=top['genename'].values.tolist()
        tail=tail['genename'].values.tolist()
        sc.tl.score_genes(adata,gene_list=top,score_name="malignant")
        sc.tl.score_genes(adata,gene_list=tail,score_name="normal")
        score1 = adata.obs['malignant']
        score2 = adata.obs['normal']

        data_matrix = np.vstack((score1, score2)).T

        gmm = GaussianMixture(n_components=2,covariance_type='full', random_state=20240403,n_init=30,max_iter=2000,tol= 1e-2)
        gmm.fit(data_matrix)
        # 获取每个样本点的归属度（属于每个聚类的概率）
        probs = gmm.predict_proba(data_matrix)
        # 找到归属度大于0.99的目标点
        interest_points = data_matrix[np.max(probs, axis=1) > 0.99]       
        interest_labels = gmm.predict(interest_points)
        cluster_centers = []
        
        for label in np.unique(interest_labels):  # 对于每个聚类标签
            cluster_points = interest_points[interest_labels == label]  # 获取属于当前聚类的所有感兴趣的点
            cluster_center = np.mean(cluster_points, axis=0)  # 计算该聚类的均值（聚类中心）
            cluster_centers.append(cluster_center)
        # 输出感兴趣的聚类中心坐标
        for i, center in enumerate(cluster_centers):
            dist = np.linalg.norm(cluster_centers[0] - cluster_centers[1])  # 计算欧氏距离 
            
        fit=0
        fit=silhouette_score(interest_points,interest_labels)+dist
        fitness.append(fit)

    print(i)
    return fitness


# In[ ]:


def select(popu,fitness,size,elite_index,elite):
    newpopu=popu
    # 先把精英拿出来
    newpopu=[newpopu[i] for i in range(len(newpopu)) if i not in elite_index]
    fitness=[fitness[i] for i in range(len(fitness)) if i not in elite_index]
    # 余下的部分进行权重抽样
    new_fitness = [item - min(fitness) for item in fitness]
    sum_weight=np.sum(new_fitness)
    weight=[x/sum_weight for x in new_fitness]
    print(len(newpopu),len(fitness),new_fitness)
    candidate_index=np.random.choice(range(len(newpopu)),size=len(newpopu),p=weight,replace=True)
    candidate=[newpopu[i] for i in candidate_index]
    
    candidate.extend(elite)
    
    newpopu=candidate
    
    return newpopu


# In[ ]:


def crossover(popu,crossr,length):
    newpopu=popu
    for i in range(0, len(popu), 2):
        if random.random()<crossr:
            cpoint=math.floor(random.random()*length)
            child1=popu[i][:cpoint]+popu[i+1][cpoint:]
            child2=popu[i+1][:cpoint]+popu[i][cpoint:]
            newpopu[i]=child1
            newpopu[i+1]=child2
    return newpopu


# In[ ]:


def mutate(popu,mutater,length,lenDEG):
    newpopu=popu
    rownum=round(lenDEG/2)
    halflen=round(length/2)
    for i in range(len(popu)):
        if random.random()<mutater:
            mpoint=math.floor(random.random()*length)
            print(f"变异行：{i}变异点：{mpoint}\n")
            if mpoint<halflen:
                newgene=math.floor(random.random()*lenDEG)
                while newgene in popu[i]:
                    newgene=math.floor(random.random()*lenDEG)
                print(f"变异原：{popu[i][mpoint]}变异后：{newgene}\n")
                popu[i][mpoint]=newgene
            else:
                newgene=lenDEG-math.floor(random.random()*lenDEG)
                while newgene in popu[i]:
                    newgene=lenDEG-math.floor(random.random()*lenDEG)
                print(f"变异原：{popu[i][mpoint]}变异后：{newgene}\n")
                popu[i][mpoint]=newgene
    return newpopu


# In[ ]:


def save_elite(popu,fitness,elitenumber):
    fitness_list = [(value, index) for index, value in enumerate(fitness)]
    fitness_list.sort(reverse=True, key=lambda x: x[0])
    elite_fitness=fitness_list[:elitenumber]
    elite_index=[tup[1] for tup in elite_fitness]
    elite = [popu[i] for i in elite_index]

    return elite_index,elite


# In[ ]:


def genetic_algorithm(length,size,DEG,adata,mutater,crossr,generation,lenDEG,toprange,elitenumber):
    #初始化种群
    popu=popinit(length,size,lenDEG,toprange)
    fitness=getfitness(adata,popu,length,size,DEG)
    #输出当前最优解
    bestsolution=popu[fitness.index(max(fitness))]
    finalsolution=bestsolution
    finalfitness=max(fitness)
    print(f"初始 Best Solution = {bestsolution}, Fitness_Silhouette = {finalfitness}\n")
    #初始精英选择
    elite,elite_index=save_elite(popu,fitness,elitenumber)
    #建一个空的list记录每次generation的maxfitness
    show_list=[]
    for gen in range(generation):
        
        newpopu=popu
        
        #选择
        newpopu=select(newpopu,fitness,size,elite,elite_index)
        print(f"选择完成{gen}\n")
        
        #突变
        newpopu=mutate(newpopu,mutater,length,lenDEG)
        print(f"变异完成{gen}\n")
        
        #交叉
        newpopu=crossover(newpopu,crossr,length)
        print(f"交叉完成{gen}\n")
        
        fitness=getfitness(adata,newpopu,length,size,DEG)
        
        #精英强制进入
        #精英选择
        elite,elite_index=save_elite(newpopu,fitness,elitenumber)
        print(f"精英选择完成{gen}\n")
        
        #输出当前最优解
        maxfitness=max(fitness)
        show_list.append(maxfitness)
        bestsolution=newpopu[fitness.index(maxfitness)]
        popu=newpopu
        
        if maxfitness>=finalfitness:
            finalfitness=maxfitness
            finalsolution=bestsolution
            print(f"新的 Generation {gen}: Best Solution = {bestsolution}, Fitness_Silhouette = {maxfitness}\n")
        else:
            print(f"稍差的 Generation {gen}: Best Solution = {bestsolution}, Fitness_Silhouette = {maxfitness}\n")
            
        print(f"循环完成{gen}\n")
    
    return finalsolution, finalfitness


# <h1>主程序

# In[ ]:


finalsolution, finalfitness=genetic_algorithm(length,size,DEG,adata,mutater,crossr,generation,lenDEG,toprange,elitenumber)
print(f"最好的基因选择list：{finalsolution},最好的silouette：{finalfitness}")

