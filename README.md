Malignant-epithelial cells defination  
==========  
## A novel method identifying malignant cells based on GMM and Swarm intelligence algorithm  
### We use ***GA***,***PSO***,***SSA*** for the SIA. And GMM for clustering.   
#### **For** using this code, you need   
- A single cell readable object: **seurat object** or **adata**
- A bulk count matrix

---
1. Use the bulk counts to get a **original** DEG gene list.
2. Use the original DEG gene list's top and tail to score each cell based on gene expression matrix.
3. Get a matrix which is 2*cell numbers.
4. Assume two different type of cells.
5. Use Findmarker function get a **prossed** DEG gene list.

---
### Then you can input sc data and prossed DEG list in this application.
