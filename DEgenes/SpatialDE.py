#!/usr/bin/env python
# coding: utf-8

# In[2]:


get_ipython().run_line_magic('pylab', 'inline')
import pandas as pd

rcParams['axes.spines.right'] = False
rcParams['axes.spines.top'] = False

import NaiveDE
import SpatialDE


# In[3]:


counts = pd.read_csv('Exampledata/MOB/Rep11_MOB_0.csv', index_col=0)
counts = counts.T[counts.sum(0) >= 3].T  # Filter practically unobserved genes
sample_info = pd.read_csv('Exampledata/MOB/MOB_sample_info.csv', index_col=0)
counts = counts.loc[sample_info.index]


# In[5]:


#correct the library size or sequencing depth of the spatial samples 
norm_expr = NaiveDE.stabilize(counts.T).T
resid_expr = NaiveDE.regress_out(sample_info, norm_expr.T, 'np.log(total_counts)').T


# In[6]:


sample_resid_expr = resid_expr.sample(n=1000, axis=1, random_state=1)
X = sample_info[['x', 'y']]
results = SpatialDE.run(X, sample_resid_expr)
results.head().T


# In[7]:


results.sort_values('qval').head(10)[['g', 'l', 'qval']]
figsize(10, 3)
for i, g in enumerate(['Kcnh3', 'Pcp4', 'Igfbp2']):
    plt.subplot(1, 3, i + 1)
    plt.scatter(sample_info['x'], sample_info['y'], c=norm_expr[g]);
    plt.title(g)
    plt.axis('equal')
    plt.colorbar(ticks=[]);


# In[8]:


#Usually set qval<0.05. A very weak qval for this example
sign_results = results.query('qval < 0.5')
#avg lengthscale ~1.5. l=1.8 for some extra spatial covariance
#C: the number of patterns
sign_results['l'].value_counts()
histology_results, patterns = SpatialDE.aeh.spatial_patterns(X, resid_expr, sign_results, C=3, l=1.8, verbosity=1)
histology_results.head()


# In[10]:


figsize(10, 3)
for i in range(3):
    plt.subplot(1, 3, i + 1)
    plt.scatter(sample_info['x'], sample_info['y'], c=patterns[i]);
    plt.axis('equal')
    plt.title('Pattern {} - {} genes'.format(i, histology_results.query('pattern == @i').shape[0] ))
    plt.colorbar(ticks=[]);


# In[9]:


#to see what the coexpressed genes determining a histological pattern 
for i in histology_results.sort_values('pattern').pattern.unique():
    print('Pattern {}'.format(i))
    print(', '.join(histology_results.query('pattern == @i').sort_values('membership')['g'].tolist()))
    print()


# In[ ]:




