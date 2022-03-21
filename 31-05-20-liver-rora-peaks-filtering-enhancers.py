#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pandas as pd
import numpy as np


# In[3]:


names = 'chromosome-segment','segment-start','segment-end', 'segment'
segments = pd.read_csv('/media/serlia/Storage 2/results-chipseq/chromhmm-run2-learnmodel/TREGWT_10_segments.bed', delimiter='\t', names=names)


# In[4]:


segments.head()


# In[5]:


segments_annotated = segments.replace(to_replace = {'E6' : 'repressed',
                                                      'E8' : 'active-enhancer',
                                                      'E7' : 'primed-enhancer-open',
                                                      'E9' : 'primed-enhancer-histone',
                                                      'E10' : 'poised-enhancer',
                                                      'E3' : 'open-active',
                                                      'E4' : 'open',
                                                      'E1' : 'weakly-active'
                                                     }
                                        )


# In[22]:


pick_segments1 = segments_annotated['segment'] == 'open-active'
pick_segments2 = segments_annotated['segment'] == 'open'
pick_segments3 = segments_annotated['segment'] == 'weakly-active'
segments_to_analyse_active = segments_annotated[pick_segments1 | pick_segments2 | pick_segments3]


# In[23]:


segments_to_analyse_active.describe()


# In[39]:


segments_to_analyse_active.head()


# In[50]:


segments_to_analyse_active = segments_to_analyse_active.sort_values(by=['chromosome-segment', 'segment-start'])


# In[51]:


segments_to_analyse_active.to_csv('treg-segments-active-open-20-may.bed', sep = '\t', encoding = 'utf8', index = False)


# In[41]:


# 20 May RORa ChIP-seq signal analysis in conservative sites between Treg and liver
# Pick also repressed regions in Treg chromatin to analyse RORa binding
# Will take E6 (named 'repressed') and also E2 and E5 which have no marks at all and may be heterochromatin
# Here I don't need only differentially bound events

pick_segments1 = segments_annotated['segment'] == 'repressed'
pick_segments2 = segments_annotated['segment'] == 'E2'
pick_segments3 = segments_annotated['segment'] == 'E5'
#significant_db = intersect_annotated['4'] <= 0.05
segments_to_analyse_inactive = segments_annotated[pick_segments1 | pick_segments2 | pick_segments3]


# In[47]:


segments_to_analyse_inactive = segments_to_analyse_inactive.sort_values(by=['chromosome-segment', 'segment-start'])


# In[49]:


segments_to_analyse_inactive.to_csv('treg-segments-inactive-20-may.bed', sep = '\t', encoding = 'utf8', index = False)


# In[48]:


segments_to_analyse_inactive.head()


# In[52]:


pick_segments1 = segments_annotated['segment'] == 'repressed'

segments_to_analyse_repressed = segments_annotated[pick_segments1]
segments_to_analyse_repressed = segments_to_analyse_repressed.sort_values(by=['chromosome-segment', 'segment-start'])
segments_to_analyse_repressed.head()


# In[53]:


segments_to_analyse_repressed.to_csv('treg-segments-repressed-20-may.bed', sep = '\t', encoding = 'utf8', index = False)


# In[58]:


pick_segments1 = segments_annotated['segment'] == 'open-active'
pick_segments2 = segments_annotated['segment'] == 'weakly-active'
segments_to_analyse_acetyl = segments_annotated[pick_segments1 | pick_segments2]
segments_to_analyse_acetyl = segments_to_analyse_acetyl.sort_values(by=['chromosome-segment', 'segment-start'])
segments_to_analyse_acetyl.head()


# In[56]:


segments_to_analyse_acetyl.to_csv('treg-segments-acetyl-20-may.bed', sep = '\t', encoding = 'utf8', index = False)


# In[6]:


# 31-may-2020
# Now will pick enhancer regions to test RORa in enhancers in a bundle with microarray

pick_segments1 = segments_annotated['segment'] == 'active-enhancer'
pick_segments2 = segments_annotated['segment'] == 'primed-enhancer-open'
pick_segments3 = segments_annotated['segment'] == 'primed-enhancer-histone'
pick_segments4 = segments_annotated['segment'] == 'poised-enhancer'
segments_to_analyse_enhancers = segments_annotated[pick_segments1 | pick_segments2 | pick_segments3 | pick_segments4]
segments_to_analyse_enhancers = segments_to_analyse_enhancers.sort_values(by=['chromosome-segment', 'segment-start'])
segments_to_analyse_enhancers.head()


# In[7]:


segments_to_analyse_enhancers.to_csv('treg-segments-enhancers-31-may.bed', sep = '\t', encoding = 'utf8', index = False)


# In[ ]:




