Cell Assembly Toolbox
===========

This repo collects the code used to assess cell assembly properties for auditory
cortical neurons.

One of the great challenges in neuroscience is to get insight from the activity of large populations of neurons. We have know that population coding can increase information transfer and protect sensory information, but there have been very few methods to make sense of the data that provide a systematic analysis coupled with appropriate statistical hypotheseis tests. 

This repo collects functions for this type of analysis. Our goal is to search for groups of neurons that are active together - they are part of a functional response network. This groups of neurons have been termed cell assemblies. For cell assembly analysis in the auditory cortex, we recorded many neurons simultaneously with multi-channel recording probes. On average, this gives 60-90 neurons, recorded over a time period of
10-15 minutes.

The goal is to determine which sets of neurons within the 90 are synchronously active. Thus, we wish to see which subnetworks are present among the 90. The first step is to
spike sort the data, then bin the spike trains into 5 ms bins. The spike trains are then z-scored to normalized by the magnitude of the response. Next, the spike trains are collected into a large matrix, the correlation matrix is obtained, and the matrix is decomposed using PCA. Finally, the individual cell assemblies are obtained by applying ICA to the significant PCA components.

This procedure finds the cell assemblies. Now we want to go further. What do the cell assemblies tell us about the stimulus? Which stimulus features trigger synchronous activity? Do assemblies transfer more information compared to non cell assembly spikes? This toolbox contains Matlab functions to perform these analyses. The toolbox does not contain code to extract the assemblies - this can be found in theoretical papers by Vitor Dos Santos.

