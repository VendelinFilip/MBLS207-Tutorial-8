# CLUSTERING EXERCISE - THE GRIFTPARK BACTERIUM

If you live in Utrecht, you have probably been to the Griftpark at least once. It is a recreational park, where people play football or frisbee, skate, or hang out next to the little lake. But do you know about the history of the Griftpark? From ~1860 till 1960, the Griftpark was the home of a municipal gas factory (figure 1). This gas factory dumped most of its waste directly into the ground, leading to major pollution of the soil, predominantly with aromatic hydrocarbons. Aromatic hydrocarbons are toxic pollutants, many of them being highly mutagenic/carcinogenic, and they are commonly found in crude oil and natural gas deposits (as well as emissions from cars and planes).

&lt;img&gt;Aerial view of the Griftpark in 1967 (panel A), and 2018 (panel B). The white outline is the approximate area that is isolated by the underground bentonite wall.&lt;/img&gt;

Figure 1. Aerial view of the Griftpark in 1967 (panel A), and 2018 (panel B). The white outline is the approximate area that is isolated by the underground bentonite wall.

---


## Page 3

Don’t worry, it is now completely safe to spend time in the Griftpark. The pollution was discovered when plans were made to turn the site into a recreational park, and the most heavily polluted soil was sealed off. For this purpose, a bentonite wall was built that reaches, on average, 65m into the ground (figure 1, white line). Additionally, there is a 1.5m thick bentonite wall on top, which is why the Griftpark is a little elevated from the streets around it.

To prevent the pollution from spreading into surrounding groundwater, the groundwater below the Griftpark is kept at constant low pressure. This is done with the Griftpark Groundwater Treatment Pipeline (GWTP, Figure 2). The GWTP is a pipeline that constantly pumps groundwater out of the Griftpark and into a small water treatment plant three kilometers away. There, the water is cleaned up, and the clean water gets treated once more with the rest of Utrecht’s wastewater. The GWTP does not clean the water by chemical means, but biologically. The GWTP only contains microorganisms that come in with the water and adds oxygen to the water to create ideal conditions for the microorganisms to degrade the pollutants. We have known for years that the GWTP is very effective at cleaning the water, but until recently, we didn’t know which microorganisms were responsible for the degradation of the aromatic hydrocarbons.

&lt;img&gt;Figure 2. A. Schematic overview of the Griftpark Groundwater Treatment Pipeline (GWTP). Groundwater from the contaminated soil (ft30) is gathered in a collection basin, which is infused with small amounts of oxygen (p2) and pumped over 3km (about three hours) to a treatment plant. The groundwater enters via the influent buffer (p4), where it is inoculated with water from the sludge thickener (p7), which is located at the end of the treatment plant. Then, the water proceeds to the aeration tank (p5), where the water is oxygenated and most of the contaminants get removed. After passing into an overflow (p6), the water then enters the settling tank, where water is left to settle. The sludge phase is pumped into the sludge thickener (p7), while the water phase undergoes treatment via sand filters (SF) and activated carbon (AC) before entering the sludge thickener as well. Finally, a portion of the sludge from the sludge thickener is pumped back into the influent buffer, while the remainder is sufficiently remediated to allow discharge into the communal wastewater influent. Seven sampled points are indicated with colored arrows.&lt;/img&gt;

In today’s Tutorial, we will work with the data from the Griftpark to get more information about different types of microorganisms that are able to degrade pollutants (i.e. break them down to non-toxic molecules), such as aromatic hydrocarbons/oil. This data is a very good example of environmental metagenomic data: samples were taken from the environment, in this case groundwater, DNA from all microorganisms was extracted, and the extracted DNA was sequenced. You will use some bioinformatics to find patterns and clusters in the data and try to understand how the samples from different points are similar or different.

The first, and one of the most important parts of any bioinformatics study is getting to know your dataset.

a. Download the file griftpark_5taxa.csv [here](griftpark_5taxa.csv) and load it into a jupyter notebook.

```python
import pandas as pd
import numpy as np
import seaborn as sns
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram

gp = pd.read_csv('griftpark_5taxa.csv')
```

b. This file contains a portion of the samples that we took from the Griftpark and, for each sample, the most abundant (=most present) bacterial taxa in the dataset. Take a look at the data frame by typing `gp`. How many rows and columns are there?

c. For which sampling point was the most DNA sequenced? Hint: you can generate a matrix of this dataset's values by doing `gp_data = gp.drop(gp.columns[0], axis=1)`, and then use `sum(axis=1)` over this new object. Does this tell you something about the sampling point?

d. Plot the data using seaborn. Which sampling points have the most similar microbiome? Can you already see which taxon the “Griftpark bacterium” belongs to?

```python
sns.heatmap(gp_data)
```

e. What does the following line do? (Hint: Print the new matrix and its row sums.)

```python
gp_norm = gp_data.div(gp_data.sum(axis=1), axis=0)
```

f. Plot a heatmap of the matrix `gp_norm`. Does it differ from the one in d?

g. Write a function that calculates the Manhattan distance between two sampling points. Use your function to calculate the Manhattan distance between sampling points:
i. p1 and p2
ii. p1 and p6
iii. p3 and p4
iv. p4 and p7

h. Execute the following line:

```python
dman = pdist(gp_norm, metric='cityblock')

---


## Page 5

to calculate the distances between the points using the Manhattan distances.

i. Are the Manhattan distances the same as the ones you calculated?

j. Now, calculate the Euclidean distances between the sampling points by changing the metric argument. Save the Euclidean distances in a new matrix. Are the Manhattan distances between the sampling point pairs always larger than the Euclidean distances? (Hint: use > expression and remember that it will return true or false for every item in a matrix/vector)

k. Take the Euclidean distance matrix. Cluster the sampling points using the function `linkage` and save the result into a new data structure with the name “cgp” (this name is important for an exercise below). Display the clustering using the `dendogram()` function. Try again by using `method='complete'` in `linkage`. Do the results differ? Write this clustering using a bracket-notation.

l. Which sampling points have the most similar microbiome? Compare the result with your answer under d, and with the Euclidean distance matrix under j. Look at Figure 2 and try to come up with a biological explanation for how the microbiomes from the sampling points inside the treatment plant cluster together. P3 and p4 are located within 20m of each other. Provide a biological explanation for why their microbiomes don’t cluster together using Figure 2.

m. Now that we know how to do the clustering in python, we can easily cluster the bacteria as well. Look at the following commands and explain per command what happens. (Hint: the function `.T` transposes a matrix.)

```python
gpt = gp_data.T
gpt_norm = gpt.div(gpt.sum(axis=1), axis=0)
dgpt_norm = pdist(gpt_norm, metric='euclidean')
cbac = linkage(dgpt_norm, method='complete')
dendrogram(cbac, color_threshold=0)
```

n. What are the two groups of bacteria that cluster closer together? Why do you think they do? Why could *Epsilonproteobacteria* be so separate in the dendrogram? Provide biological explanations using only the data you have (Hint: Groundwater is typically anaerobic).

o. Finally, we can make a bi-clustered heatmap of this same data as follows. Look at the clustering of rows and columns. Do they look familiar? Compare this with the heatmaps and clusterings you made above.

```python
sns.clustermap(gp_norm, cmap='coolwarm', row_linkage=cgp, col_linkage=cbac, figsize=(8,8))
```

Today, you have worked with real data from a research project done on the Griftpark. You have calculated different distance measures between samples and seen that the samples from the treatment plant have short distances to each other (=are more similar) and large

---


## Page 6

distances to the samples from the park. You have made a heatmap and found that the “Griftpark bacterium” most likely belongs to the class *Betaproteobacteria*, which has a low abundance in the park but dominates the microbial community inside the treatment plant (it would be more correct to say “Griftpark bacteria”, because it is actually a number of different *Betaproteobacteria* that degrade the pollutants). You have seen that *Epsilonproteobacteria*, which are the dominating class of bacteria in the park, are farthest away from the other taxa in the dendrogram. In the beginning of the Griftpark research project, we took the exact same steps to try and find the Griftpark bacterium. Some of the *Betaproteobacteria* that we found inside the treatment plant were already well-known, but some of them were completely new!
