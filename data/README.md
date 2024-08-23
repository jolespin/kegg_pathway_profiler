# Database
* Version: v2024.4.27
* Filename: database.pkl.gz
* Description: Nested Python dictionary with keys as KEGG pathway modules
* Example: 

```python
from kegg_pathway_profiler.utils import read_pickle
database = read_pickle("data/database.pkl.gz")
database["M00001"]
{'graph': <networkx.classes.multidigraph.MultiDiGraph at 0x1ab02ccd0>,
 'ko_to_nodes': {'K00844': [[0, 2]],
  'K12407': [[0, 2]],
  'K00845': [[0, 2]],
  'K25026': [[0, 2]],
  'K00886': [[0, 2]],
  'K08074': [[0, 2]],
  'K00918': [[0, 2], [3, 4]],
  'K01810': [[2, 3]],
  'K06859': [[2, 3]],
  'K13810': [[2, 3]],
  'K15916': [[2, 3]],
  'K00850': [[3, 4]],
  'K16370': [[3, 4]],
  'K21071': [[3, 4]],
  'K01623': [[4, 5]],
  'K01624': [[4, 5]],
  'K11645': [[4, 5]],
  'K16305': [[4, 5]],
  'K16306': [[4, 5]],
  'K01803': [[5, 6]],
  'K00134': [[6, 8]],
  'K00150': [[6, 8]],
  'K00927': [[8, 7]],
  'K11389': [[6, 7]],
  'K01834': [[7, 9]],
  'K15633': [[7, 9]],
  'K15634': [[7, 9]],
  'K15635': [[7, 9]],
  'K01689': [[9, 10]],
  'K00873': [[10, 1]],
  'K12406': [[10, 1]]},
 'optional_kos': set(),
 'name': 'Glycolysis (Embden-Meyerhof pathway), glucose => pyruvate',
 'classes': 'Pathway modules; Carbohydrate metabolism; Central carbohydrate metabolism',
 'definition': '(K00844,K12407,K00845,K25026,K00886,K08074,K00918) (K01810,K06859,K13810,K15916) (K00850,K16370,K21071,K00918) (K01623,K01624,K11645,K16305,K16306) K01803 ((K00134,K00150) K00927,K11389) (K01834,K15633,K15634,K15635) K01689 (K00873,K12406)'}
```

