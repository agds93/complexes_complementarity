# Complementarità di un complesso proteico nella regione di legame  
Di seguito riporto la procedura per trovare la zona di contatto in un complesso di due proteine e stimare il grado di similarità tra loro.  
Il `testo` scritto in questa maniera rappresenta le variabili del codice usato, visibile in appendice.  

## Appendice
### Librerie e moduli
Il codice scritto è stato eseguito con <a href="https://jupyterlab.readthedocs.io/en/stable/" target="_blank">JupyterLab</a> utilizzando `python 3.8`.  
I moduli python usati, compreso `jupyterlab`, installati tramite 
<a href="https://pip.pypa.io/en/stable/" target="_blank">pip</a>, sono elencati sotto.
```python
import os, sys
import numpy as np
import matplotlib.pyplot as mpl
import scipy as sp
import pandas as pd
```
```python
from mayavi import mlab
```
Il modulo `mayavi`, in particolare `mlab`, è necessario per visualizzare le superfici 3D in una finestra Qt.  
Mentre le librerie di base sono
```python
sys.path.append("./bin/")
import ZernikeFunc as ZF
import SurfaceFunc as SF
```
scritte da <a href="https://scholar.google.it/citations?user=hjkTN0YAAAAJ&hl=it" target="_blank">Mattia Miotto</a>.
