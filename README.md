**Autore**: Alessandro Giudice    
**Collaboratore**: Samuel Santhosh Gomez  

# Complementarità di un complesso proteico nella regione di legame  
Di seguito riporto la procedura per trovare la zona di contatto in un complesso di due proteine e stimare la loro similarità.  
Il `testo` scritto in questa maniera rappresenta le variabili del codice usato, visibile in appendice.  

## Ricerca della regione di legame
L'intera superficie del complesso proteico studiato è visibile in Figura 0.
<p align="center"><img src="img/two_proteins_01.png" width=700px></p>
<p align="center"><img src="img/two_proteins_02.png" width=700px></p>
<p align="center"><i>Figura 0</i>: Proteina A (blu) e proteina B (rosso) da due punti di vista.</p>

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
### Caricare le superfici
Per caricare i punti della superficie della proteina A ed inizializzare l'oggetto di classe `Surface` si usa
```python
surf_name_a = "./data/3B0F_A_min.dms"
surf_a_ = pd.read_csv(surf_name_a)
l_a = len(surf_a_["x"])
print("Npoints Protein A =", l_a)
surf_a = np.zeros((l_a, 6))
surf_a[:,:] = surf_a_[["x", "y", "z", "Nx", "Ny", "Nz"]]
surf_a_obj = SF.Surface(surf_a[:,:], patch_num = 0, r0 = Rs, theta_max = 45)
```
Per caricare i punti della superficie della proteina B ed inizializzare l'oggetto di classe `Surface` si usa
```python
surf_name_b = "./data/3B0F_B_min.dms"
surf_b_ = pd.read_csv(surf_name_b) 
l_b = len(surf_b_["x"])
print("Npoints Protein B =", l_b)
surf_b = np.zeros((l_b, 6))
surf_b[:,:] = surf_b_[["x", "y", "z", "Nx", "Ny", "Nz"]]
surf_b_obj = SF.Surface(surf_b[:,:], patch_num = 0, r0 = Rs, theta_max = 45)
```
### Utilità
Dato un insieme di punti `points` e un altro punto `P`, la seguente funzione restituisce l'elemento di `points` più vicino `P`.
```python
def PointNearPoint(points, P) :
    dist_points_P = np.sqrt( (points[:,0]-P[0])**2 + (points[:,1]-P[1])**2 + (points[:,2]-P[2])**2 )
    min_dist = np.amin(dist_points_P)
    point_near_P = int( np.where(dist_points_P == min_dist)[0] )
    return point_near_P
```
Questa funzione sarà usata in `GroupNearPoints` per trovare il punto della zona di contatto più vicino al centro di massa.  
### Ricerca zona di contatto
Dati due oggetti superficie (`surf_a_obj` e `surf_b_obj`) e una distanza di soglia `Daa`, la seguente funzione restituisce:  
* l'indice del punto della superficie A più vicino al centro di massa della propria regione di contatto.
* la matrice dei punti della zona di contatto sulla superficie A.
* l'indice del punto della superficie B più vicino al centro di massa della propria regione di contatto.
* la matrice dei punti della zona di contatto sulla superficie B.
```python
def GroupNearPoints(Daa, surf_a_obj, surf_b_obj) :
    
    prot_A = surf_a_obj.surface[:,:3]
    prot_B = surf_b_obj.surface[:,:3]
    print("Research of contact points. This step requires time...")
    patch_prot_a, patch_prot_b = SF.ContactPoints(prot_A, prot_B, Daa)
    print("Research complete.")
    
    cm_a = np.mean(patch_prot_a[:,:3], axis=0)
    print("CM of protein A group =", cm_a)
    center_a = PointNearPoint(patch_prot_a[:,:3], cm_a)
    print("Patch protein A: Center = {} with coord = {}".format(center_a, patch_prot_a[center_a,:3]))
    
    cm_b = np.mean(patch_prot_b[:,:3], axis=0)
    print("CM of protein B group =", cm_b)
    center_b = PointNearPoint(patch_prot_b[:,:3], cm_b)
    print("Patch protein B: Center = {} with coord = {}".format(center_b, patch_prot_b[center_b,:3]))
    
    return center_a, patch_prot_a[:,:3], center_b, patch_prot_b[:,:3]
```
