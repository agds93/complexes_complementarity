**Autore**: Alessandro Giudice    
**Collaboratore**: Samuel Santhosh Gomez  

# Complementarità nella regione di legame di un complesso proteico
Di seguito riporto la procedura per trovare la zona di contatto in un complesso di due proteine e stimare la loro similarità.  
Il `testo` scritto in questa maniera rappresenta le variabili del codice usato, che è visibile in appendice o qui.  
I metodi per selezionare una patch, calcolare media e varianza per ogni pixel in disco unitario con due metodi, e produrre i rispettivi grafici sono riportati <a href="https://github.com/agds93/percentage_non_functionality/" target="_blank">qui</a>.  

## Ricerca della regione di legame
L'intera superficie del complesso proteico studiato è visibile in Figura 0.

<p align="center"><img src="img/two_proteins_01.png" width=700px /></p>
<p align="center"><img src="img/two_proteins_02.png" width=700px /></p>
<p align="center"><i>Figura 0</i>: Proteina A (blu) e proteina B (rosso) da due punti di vista.</p>

Per trovare la zona di contatto tra le due superfici si utilizza la funzione `GroupNearPoints` in base alla distanza di soglia `Daa` scelta. Tale funzione restituisce le due parti della zona di contatto:
* `patch_prot_a`, cioè la zona di contatto sulla superficie A (punti blu in Figura 0).
* `patch_prot_b`, cioè la zona di contatto sulla superficie B (punti rossi in Figura 0).

Tale funzione fornisce anche gli indici `center_a` e `center_b` più vicini al centro di massa (CoM) della zona di contatto, rispettivamente sulla superficie A e sulla superficie B. Tali punti insieme alle due parti della zona di contatto sono visibili in Figura 1.

<p align="center"><img src="img/contact_zone_01.png" width=700px /></p>
<p align="center"><img src="img/contact_zone_02.png" width=700px /></p>
<p align="center"><i>Figura 1</i>: Zona di contatto sulle superfici A (sopra) e B (sotto). I punti rossi sono i più vicini al CoM della relativa zona.</p>

I grafici in Figura 2 rappresentano la media della patch con centro `center_a` ottenuta con due metodi, dove, nei grafici in alto, è presente la media originale. Nella stessa figura, ma nei grafici in basso, è presente la media processata della patch, in cui i pixels sono incrementati e riempiti per rimuovere le aree vuote presenti nella media originale.

<p align="center"><img src="img/ProteinA_Point6351.png" width=700px></p>
<p align="center"><img src="img/ProteinA_Point6351_processed.png" width=700px></p>
<p align="center"><i>Figura 2</i>: Media originale (in alto) e processata (in basso) di una patch della superficie A.</p>

Le stesse cose sono rappresentate nella Figura 3 ma riferite alla patch con centro `center_b`.

<p align="center"><img src="img/ProteinB_Point6026.png" width=700px></p>
<p align="center"><img src="img/ProteinB_Point6026_processed.png" width=700px></p>
<p align="center"><i>Figura 3</i>: Media originale (in alto) e processata (in basso) di una patch della superficie B.</p>

## Complementarietà  
Date due patches nella zona di contatto (la prima nella superficie A e la seconda nella superficie B) si vuole sapere quanto sono simili tra loro.  
Per tale scopo:
* Le patches devono avere i versori normali rivolti in direzioni opposte, così da essere confrontabili.
* Le patches non devono avere pixels vuoti (o isole) all'interno del cerchio unitario in cui è definito il formalismo di Zernike.

La prima condizione è soddisfatta se si generano le medie tramite la funzione `PatchesMethods`, dove la patch della prima superficie è rivolta verso l'alto mentre l'altra verso il basso. Invece la seconda condizione è rispettata se si utilizzano per il calcolo dei coefficienti di Zernike la versione processata della media originale, cioè con i pixel incrementati e senza isole, come quelli della parte bassa della Figura 2-3, fornita dalla funzione `ZernikeCoeff` insieme ai coefficienti e ai loro moduli. L'idea è di riempire i pixels vuoti con la media dei pixels vicini.  
La complementarietà delle patches si stima tramite `ZernikeCoeff_Distance` in particolare dalla differenza `c_inv_diff` dei moduli dei coefficienti dell'espansione di Zernike tra i rispettivi piani processati delle due patch. Le patches da cui si ricavano i grafici in Figura 2 e Figura 3 hanno come centro rispettivamente `center_a` e `center_b`, cioè il punto più vicino al centro di massa di tale zona. Di conseguenza tali patch hanno una buona similarità, infatti il valore della differenza `c_inv_diff` tra le due rispettive liste di coefficienti di Zernike è pari a un numero vicino a uno.  
In Figura 4 è visibile la differenza degli invarianti `c_inv_diff` tra il metodo Weights e il metodo Projections per cento punti scelti casualmente su `patch_prot_a` e `patch_prot_b` in funzione della differenza di percentuale di non funzionalità tra i due metodi.

<p align="center"><img src="img/inv_vs_perc.png" width=800px></p>
<p align="center"><i>Figura 4</i>: Differenza degli invarianti tra i due metodi per cento punti su ogni zona di contatto in funzione della differenza di percentuale.</p>

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
Il modulo `mayavi`, in particolare `mlab`, è necessario per visualizzare le superfici 3D in una finestra Qt, così da produrre la Figura 0-1.  
Mentre le librerie di base sono
```python
sys.path.append("./bin/")
import ZernikeFunc as ZF
import SurfaceFunc as SF
```
scritte da <a href="https://scholar.google.it/citations?user=hjkTN0YAAAAJ&hl=it" target="_blank">Mattia Miotto</a>.
### Parametri
I valori dei parametri usati per selezionare una patch e produrre il piano di fit sono
```python
Npixel = 25    # il lato del piano in pixel
Dpp = 0.5      # la distanza tra i punti della stessa patch
Rs = 6         # il raggio della sfera che include la patch
threshold = 5  # valore soglia per stabilire se la varianza è alta
ZOrder = 20    # ordine dell'espansione di Zernike
Daa = 3        # distanza di soglia per trovare punti di contatto tra le superfici
```
I valori sono in ångström, tranne `Npixel` e `ZOrder`.
### Caricare le superfici
Per caricare i punti della superficie della proteina A (disponibile 
<a href="data/3B0F_A_min.dms" target="_blank">qui</a>
) ed inizializzare l'oggetto di classe `Surface` si usa
```python
surf_name_a = "./data/3B0F_A_min.dms"
surf_a_ = pd.read_csv(surf_name_a)
l_a = len(surf_a_["x"])
print("Npoints Protein A =", l_a)
surf_a = np.zeros((l_a, 6))
surf_a[:,:] = surf_a_[["x", "y", "z", "Nx", "Ny", "Nz"]]
surf_a_obj = SF.Surface(surf_a[:,:], patch_num = 0, r0 = Rs, theta_max = 45)
```
Per caricare i punti della superficie della proteina B (disponibile 
<a href="data/3B0F_B_min.dms" target="_blank">qui</a>
) ed inizializzare l'oggetto di classe `Surface` si usa
```python
surf_name_b = "./data/3B0F_B_min.dms"
surf_b_ = pd.read_csv(surf_name_b) 
l_b = len(surf_b_["x"])
print("Npoints Protein B =", l_b)
surf_b = np.zeros((l_b, 6))
surf_b[:,:] = surf_b_[["x", "y", "z", "Nx", "Ny", "Nz"]]
surf_b_obj = SF.Surface(surf_b[:,:], patch_num = 0, r0 = Rs, theta_max = 45)
```
I grafici in Figura 0 sono prodotti da
```python
res1, c = SF.ConcatenateFigPlots([surf_a_obj.surface[:,:3],surf_b_obj.surface[:,:3]])
SF.Plot3DPoints(res1[:,0], res1[:,1], res1[:,2], c, 0.3)
```
### Utilità
Dato un insieme di punti `points` e un altro punto `P`, la seguente funzione restituisce l'elemento di `points` più vicino a `P`.
```python
def PointNearPoint(points, P) :
    dist_points_P = np.sqrt( (points[:,0]-P[0])**2 + (points[:,1]-P[1])**2 + (points[:,2]-P[2])**2 )
    min_dist = np.amin(dist_points_P)
    point_near_P = int( np.where(dist_points_P == min_dist)[0] )
    return point_near_P
```
Questa funzione sarà usata in `GroupNearPoints` per trovare il punto più vicino al centro di massa.  
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
Tale funzione viene utilizzata con
```python
center_a, patch_prot_a, center_b, patch_prot_b = GroupNearPoints(Daa, surf_a_obj, surf_b_obj)
```
Però gli indici dei punti `center_a` e `center_b` dalla funzione `GroupNearPoints` sono riferiti alle zone di contatto `patch_prot_a` e `patch_prot_b`. Per ottenere gli indici rispetto alle superfici intere si utilizza
```python
center_a_true = PointNearPoint(surf_a[:,:3], patch_prot_a[center_a])
center_b_true = PointNearPoint(surf_b[:,:3], patch_prot_b[center_b])
center_a = center_a_true
center_b = center_b_true
```
I grafici in Figura 1 sono prodotti rispettivamente da
```python
cm = patch_prot_a[center_a]
res1, c = SF.ConcatenateFigPlots([patch_prot_a,np.row_stack([cm,cm])])
SF.Plot3DPoints(res1[:,0], res1[:,1], res1[:,2], c, 0.3)
```
```python
cm = patch_prot_b[center_b]
res1, c = SF.ConcatenateFigPlots([patch_prot_b,np.row_stack([cm,cm])])
SF.Plot3DPoints(res1[:,0], res1[:,1], res1[:,2], c, 0.3)
```
### Creazione delle patch
La seguente funzione genera le medie di due patch (una per superficie) con orientazioni opposte per ogni metodo.
```python
def PatchesMethods(Npixel, surf_a_obj, c_a, surf_b_obj, c_b, Dpp) :
    patch_a, _ = surf_a_obj.BuildPatch(point_pos=c_a, Dmin=Dpp)
    rot_patch_a, rot_patch_nv_a = surf_a_obj.PatchReorientNew(patch_a, +1)
    z_pa = surf_a_obj.FindOrigin(rot_patch_a)
    plane_W_a, _, _, _ = CreatePlane_Weigths("mean", patch=rot_patch_a, z_c=z_pa, Np=Npixel)
    plane_P_a, _, _, _ = CreatePlane_Projections("mean", patch=rot_patch_a, z_c=z_pa, Np=Npixel)
    patch_b, _ = surf_b_obj.BuildPatch(point_pos=c_b, Dmin=Dpp)
    rot_patch_b, rot_patch_nv_b = surf_b_obj.PatchReorientNew(patch_b, -1)
    z_pb = surf_b_obj.FindOrigin(rot_patch_b)
    plane_W_b, _, _, _ = CreatePlane_Weigths("mean", patch=rot_patch_b, z_c=z_pb, Np=Npixel)
    plane_P_b, _, _, _ = CreatePlane_Projections("mean", patch=rot_patch_b, z_c=z_pb, Np=Npixel)
    return plane_W_a, plane_P_a, plane_W_b, plane_P_b
```
Le medie processate si ottengono a partire dalle medie originali tramite
```python
plane_W_a_proc = surf_a_obj.EnlargePixels( surf_a_obj.FillTheGap_everywhere(plane_=plane_W_a) )
plane_P_a_proc = surf_a_obj.EnlargePixels( surf_a_obj.FillTheGap_everywhere(plane_=plane_P_a) )
plane_W_b_proc = surf_a_obj.EnlargePixels( surf_a_obj.FillTheGap_everywhere(plane_=plane_W_b) )
plane_P_b_proc = surf_a_obj.EnlargePixels( surf_a_obj.FillTheGap_everywhere(plane_=plane_P_b) )
```
Comunque tale operazione è compresa nella funzione `ZernikeCoeff`.
### Coefficienti di Zernike
La seguente funzione restituisce il piano processato `plane_proc`, la lista dei coefficienti di Zernike `coeff` e i loro moduli `coeff_inv`.
```python
def ZernikeCoeff(ZOrder, surf_a_obj, plane) :
    plane_proc = surf_a_obj.EnlargePixels( surf_a_obj.FillTheGap_everywhere(plane_=plane) )
    zernike_env = ZF.Zernike2d(plane_proc)
    coeff = zernike_env.ZernikeDecomposition(order=ZOrder)  # coeff is a list
    coeff_inv = np.absolute(coeff)
    return plane_proc, coeff, coeff_inv
```
La seguente funzione che calcola la differenza `c_inv_diff` tra due liste di moduli di coefficienti di Zernike.
```python
def ZernikeCoeff_Distance(ZOrder, surf_a_obj, plane_1, surf_b_obj, plane_2) :
    _, _, c_inv_1 = ZernikeCoeff(ZOrder, surf_a_obj, plane_1)
    _, _, c_inv_2 = ZernikeCoeff(ZOrder, surf_b_obj, plane_2)
    c_inv_diff = np.sqrt( sum( (c_inv_1[:]-c_inv_2[:])**2 ) )
    return c_inv_diff
```
### Complementarietà
La differenza tra i due metodi dei moduli dei coefficienti di Zernike per il punto più vicino al centro di massa (CoM) di ogni zona di contatto, con la relativa differenza di percentuale di non funzionalità, è data da 
```python
center_1 = center_a
center_2 = center_b
plane_W_1, plane_P_1, plane_W_2, plane_P_2 = MU.PatchesMethods(Npixel, surf_a_obj, center_1, surf_b_obj, center_2, Dpp)
coeff_diff_a = MU.ZernikeCoeff_Distance(ZOrder, surf_a_obj, plane_W_1, surf_a_obj, plane_P_1)
coeff_diff_b = MU.ZernikeCoeff_Distance(ZOrder, surf_b_obj, plane_W_2, surf_b_obj, plane_P_2)
_, perc_W_a = MU.PercHigherVariance_Weights("var", Npixel, surf_a_obj, center_1, Dpp, threshold)
_, perc_P_a = MU.PercHigherVariance_Projections("var", Npixel, surf_a_obj, center_1, Dpp, threshold)
_, perc_W_b = MU.PercHigherVariance_Weights("var", Npixel, surf_b_obj, center_2, Dpp, threshold)
_, perc_P_b = MU.PercHigherVariance_Projections("var", Npixel, surf_b_obj, center_2, Dpp, threshold)
perc_a = np.absolute( perc_W_a - perc_P_a )
perc_b = np.absolute( perc_W_b - perc_P_b )
print("Protein A: Patch center = {}".format(center_1))
print("Protein B: Patch center = {}".format(center_2))
print("Patch A: Difference Zernike coefficients = {},  perc_a = {}".format(coeff_diff_a,perc_a))
print("Patch B: Difference Zernike coefficients = {},  perc_b = {}".format(coeff_diff_b,perc_b))
```
La differenza degli invarianti e della percentuale tra i due metodi per cento punti della zona A e della zona B è data da
```python
with open("inv_vs_perc.txt", "w") as last_file :
    for i in range(100) :
        center_1 = np.random.randint(5000,7000)
        center_2 = np.random.randint(5000,7000)
        plane_W_a, plane_P_a, plane_W_b, plane_P_b = MU.PatchesMethods(Npixel, surf_a_obj, center_1, surf_b_obj, center_2, Dpp)
        coeff_diff_a = MU.ZernikeCoeff_Distance(ZOrder, surf_a_obj, plane_W_a, surf_a_obj, plane_P_a)
        coeff_diff_b = MU.ZernikeCoeff_Distance(ZOrder, surf_b_obj, plane_W_b, surf_b_obj, plane_P_b)
        _, perc_W_a = MU.PercHigherVariance_Weights("var", Npixel, surf_a_obj, center_1, Dpp, threshold)
        _, perc_P_a = MU.PercHigherVariance_Projections("var", Npixel, surf_a_obj, center_1, Dpp, threshold)
        _, perc_W_b = MU.PercHigherVariance_Weights("var", Npixel, surf_b_obj, center_2, Dpp, threshold)
        _, perc_P_b = MU.PercHigherVariance_Projections("var", Npixel, surf_b_obj, center_2, Dpp, threshold)
        perc_a = np.absolute( perc_W_a - perc_P_a )
        perc_b = np.absolute( perc_W_b - perc_P_b )
        last_file.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(center_1, center_2, coeff_diff_a, coeff_diff_b, perc_a, perc_b))
```
Il file generato è disponibile <a href="data/inv_vs_perc.txt" target="_blank">qui</a>.
### Grafici
La seguente funzione grafica la media di una patch prodotta con due diversi metodi: `CreatePlane_Weigths` e `CreatePlane_Projections`. Tale funzione produce i grafici in Figura 2-3. Gli input sono:
* il nome della proteina da inserire nel titolo.
* il numero di pixel `Npixel`. 
* il raggio `Rs` della sfera che include la patch.
* il piano generato da `CreatePlane_Weigths`.
* il piano generato da `CreatePlane_Projections`.
* l'indice `center` scelto come centro della patch.
* la distanza tra i punti della patch `Dpp`.
* il valore `Daa` della distanza di soglia.
* le mappe dei <a href="https://matplotlib.org/stable/tutorials/colors/colormaps.html" target="_blank">colori</a> da utilizzare.
* il nome del file di output con opportuna estensione.
```python
def PlotPatchesComparison(obj_name, Npixel, Rs, p_W, p_P, center, Dpp, Daa, color_map, name) :
    if obj_name == "" :
        obj_name = "Unknown"
    if len(color_map) != 1 :
        color_map = "Greens"
    matrix = [ p_W, p_P ]
    s0 = "Protein {} with threshold distance = {}\nPatch of center = {}, distance between points = {}, radius = {}".format(obj_name, Daa, center, Dpp, Rs)
    s1 = "Mean with Weights Method in function of position\n\n"
    s2 = "Mean with Projections Method in function of position\n\n"
    titles = [ s1, s2 ]  
    fig, ax = mpl.subplots(nrows=1, ncols=2, figsize=(8,4), dpi=200, facecolor="white")  # dpi=200 per compensare rasterized
    fig.suptitle(s0, fontsize="9")
    for row in range(1) :
        for col in range(2):
            data = matrix[col]
            ax[col].set_title(titles[col], fontsize="8")
            ax[col].set_xlabel("x", fontsize="8")
            ax[col].set_ylabel("y", fontsize="8")
            ax[col].tick_params(axis="both", width ="0.30", color="black", labelsize="8")
            for side in ax[col].spines.keys():  # 'top', 'bottom', 'left', 'right'
                ax[col].spines[side].set_linewidth(0.30)
                ax[col].spines[side].set_color("black")
            im = ax[col].pcolormesh(data, cmap=color_map, rasterized=True)   # senza rasterized il file è troppo grande
            ticks_list = [np.amin(data), np.amax(data)]
            cb = mpl.colorbar(im, ax=ax[col], ticks=ticks_list)
            cb.ax.tick_params(axis="both", width ="0.30", color="black", labelsize="8")
            for side in cb.ax.spines.keys():  # 'top', 'bottom', 'left', 'right'
                cb.ax.spines[side].set_linewidth(0.30)
                cb.ax.spines[side].set_color("black")
    fig.tight_layout()
    if name != "" or name == "default" :
        if name == "default" :
            n = "Mean_Patch{}.pdf".format(center, Dpp)
        else :
            n = name
        mpl.savefig("{}".format(n))
        print("The figure was generated.")
```
La Figura 4 è prodotta da
```python
coeff_a = np.loadtxt("./risultati/inv_vs_perc.txt", usecols=2, unpack=True)
coeff_b = np.loadtxt("./risultati/inv_vs_perc.txt", usecols=3, unpack=True)
perc_a = np.loadtxt("./risultati/inv_vs_perc.txt", usecols=4, unpack=True)
perc_b = np.loadtxt("./risultati/inv_vs_perc.txt", usecols=5, unpack=True)

fig, ax = mpl.subplots(nrows=1, ncols=1, figsize=(8,4), facecolor="white", dpi=200)
ax.set_xlabel("Difference of percentage", fontsize="8")
ax.set_ylabel("Difference of invariants", fontsize="8")
ax.tick_params(axis="both", width ="0.60", color="black", labelsize="6")
ax.locator_params(axis="x", nbins=21)
ax.locator_params(axis="y", nbins=21)
for side in ax.spines.keys():  # 'top', 'bottom', 'left', 'right'
    ax.spines[side].set_linewidth(0.60)
    ax.spines[side].set_color("black")
ax.plot(perc_a, coeff_a, "o", markersize="1", label="Zone A", color="blue", rasterized=True)
ax.plot(perc_b, coeff_b, "o", markersize="1", label="Zone B", color="red", rasterized=True)
ax.legend(fontsize="7")
fig.tight_layout()
mpl.savefig("inv_vs_perc.png")
```
