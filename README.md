**Author**: Alessandro Giudice    
**Contributor**: Samuel Santhosh Gomez  

# Complementarity in the binding region of a protein complex
Below I report the procedure for finding the contact zone in a complex of two proteins and estimating their similarity.  
The `text` written in this manner represents the variables of the code used, which is visible in the appendix or <a href="https://github.com/agds93/complexes_complementarity/tree/main/code" target="_blank">here</a>.    
Procedures regarding the patches on a protein surface can be found <a href="https://github.com/agds93/percentage_non_functionality/" target="_blank">here</a>.

## Search for the binding region
The entire surface of the studied protein complex is visible in Figure 0.

<p align="center"><img src="https://github.com/agds93/complexes_complementarity/blob/main/img/two_proteins_01.png" width=600px /></p><p align="center"><img src="https://github.com/agds93/complexes_complementarity/blob/main/img/two_proteins_02.png" width=600px /></p>
<p align="center"><i>Figure 0</i>: Protein A (blue) and Protein B (red) from two points of view.</p>

To find the contact zone between the two surfaces, use the `GroupNearPoints` function based on the selected threshold distance `Daa`. This function returns the two parts of the contact area:
* `patch_prot_a`, i.e. the contact patch on surface A (blue points in Figure 0).
* `patch_prot_b`, i.e. the contact patch on surface B (red points in Figure 0).

This function also provides the `center_a` and `center_b` indices closest to the center of mass (CoM) of the contact zone, respectively on the surface A and on the surface B. These points together with the two parts of the contact zone are visible in Figure 1.  

<p align="center"><img src="https://github.com/agds93/complexes_complementarity/blob/main/img/contact_zone_01.png" width=600px /></p>
<p align="center"><img src="https://github.com/agds93/complexes_complementarity/blob/main/img/contact_zone_02.png" width=600px /></p>
<p align="center"><i>Figure 1</i>: Contact area on surfaces A (above) and B (below). The red points are the closest ones to the CoM of the relative area.</p>

The Figure 2 represents the mean of the patch with center `center_a` obtained by the *Weights* method (left) and the *Projections* method (right). In the top row of  the figure there are the two original means created in the same way explained in the study of the percentage of non-functionality. Instead, in the bottom row of the figure, there are the relative processed means of the patch, in which number of pixels is increased and every empty pixel of original mean is filled with the average value of the neighboring pixels. This new mean returned by the function `ZernikeCoeff`.

<p align="center"><img src="https://github.com/agds93/complexes_complementarity/blob/main/img/ProteinA_Point6351.png" width=700px></p>
<p align="center"><img src="https://github.com/agds93/complexes_complementarity/blob/main/img/ProteinA_Point6351_processed.png" width=700px></p>
<p align="center"><i>Figure 2</i>: Original (top) and processed (bottom) mean of a patch of surface A with the two methods. Its center is the point 6351.</p>

The same things are represented in Figure 3 but referring to the patch with center `center_b`.

<p align="center"><img src="https://github.com/agds93/complexes_complementarity/blob/main/img/ProteinB_Point6026.png" width=700px></p>
<p align="center"><img src="https://github.com/agds93/complexes_complementarity/blob/main/img/ProteinB_Point6026_processed.png" width=700px></p>
<p align="center"><i>Figure 3</i>: Original (top) and processed (bottom) mean of a patch of surface B with the two methods. Its center is the point 6026.</p>

## Complementarity
Given two patches in the contact zone (the first in the surface A and the second in the surface B) to know how similar one is to the other there are two conditions:
* Patches must have normal versors facing in opposite directions so to be comparable.
* Patches must not have empty pixels (or islands) within the unit circle where the Zernike formalism is defined.

The first condition is satisfied if the matrix of original means (top row of Figure 2-3) are generated via the `PatchesMethods` function, where the first surface patch is facing up and the other is facing down. Instead, the second condition is met if the matrix of the processed mean (bottom row of Figure 2-3) returned by the function `ZernikeCoeff` is used.    
Hence the complementarity of two patches is estimated through `ZernikeCoeff_Distance` in particular by the difference `c_inv_diff` of the coefficient modules of the Zernike expansion between the two processed plans of patches. As is visible in Figure 2 and Figure 3, the two patches have a good similarity in fact the value of `c_inv_diff` is approximately equal to one.  
In Figure 4 are the visible values of `c_inv_diff` for one hundred points chosen randomly on `patch_prot_a` and `patch_prot_b` as a function of the difference in percentage of non-functionality between the two methods. 

<p align="center"><img src="https://github.com/agds93/complexes_complementarity/blob/main/img/inv_vs_perc.png" width=800px></p>
<p align="center"><i>Figure 4</i>: Difference of the invariants between the two methods for one hundred random points on each contact zone as a function of the relative difference of percentage of non-functionality.</p>

## Appendix
### Libraries and modules
The code written was executed with <a href="https://jupyterlab.readthedocs.io/en/stable/" target="_blank">JupyterLab</a> using `python 3.8`.  
The python modules used, including `jupyterlab`, installed via <a href="https://pip.pypa.io/en/stable/" target="_blank">pip</a>, are listed below.
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
The `mayavi` module, specifically `mlab`, is needed to display 3D surfaces in a Qt window, so as to produce Figure 0-1.  
While the basic libraries are
```python
sys.path.append("./bin/")
import ZernikeFunc as ZF
import SurfaceFunc as SF
```
written by <a href="https://scholar.google.it/citations?user=hjkTN0YAAAAJ&hl=it" target="_blank">Mattia Miotto</a>.
### Parameters
The parameter values used to select a patch and produce the fit plan are
```python
Npixel = 25    # the side of the plane in pixels
Dpp = 0.5      # the distance between points of the same patch
Rs = 6         # the radius of the sphere that includes the patch
threshold = 5  # threshold value to determine if the variance is high
ZOrder = 20    # order of the Zernike expansion
Daa = 3        # threshold distance to find contact points between surfaces
```
Values are in ångström, except `Npixel` and` ZOrder`.
### Load surfaces
To load the surface points of Protein A (available <a href="https://github.com/agds93/complexes_complementarity/blob/main/data/3B0F_A_min.dms" target="_blank">here</a>) and initialize the `Surface` class object is used
```python
surf_name_a = "./data/3B0F_A_min.dms"
surf_a_ = pd.read_csv(surf_name_a)
l_a = len(surf_a_["x"])
print("Npoints Protein A =", l_a)
surf_a = np.zeros((l_a, 6))
surf_a[:,:] = surf_a_[["x", "y", "z", "Nx", "Ny", "Nz"]]
surf_a_obj = SF.Surface(surf_a[:,:], patch_num = 0, r0 = Rs, theta_max = 45)
```
To load the surface points of Protein B (available <a href="https://github.com/agds93/complexes_complementarity/blob/main/data/3B0F_B_min.dms" target="_blank">here</a>) and initialize the `Surface` class object is used
```python
surf_name_b = "./data/3B0F_B_min.dms"
surf_b_ = pd.read_csv(surf_name_b) 
l_b = len(surf_b_["x"])
print("Npoints Protein B =", l_b)
surf_b = np.zeros((l_b, 6))
surf_b[:,:] = surf_b_[["x", "y", "z", "Nx", "Ny", "Nz"]]
surf_b_obj = SF.Surface(surf_b[:,:], patch_num = 0, r0 = Rs, theta_max = 45)
```
The graphs in Figure 0 are produced by
```python
res1, c = SF.ConcatenateFigPlots([surf_a_obj.surface[:,:3],surf_b_obj.surface[:,:3]])
SF.Plot3DPoints(res1[:,0], res1[:,1], res1[:,2], c, 0.3)
```
### Utilities
Given a set of points `points` and another point` P`, the following function returns the closest element of `points` to` P`.
```python
def PointNearPoint(points, P) :
    dist_points_P = np.sqrt( (points[:,0]-P[0])**2 + (points[:,1]-P[1])**2 + (points[:,2]-P[2])**2 )
    min_dist = np.amin(dist_points_P)
    point_near_P = int( np.where(dist_points_P == min_dist)[0] )
    return point_near_P
```
This function will be used in `GroupNearPoints` to find the closest point to the center of mass.
### Contact area search
Given two surface objects (`surf_a_obj` and` surf_b_obj`) and a threshold distance `Daa`, the following function returns:  
* the index of the point of surface A closest to the center of mass of its contact region.
* the matrix of the points of the contact zone on the surface A.
* the index of the point on surface B closest to the center of mass of its contact region.
* the matrix of the points of the contact zone on the surface B.

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
This function is used with
```python
center_a, patch_prot_a, center_b, patch_prot_b = GroupNearPoints(Daa, surf_a_obj, surf_b_obj)
```
However, the point indices `center_a` and` center_b` from the `GroupNearPoints` function refer to the` patch_prot_a` and `patch_prot_b` contact zones. To obtain the indices with respect to whole surface is used
```python
center_a_true = PointNearPoint(surf_a[:,:3], patch_prot_a[center_a])
center_b_true = PointNearPoint(surf_b[:,:3], patch_prot_b[center_b])
center_a = center_a_true
center_b = center_b_true
```
The graphs in Figure 1 are produced by respectively
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
### Creating Patches
The following function generates the averages of two patches (one per surface) with opposite orientations for each method.
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
Processed averages are obtained from the original averages via
```python
plane_W_a_proc = surf_a_obj.EnlargePixels( surf_a_obj.FillTheGap_everywhere(plane_=plane_W_a) )
plane_P_a_proc = surf_a_obj.EnlargePixels( surf_a_obj.FillTheGap_everywhere(plane_=plane_P_a) )
plane_W_b_proc = surf_a_obj.EnlargePixels( surf_a_obj.FillTheGap_everywhere(plane_=plane_W_b) )
plane_P_b_proc = surf_a_obj.EnlargePixels( surf_a_obj.FillTheGap_everywhere(plane_=plane_P_b) )
```
However this operation is included in the `ZernikeCoeff` function.
### Zernike coefficients
The following function returns the `plane_proc` processed plane, the list of Zernike coefficients` coeff` and their `coeff_inv` modules.
```python
def ZernikeCoeff(ZOrder, surf_a_obj, plane) :
    plane_proc = surf_a_obj.EnlargePixels( surf_a_obj.FillTheGap_everywhere(plane_=plane) )
    zernike_env = ZF.Zernike2d(plane_proc)
    coeff = zernike_env.ZernikeDecomposition(order=ZOrder)  # coeff is a list
    coeff_inv = np.absolute(coeff)
    return plane_proc, coeff, coeff_inv
```
The following function calculates the `c_inv_diff` difference between two lists of Zernike coefficient modules.
```python
def ZernikeCoeff_Distance(ZOrder, surf_a_obj, plane_1, surf_b_obj, plane_2) :
    _, _, c_inv_1 = ZernikeCoeff(ZOrder, surf_a_obj, plane_1)
    _, _, c_inv_2 = ZernikeCoeff(ZOrder, surf_b_obj, plane_2)
    c_inv_diff = np.sqrt( sum( (c_inv_1[:]-c_inv_2[:])**2 ) )
    return c_inv_diff
```
### Complementarity
The difference between the two methods of the modules of the Zernike coefficients for the point closest to the center of mass (CoM) of each contact zone, with the relative difference in percentage of non-functionality, is given by 
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
The difference of the invariants and of the percentage between the two hundred point methods of zone A and zone B is given by
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
The generated file is available <a href="https://github.com/agds93/complexes_complementarity/blob/main/data/inv_vs_perc.txt" target="_blank">here</a>.
### Charts
The following graphical function averages a patch produced with two different methods: `CreatePlane_Weigths` and` CreatePlane_Projections`. This function produces the graphs in Figure 2-3. The inputs are:
* the name of the protein to be included in the title.
* the number of `Npixel` pixels. 
* the radius `Rs` of the sphere that includes the patch.
* the plan generated by `CreatePlane_Weigths`.
* the plan generated by `CreatePlane_Projections`.
* the `center` index chosen as the center of the patch.
* the distance between the points of the `Dpp` patch.
* the `Daa` value of the threshold distance.
* the color maps to use.
* the name of the output file with the appropriate extension.

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
Figure 4 is produced by
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
