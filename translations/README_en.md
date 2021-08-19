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
<p align="center"><i>Figure 1</i>: Contact area on surfaces A (above) and B (below). The red points are the closest to the CoM of the relative area.</p>

The graphs in Figure 2 represent the average of the patch with center `center_a` obtained by two methods, where, in the graphs above, the original average is present. In the same figure, but in the graphs below, there is the processed average of the patch, in which the pixels are increased and filled to remove the empty areas present in the original media.

<p align="center"><img src="https://github.com/agds93/complexes_complementarity/blob/main/img/ProteinA_Point6351.png" width=700px></p>
<p align="center"><img src="https://github.com/agds93/complexes_complementarity/blob/main/img/ProteinA_Point6351_processed.png" width=700px></p>
<p align="center"><i>Figure 2</i>: Original (top) and processed (bottom) average of a patch of surface A.</p>

The same things are represented in Figure 3 but referring to the patch with center `center_b`.

<p align="center"><img src="https://github.com/agds93/complexes_complementarity/blob/main/img/ProteinB_Point6026.png" width=700px></p>
<p align="center"><img src="https://github.com/agds93/complexes_complementarity/blob/main/img/ProteinB_Point6026_processed.png" width=700px></p>
<p align="center"><i>Figure 3</i>: Original (top) and processed (bottom) average of a patch of surface B.</p>

## Complementarietà
Given two patches in the contact zone (the first in the surface A and the second in the surface B) you want to know how similar they are to each other.  
For this purpose:
* Patches must have normal versors facing in opposite directions so that they are comparable.
* Patches must not have empty pixels (or islands) within the unit circle where the Zernike formalism is defined.

The first condition is satisfied if averages are generated via the `PatchesMethods` function, where the first surface patch is facing up and the other down.  
Instead, the second condition is met if the matrix of the mean returned by the function `ZernikeCoeff` is used. This new matrix has

