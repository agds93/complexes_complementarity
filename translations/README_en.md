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
