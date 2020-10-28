# Bifurcations_Savanna
Code for PNAS paper "On the complex dynamics of savanna landscapes" by Touboul, Staver and Levin

General Remarks:
  - Codimension 2 bifurcation diagrams generated with Matcont (https://sourceforge.net/projects/matcont/). Code is provided in the folder "Matcont Code Codim 2". There, two files define the system: System.m and System.mat; the curves computed for each figure are saved in a .mat file labeled FigureNumberAndPanel.mat.
  - Codimension 1 bifurcation diagrams and trajectories are generated using XPP Aut (http://www.math.pitt.edu/~bard/xpp/xpp.html), and all files associated are in the folder "XPP Files". A system file .ode is defined, and set files are saved for each panel with precise parameter sets used in the simulations for each panel set. These set files are labeled System_FigureNumber_PanelLabel.set and can be loaded directly in XPP Aut. 
  - Stochastic simulations were generated using a custom code based on Euler-Maruyama scheme. Files are stored in the folder Matlab - Stochastic Simulations. 
  - Some curves were obtained using formal calculations, particularly Figs. S1, S2, S3. Code is in the Maple folder.
  
Code associated with each figure:

Figure 1: Matcont system "Trees2d". 

Figure 2, 3, S5, S6, S7 : Codimension 2 diagram (e.g., Fig.2A) generated with Matcont System "Trees". Codimension 1 bifurcation diagrams (Fig.2 panels B and C and Fig. 3) were generated with XPP system Trees.ode, and set files associated are stored as indicated in the general nomenclature above (e.g., Trees.Fig2_C1_3.set for the third plot of panel C1 in Figure 2). 

Figure 4 and S8: Codimension 2 Matcont file Trees_Modified.m, codim 1 bifurcation diagrams and time traces generated with XPP Aut Trees_Modified.ode (set files Trees_Modified.Fig4_*panel label*.set)

Figure 5, 6, S9, S10, S11 were generated with a custom code available in the folder Matlab-Stochastic Simulations-Fig5_6_S11. The figure label is indicated in the name of the file. Figure S9 and S10 is generate with the same code as Fig.5.
