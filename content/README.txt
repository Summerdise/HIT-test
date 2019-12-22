The source code and data in this folder reproduces the learning experiments of hybrid image templates.



- To learn and display hybrid image templates, use "run_demo.m" in matlab. 
  The display will be paused after each template is learned; press any key to continue.
  The hybrid image templates are stored in the root folder as .png images.

- To reproduce Figure 3 of the paper, use "display_hedgehog.m", "display_pighead.m", "display_pigeon.m".
  The result figure is saved as "matched_hedgehog.png" and similar names for the other two categories.

- To reproduce the contrast templates in Figure 16, use "illustrateContrastTemplate.m".
  The result figures are saved as "constrast_basis_<category1>_<category2>_<pos or neg>.png" for 
  contrast templates obtained by active basis model; 
  and "constrast_svm_<category1>_<category2>_<pos or neg>.png" for contrast templates obtained by SVM.


