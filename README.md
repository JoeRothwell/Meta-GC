# The Meta-GC project : untargeted metabolomics and gastric cancer
This repo contains the code for data analysis for the Meta-GC project.

---

Meta-GC is a case-control study nested in the EPIC prospective cohort. It includes 437 GC cases matched to 437 cancer-free controls. Plasma samples from participants were analysed by untargeted metabolomics using a Thermo Q-Exactive mass spectrometer. After alignment and filtering, over 2000 spectral features were found in positive and negative ionisation modes.

---

This code performs multivariable conditional logistic regressions on each feature in turn, using several different models. The resulting p-values are adjusted by controlling the false discovery rate (FDR) and those significant at FDR p < 0.05 retained. The data are also visualised on smile plots of OR vs p-value.

---

Different subgroups of data are also tested, and H. pylori status taken into account.

