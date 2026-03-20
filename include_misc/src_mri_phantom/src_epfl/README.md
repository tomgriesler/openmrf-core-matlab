// README for the code in this folder
//
// The package is available online at http://bigwww.epfl.ch/algorithms/mriphantom/
//
// This source code is related to the paper:
//	Realistic Analytical Phantoms for Parallel Magnetic Resonance Imaging
//	by M. Guerquin-Kern, L. Lejeune, K. P. Pruessmann, and M. Unser,
//	accepted for publication in IEEE Transactions on Medical Imaging.
//
// It extends to code that was provided with the paper:
//	ANALYTICAL FORM OF SHEPP-LOGAN PHANTOM FOR PARALLEL MRI
//	from M. Guerquin-Kern, F. I. Karahanoglu, D. Van De Ville, K. P. Pruessmann and M. Unser,
//      Proceedings of the Seventh IEEE International Symposium on Biomedical Imaging: From Nano to Macro (ISBI'10),
//      Rotterdam, The Netherlands, April 14-17, 2010, pp. 261-264.
//
// CONDITIONS OF USE:
// You are free to use this software for research purposes, but you should not redistribute it without our consent.
// It would also be courteous for you to include the aforementioned citation and acknowledgment whenever you present
// or publish results that are based on it. EPFL makes no warranties of any kind on this software and shall in no
// event be liable for damages of any kind in connection with the use and exploitation of this technology.
//
// Matthieu Guerquin-Kern, Biomedical Imaging Group - EPFL, 2011-05-04 (yyyy-mm-dd)
//
// Version 0.8

HISTORY:
* Versions 0.1 to 0.3: 	initial code with code efficiency improvements (cf. ISBI paper, 2010)
* Version 0.4: 		added implementation for the sinusoidal sensitivity model, the experiments of the ISBI paper are reproduced with that model
* Version 0.5:		added GUI to generate phantom composed of ellipses, polygons and quadratic Bezier curves, the rasterization of such phantoms is supported.
* Version 0.6:      	added SVG and PDF export for the phantom
* Version 0.7:      	added support for polygon and Bezier curves is yet to be added (cf. Laurent Lejeune's project).
* Version 0.8:		added experiments

FILES:
To start with the code, you can run and look into 'demo.m'
The files 'demo_SL.m', 'isbi_exp2.m' and 'isbi_exp3.m' shall reproduce the results presented in the aforementionned ISBI paper.
New files that add experiments with the sinusoidal model and the new class of phantoms:
	* 'exp_sens_models.m'		comparison of the two sensitivity models
	* 'exp_Gibbs.m'			illustration of the aliasing and Gibbs phenomenon depending on the simulation method
	* 'exp_rectangle.m'		validation of the theory and implementation in the case of homogeneous sensitivity
	* 'exp_validation_sens.m'	validation of the implementation with non-homogeneous sensitivity
	* 'exp_recons'			SENSE reconstructions out of data simulated with different method, illustration of
					the impact of rasterized simulations on reconstruction performances.

INSTALLATION:
Add the folder containing this code in your matlab path.
Some C++ sources are distributed with the code. To compile them run 'make' in matlab.

DEPENDENCIES:
This code was tested on Matlab R2011a.
Some functions require the IRT toolbox of Jeff Fessler for NUFT routines (cf. http://www.eecs.umich.edu/~fessler/irt/fessler.tgz). Install this toolbox in your matlab path.

ACKNOWLEGMENTS
* The code for error function of a complex variable is a homemade MEX/C++ implementation of Marcel Leutenegger's erfz.m ( (c) January 2008, LGPL version 2.1)
(cf. https://documents.epfl.ch/users/l/le/leuteneg/www/MATLABToolbox/ErrorFunction.html).
* For multinomial sums we rely on Matt Fig's npermutek.m implementation.
(cf. http://www.mathworks.co.kr/matlabcentral/fileexchange/11462-npermutek)
* Rasterization and design of phantoms rely on Bruno Luong's insidepoly MEX implementation rather than matlab's internal inpolygon
(cf http://www.mathworks.com/matlabcentral/fileexchange/27840-2d-polygon-interior-detection)

CONTACT:
For any information or comment, please contact Matthieu:
matthieu.guerquin-kern AT epfl DOT ch
