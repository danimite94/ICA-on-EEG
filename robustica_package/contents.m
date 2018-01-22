% RobustICA algorithm for independent component analysis (Release 1 - March 31, 2008)
% -----------------------------------------------------------------------------------
% 
% RobustICA is based on the normalized kurtosis contrast function, which is optimized by a
% computationally efficient gradient-descent technique. This technique computes algebraically
% the step size (adaption coefficient) globally optimizing the contrast in the search direction
% at each iteration. Any independent component with non-zero kurtosis can be extracted in this
% manner.
%
% The present implementation performs the deflationary separation of statistically independent
% sources under the instantaneous linear mixture model. Full separation is achieved if at most
% one source has zero kurtosis.
%
% Some advantages of RobustICA are:
%
% - The optimal step-size technique provides some robustness to the presence of saddle points and
%   spurious local extrema in the contrast function.
%
% - The method shows a very high convergence speed measured in terms of source extraction quality
%   versus number of operations.
%
% - Real- and complex-valued signals are treated by exactly the same algorithm. Both type of source
%   signals can be present simultaneously in a given mixture. Complex sources need not be circular.
%   The mixing matrix coefficients may be real or complex, regardless of the source type.
%
% - Sequential extraction (deflation) can be performed via linear regression. As a result, prewhitening
%   and the performance limitations it imposes can be avoided. This feature may prove especially
%   beneficial in ill-conditioned scenarios, the convolutive case and underdetermined mixtures.
%
% - Optionally, the algorithm can target sub-Gaussian or super-Gaussian sources in the order defined
%   by a kurtosis-sign vector provided by the user. 
%
%
% The package is composed of the following M-files:
%
%  - '<a href = "matlab:doc robustica">robustica.m</a>':             implements the algorithm itself.
%
%  - '<a href = "matlab:doc kurt_gradient_optstep">kurt_gradient_optstep.m</a>': computes the optimal step-size of the normalized kurtosis contrast
%                               using the gradient vector as search direction. 
%
%  - '<a href = "matlab:doc deflation_regression">deflation_regression.m</a>':  performs deflation via linear regression.
%
%  - '<a href = "matlab:doc robustica_demo">robustica_demo.m</a>':        a simple demonstration illustrating the performance of RobustICA
%                               on synthetic mixtures.
% 
%
% More details about the RobustICA algorithm can be found in the references below:
%
% - V. Zarzoso and P. Comon, <a href = "http://www.i3s.unice.fr/~zarzoso/biblio/ica-2007.pdf.gz">"Comparative Speed Analysis of FastICA"</a>, 
%   in: Proceedings ICA-2007, 7th International Conference on Independent Component Analysis
%   and Signal Separation, London, UK, September 9-12, 2007, pp. 293-300.
% 
% - V. Zarzoso, P. Comon and M. Kallel, <a href = "http://www.i3s.unice.fr/~zarzoso/biblio/eusipco-2006.pdf.gz">"How Fast is FastICA?"</a>, 
%   in: Proceedings EUSIPCO-2006, XIV European Signal Processing Conference, 
%   Florence, Italy, September 4-8, 2006. 
%
%
% Conditions of use:
%
%    This is open-source code. Please, feel free to edit the M-files and modify them as you wish.
%    It is only requested that original authorship be acknowledged and modifications be clearly
%    indicated in the code (history section and elsewhere as appropriate).
%      
%
% Compatibility:
%
%   The package has been tested on Matlab 6.5, 7.3 and 7.5. Hyperlinks in the help section of some
%   files do not work properly under version 6.5.
%
%
% Please, report any bugs, comments or suggestions to <a href = "mailto:zarzoso@i3s.unice.fr">zarzoso(a)i3s.unice.fr</a>.
%
%
% (c) 2008 V. Zarzoso, P. Comon.

