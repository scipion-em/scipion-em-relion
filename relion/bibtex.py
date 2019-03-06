# coding: latin-1
# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
"""

@article{Scheres2012a,
title = "A Bayesian View on Cryo-EM Structure Determination ",
journal = "JMB",
volume = "415",
number = "2",
pages = "406 - 418",
year = "2012",
issn = "0022-2836",
doi = "http://dx.doi.org/10.1016/j.jmb.2011.11.010",
url = "http://www.sciencedirect.com/science/article/pii/S0022283611012290",
author = "Scheres, Sjors H.W.",
keywords = "cryo-electron microscopy, three-dimensional reconstruction, maximum a posteriori estimation "
}

@article{Scheres2012b,
title = "RELION: Implementation of a Bayesian approach to cryo-EM structure determination ",
journal = "JSB",
volume = "180",
number = "3",
pages = "519 - 530",
year = "2012",
issn = "1047-8477",
doi = "http://dx.doi.org/10.1016/j.jsb.2012.09.006",
url = "http://www.sciencedirect.com/science/article/pii/S1047847712002481",
author = "Scheres, Sjors H.W.",
keywords = "Electron microscopy, Single-particle analysis, Maximum likelihood, Image processing, Software development "
}

@article{Kimanius2016,
title = "Accelerated cryo-EM structure determination with parallelisation using GPUs in RELION-2",
journal = "eLife",
volume = "5",
pages = "e18722",
year = "2016",
doi = "http://dx.doi.org/10.7554/eLife.18722",
author = "Kimanius, Dari and Forsberg, Björn O. and Scheres, Sjors H.W. and Lindahl, Erik",
}

@article{Scheres2015,
title = "Semi-automated selection of cryo-EM particles in RELION-1.3",
journal = "JSB",
volume = "189",
issue = "2",
pages = "114-122",
year = "2015",
doi = "http://dx.doi.org/10.1016/j.jsb.2014.11.010",
author = "Scheres, Sjors H.W.",
}

@article{Chen2012,
title = "Prevention of overfitting in cryo-EM structure determination",
journal = "Nat. Meth.",
volume = "9",
number = "3",
pages = "853 - 854",
year = "2012",
doi = "http://dx.doi.org/10.1038/nmeth.2115",
author = "Chen, Shaoxia and Scheres, Sjors H.W.",
}

@article{Chen2013,
title = "High-resolution noise substitution to measure overfitting and validate resolution in 3D structure determination by single particle electron cryomicroscopy.",
journal = "JSB",
volume = "135",
number = "3",
pages = "24-35",
year = "2013",
doi = "http://dx.doi.org/10.1016/j.ultramic.2013.06.004",
author = "Chen, Shaoxia and McMullan, Greg and Faruqi, Abdul R. and Murshudov, Garib N. and Short, Judith M. and Scheres, Sjors H.W. and Henderson, Richard",
}

@article{Zivanov2018,
title = "New tools for automated high-resolution cryo-EM structure determination in RELION-3",
journal = "eLife",
volume = "7",
pages = "e42166",
year = "2018",
doi = "http://dx.doi.org/10.7554/eLife.42166",
author = "Zivanov, Jasenko and Nakane, Takanori and Forsberg, Björn O and Kimanius, Dari and Hagen, Wim JH and Lindahl, Erik and Scheres, Sjors H.W.",
}

@article{Nakane2018,
title = "Characterisation of molecular motions in cryo-EM single-particle data by multi-body refinement in RELION",
journal = "eLife",
volume = "7",
pages = "e36861",
year = "2018",
doi = "http://dx.doi.org/10.7554/eLife.36861",
author = "Nakane, Takanori and Kimanius, Dari and Lindahl, Erik and Scheres, Sjors H.W.",
}

@article{Zivanov2019,
author = "Zivanov, Jasenko and Nakane, Takanori and Scheres, Sjors H. W.",
title = "{A Bayesian approach to beam-induced motion correction in cryo-EM single-particle analysis}",
journal = "IUCrJ",
year = "2019",
volume = "6",
number = "1",
pages = "5-17",
doi = {http://dx.doi.org/10.1107/S205225251801463X},
abstract = {A new method to estimate the trajectories of particle motion and the amount of cumulative beam damage in electron cryo-microscopy (cryo-EM) single-particle analysis is presented. The motion within the sample is modelled through the use of Gaussian process regression. This allows a prior likelihood that favours spatially and temporally smooth motion to be associated with each hypothetical set of particle trajectories without imposing hard constraints. This formulation enables the {\it a posteriori} likelihood of a set of particle trajectories to be expressed as a product of that prior likelihood and an observation likelihood given by the data, and this {\it a posteriori} likelihood to then be maximized. Since the smoothness prior requires three parameters that describe the statistics of the observed motion, an efficient stochastic method to estimate these parameters is also proposed. Finally, a practical algorithm is proposed that estimates the average amount of cumulative radiation damage as a function of radiation dose and spatial frequency, and then fits relative {\it B} factors to that damage in a robust way. The method is evaluated on three publicly available data sets, and its usefulness is illustrated by comparison with state-of-the-art methods and previously published results. The new method has been implemented as Bayesian polishing in {\it RELION}-3, where it replaces the existing particle-polishing method, as it outperforms the latter in all tests conducted.},
keywords = {Bayesian particle polishing, beam-induced motion correction, cryo-EM, single-particle analysis, electron cryo-microscopy},
}

"""
