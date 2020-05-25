# coding: latin-1
# **************************************************************************
# *
# *  Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology (MRC-LMB)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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

@article{Grant2018,
title = "cisTEM, user-friendly software for single-particle image processing",
journal = "eLife",
volume = "7",
pages = "e35383",
year = "2018",
doi = "http://doi.org/10.7554/eLife.35383.001",
author = "Grant, T. and Rohou, A. and Grigorieff, N.",
}

@article{Campbell2012,
title = "Movies of ice-embedded particles enhance resolution in electron cryo-microscopy",
journal = "Structure",
volume = "20",
year = "2012",
month = "9/2012",
pages = "1823 - 1828",
issn = "09692126",
doi = "http://dx.doi.org/10.1016/j.str.2012.08.026",
author = "Campbell, M. G. and Cheng, A. and Brilot, A. F. and Moeller, A. and Lyumkis, D. and Veesler, D. and Pan, J. and Harrison, S. C. and Potter, C. S. and Carragher, B. and Grigorieff, N.",
}

@article{Grant2015b,
title = "Measuring the optimal exposure for single particle cryo-EM using a 2.6A reconstruction of rotavirus VP6",
journal = "eLife",
volume = "4",
year = "2015",
month = "06/2015",
pages = "1-19",
doi = "http://dx.doi.org/10.7554/eLife.06980",
url = "https://elifesciences.org/articles/06980",
author = "Grant, T. and Grigorieff, N.",
}

@article{Mindell2003,
title = "Accurate determination of local defocus and specimen tilt in electron microscopy ",
journal = "JSB ",
volume = "142",
number = "3",
pages = "334 - 347",
year = "2003",
note = "",
issn = "1047-8477",
doi = "http://dx.doi.org/10.1016/S1047-8477(03)00069-8",
url = "http://www.sciencedirect.com/science/article/pii/S1047847703000698",
author = "Mindell, Joseph A. and Grigorieff, Nikolaus",
keywords = "Electron microscopy, Contrast transfer function, Algorithm, Tilt determination ",
}

@Article{Rohou2015,
author="Rohou, A.  and Grigorieff, N. ",
title="{{C}{T}{F}{F}{I}{N}{D}4: {F}ast and accurate defocus estimation from electron micrographs}",
journal="J. Struct. Biol.",
year="2015",
volume="192",
number="2",
pages="216 - 221",
month="Nov",
doi = "http://dx.doi.org/10.1016/j.jsb.2015.08.008",
url = "http://www.sciencedirect.com/science/article/pii/S1047847715300460",
}

@article{Sigworth2004,
title = "Classical detection theory and the cryo-EM particle selection problem",
journal = "Journal of Structural Biology",
volume = "145",
number = "1",
pages = "111 - 122",
year = "2004",
issn = "1047-8477",
doi = "http://doi.org/10.1016/j.jsb.2003.10.025",
url = "http://www.sciencedirect.com/science/article/pii/S1047847703002405",
author = "Fred J. Sigworth",
}

@article{Sigworth1998,
title = "A Maximum-Likelihood Approach to Single-Particle Image Refinement",
journal = "Journal of Structural Biology",
volume = "122",
number = "3",
pages = "328 - 339",
year = "1998",
issn = "1047-8477",
doi = "http://doi.org/10.1006/jsbi.1998.4014",
url = "http://www.sciencedirect.com/science/article/pii/S104784779894014X",
author = "F.J. Sigworth",
keywords = "electron microscopy, maximum likelihood, single-particle alignment.",
}

@article{Scheres2005,
title = "Maximum-likelihood Multi-reference Refinement for Electron Microscopy Images",
journal = "Journal of Molecular Biology",
volume = "348",
number = "1",
pages = "139 - 149",
year = "2005",
issn = "0022-2836",
doi = "http://doi.org/10.1016/j.jmb.2005.02.031",
url = "http://www.sciencedirect.com/science/article/pii/S0022283605001932",
author = "Sjors H.W. Scheres and Mikel Valle and Rafael Nunez and Carlos O.S. Sorzano and Roberto Marabini and Gabor T. Herman and Jose-Maria Carazo",
keywords = "maximum-likelihood, multi-reference refinement, single-particles, 2D-alignment, classification",
}

@article{Liang2015,
doi = "10.1016/j.cell.2015.06.018",
url = "http://doi.org/10.1016/j.cell.2015.06.018",
year = "2015",
month = "jul",
volume = "162",
number = "2",
pages = "314-327",
author = "Bo Liang and Zongli Li and Simon Jenni and Amal A. Rahmeh and Benjamin M. Morin and Timothy Grant and Nikolaus Grigorieff and Stephen C. Harrison and Sean P.J. Whelan",
title = "Structure of the L Protein of Vesicular Stomatitis Virus from Electron Cryomicroscopy",
journal = "Cell"
}

@article{Grigorieff2016,
title = "Chapter Eight - Frealign: An Exploratory Tool for Single-Particle Cryo-EM",
journal = "Methods in Enzymology",
volume = "579",
pages = "191 - 226",
year = "2016",
booktitle = "The Resolution Revolution: Recent Advances In cryoEM",
issn = "0076-6879",
doi = "http://doi.org/10.1016/bs.mie.2016.04.013",
url = "http://www.sciencedirect.com/science/article/pii/S0076687916300313",
author = "N. Grigorieff",
}

"""
