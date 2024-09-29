=============
cisTEM plugin
=============

This plugin provides wrappers for several programs of `cisTEM <https://cistem.org>`_ software suite.

.. image:: https://img.shields.io/pypi/v/scipion-em-cistem.svg
        :target: https://pypi.python.org/pypi/scipion-em-cistem
        :alt: PyPI release

.. image:: https://img.shields.io/pypi/l/scipion-em-cistem.svg
        :target: https://pypi.python.org/pypi/scipion-em-cistem
        :alt: License

.. image:: https://img.shields.io/pypi/pyversions/scipion-em-cistem.svg
        :target: https://pypi.python.org/pypi/scipion-em-cistem
        :alt: Supported Python versions

.. image:: https://img.shields.io/sonar/quality_gate/scipion-em_scipion-em-cistem?server=https%3A%2F%2Fsonarcloud.io
        :target: https://sonarcloud.io/dashboard?id=scipion-em_scipion-em-cistem
        :alt: SonarCloud quality gate

.. image:: https://img.shields.io/pypi/dm/scipion-em-cistem
        :target: https://pypi.python.org/pypi/scipion-em-cistem
        :alt: Downloads

Installation
------------

You will need to use 3.0+ version of Scipion to be able to run these protocols. To install the plugin, you have two options:

a) Stable version

.. code-block::

    scipion installp -p scipion-em-cistem

b) Developer's version

    * download repository

    .. code-block::

        git clone -b devel https://github.com/scipion-em/scipion-em-cistem.git

    * install

    .. code-block::

        scipion installp -p /path/to/scipion-em-cistem --devel

**Important:** Starting from plugin v3.9, the config variables have been renamed. See: `scipion3 config -p cistem`

cisTEM and CTFFind binaries will be installed automatically with the plugin, but you can also link an existing installation.

    * Default cisTEM installation path assumed is ``software/em/cistem-1.0.0-beta``, if you want to change it, set *CISTEM_HOME* in ``scipion.conf`` file pointing to the folder where the cisTEM is installed.
    * Default CTFFind installation path assumed is ``software/em/ctffind-5.0.2``,if you want to change it, set *CTFFIND_HOME* in ``scipion.conf`` file pointing to the folder where CTFFind is installed.

A complete list of tests can be seen by executing ``scipion test --show --grep cistem``

Supported versions
------------------

* 1.0.0-beta (cistem)
* 4.1.14 (ctffind4)
* 5.0, 5.0.2 (ctffind5)

Licenses
--------

* cistem binaries are distributed under `The Janelia Research Campus Software License 1.2 license <http://license.janelia.org/license/janelia_license_1_2.html>`_
* ctffind4 binaries are distributed under `The Janelia Research Campus Software License 1.1 license <https://www.janelia.org/node/47808>`_
* ctffind5 binaries are distributed under `The Janelia Research Campus Software License 1.2 license <http://license.janelia.org/license/janelia_license_1_2.html>`_

Protocols
---------

* ctffind
* import tomo CTFs
* tilt-series ctffind
* unblur
* find particles
* classify 2D

References
----------

1. Johannes Elferich, Lingli Kong, Ximena Zottig, Nikolaus Grigorieff. (2024) CTFFIND5 provides improved insight into quality, tilt and thickness of TEM samples. bioRxiv 2024.02.26.582023; doi: https://doi.org/10.1101/2024.02.26.582023
2. Timothy Grant and Alexis Rohou and Nikolaus Grigorieff. (2018) cisTEM, user-friendly software for single-particle image processing. eLife 7:e35383.
