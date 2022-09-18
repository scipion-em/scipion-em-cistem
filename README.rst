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

cisTEM binaries will be installed automatically with the plugin, but you can also link an existing installation.

    * Default installation path assumed is ``software/em/cistem-1.0.0-beta``, if you want to change it, set *CISTEM_HOME* in ``scipion.conf`` file pointing to the folder where the cisTEM is installed.
    * It's possible to use CTFFIND4 installed separately from cisTEM by defining *CTFFIND4_HOME* variable in ``scipion.conf``. Setting this value will have the priority over Ctffind inside cisTEM.

A complete list of tests can be seen by executing ``scipion test --show --grep cistem``

Supported versions
------------------

* 1.0.0-beta (cistem)
* 4.1.13, 4.1.14 (ctffind4)

Protocols
---------

* ctffind4
* import tomo CTFs
* tilt-series ctffind4
* unblur
* find particles
* classify 2D

References
----------

1. Timothy Grant and Alexis Rohou and Nikolaus Grigorieff. (2018) cisTEM, user-friendly software for single-particle image processing. eLife 7:e35383.
