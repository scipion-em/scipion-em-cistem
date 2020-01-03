=============
cisTEM plugin
=============

**ATTENTION: This plugin is still in development and will eventually replace grigoriefflab plugin.**

This plugin provide wrappers around several programs of `cisTEM <https://cistem.org>`_ software suite.

.. figure:: http://scipion-test.cnb.csic.es:9980/badges/cistem_devel.svg
   :align: left
   :alt: build status

Installation
------------

You will need to use `2.0 <https://github.com/I2PC/scipion/releases/tag/V2.0.0>`_ version of Scipion to be able to run these protocols. To install the plugin, you have two options:

a) Stable version

.. code-block::

    scipion installp -p scipion-em-cistem

b) Developer's version

    * download repository

    .. code-block::

        git clone https://github.com/scipion-em/scipion-em-cistem.git

    * install

    .. code-block::

        scipion installp -p path_to_scipion-em-cistem --devel

cisTEM binaries will be installed automatically with the plugin, but you can also link an existing installation.

    * Default installation path assumed is ``software/em/cistem-1.0.0-beta``, if you want to change it, set *CISTEM_HOME* in ``scipion.conf`` file pointing to the folder where the cisTEM is installed.
    * It's possible to use CTFFIND4 installed separately from cisTEM by defining *CTFFIND4_HOME* variable in ``scipion.conf``.

To check the installation, simply run one of the following Scipion tests:

.. code-block::

   scipion test cistem.tests.test_protocols_cistem_movies.TestUnblur
   scipion test cistem.tests.test_protocols_cistem.TestCtffind4
   scipion test cistem.tests.test_protocols_cistem.TestFindParticles
   scipion test cistem.tests.test_protocols_cistem.TestRefine2D

A complete list of tests can also be seen by executing ``scipion test --show --grep cistem``

Supported versions
------------------

1.0.0-beta


Protocols
---------

* ctffind4
* unblur
* find particles
* classify 2D

References
----------

1. Timothy Grant and Alexis Rohou and Nikolaus Grigorieff. (2018) cisTEM, user-friendly software for single-particle image processing. eLife 7:e35383.
