=============
Relion plugin
=============

This plugin provides wrappers for several programs of `RELION <https://relion.readthedocs.io/en/latest/index.html>`_ software suite.

.. image:: https://img.shields.io/pypi/v/scipion-em-relion.svg
        :target: https://pypi.python.org/pypi/scipion-em-relion
        :alt: PyPI release

.. image:: https://img.shields.io/pypi/l/scipion-em-relion.svg
        :target: https://pypi.python.org/pypi/scipion-em-relion
        :alt: License

.. image:: https://img.shields.io/pypi/pyversions/scipion-em-relion.svg
        :target: https://pypi.python.org/pypi/scipion-em-relion
        :alt: Supported Python versions

.. image:: https://img.shields.io/sonar/quality_gate/scipion-em_scipion-em-relion?server=https%3A%2F%2Fsonarcloud.io
        :target: https://sonarcloud.io/dashboard?id=scipion-em_scipion-em-relion
        :alt: SonarCloud quality gate

.. image:: https://img.shields.io/pypi/dm/scipion-em-relion
        :target: https://pypi.python.org/pypi/scipion-em-relion
        :alt: Downloads


**IMPORTANT NOTES!**

    1. If you have imported movies with a gain file in **DM4** format, you need to **flip the gain reference upside-down** in the motion correction protocol! (`bug details <https://github.com/I2PC/xmippCore/issues/39>`_)
    2. If you have provided a gain reference or defects file during movie import or motion correction and plan to run bayesian polishing, please **make sure to run first "assign optics groups" protocol for aligned movies**, specifying the gain file etc. Currently, Scipion has no other way of knowing if you have e.g. rotated the gain during motion correction. Output movies then can be used in the polishing protocol.
    3. When importing EER movies, you should specify dose per single EER frame during import step. Pixel size should be the physical pixel size unless you plan to render EER on 8K grid. In the Motioncor protocol choose binning, EER fractionation and upsampling (default=1 is for 4K grid).

Installation
------------

You will need to use 3.0+ version of Scipion to be able to run these protocols. To install the plugin, you have two options:

a) Stable version

   .. code-block::

      scipion installp -p scipion-em-relion

b) Developer's version

   * download repository

   .. code-block::

      git clone -b devel https://github.com/scipion-em/scipion-em-relion.git

   * install

   .. code-block::

      scipion installp -p /path/to/scipion-em-relion --devel

- RELION sources will be downloaded and compiled automatically with the plugin, but you can also link an existing installation. Default installation path assumed is ``software/em/relion-5.0``, if you want to change it, set *RELION_HOME* in ``scipion.conf`` file to the folder where the RELION is installed.
- If you need to use CUDA different from the one used during Scipion installation (defined by *CUDA_LIB*), you can add *RELION_CUDA_LIB* variable to the config file. Optionally, you can also specify *RELION_CUDA_BIN* path for nvcc.
- If you have to use a MPI for Relion different from Scipion MPI, you can set *RELION_MPI_BIN* and *RELION_MPI_LIB* variables in the config file.
- To add support for Python modules (e.g. Blush, ModelAngelo and DynaMight) you will have to setup a Python environment with dependencies. You need to set *RELION_ENV_ACTIVATION* that points to a conda environment. Default = ``conda activate relion-5.0``
- You might want to set *TORCH_HOME* pointing to the path with downloaded models. Default = ``software/em/modelangelomodels-1.0`` (shared with ModelAngelo plugin).
- If you want to use SIDESPLITTER, you need the `sidesplitter <https://github.com/scipion-em/scipion-em-sidesplitter>`_ plugin installed and *SIDESPLITTER_HOME* set properly. After that, you only need to add *--external_reconstruct* into Additional arguments field of a 3D auto-refine or multi-body job.

To check the installation, simply run one of the tests. A complete list of tests can be displayed by executing ``scipion test --show --grep relion``

Supported versions
------------------

4.0, 5.0

Protocols
---------

* 2D class ranker
* 2D classification         
* 3D auto-refine            
* 3D classification         
* 3D initial model          
* 3D multi-body
* assign optics groups
* auto-picking (reference-based)
* auto-picking LoG          
* bayesian polishing
* calculate fsc
* center averages
* clean project
* compress movies
* create 3d mask
* crop / resize volumes
* ctf refinement
* DynaMight flexibility
* estimate gain to compress
* expand symmetry
* export coordinates
* export ctf                
* export particles
* local resolution
* motion correction
* particles extraction
* post-processing           
* preprocess particles      
* reconstruct
* remove preferential views
* subtract projection
* symmetrize volume

References
----------

1. Scheres et al., JMB, 2012 
2. Scheres et al., JSB, 2012 
3. Kimanius et al., eLife, 2016 
4. Zivanov et al., eLife, 2018
5. Kimanius et al., Biochemical Journal, 2021
