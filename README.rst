=============
Relion plugin
=============

This plugin provide wrappers around several programs of `RELION <https://www3.mrc-lmb.cam.ac.uk/relion/index.php/Main_Page>`_ software suite.

+------------------+------------------+
| stable: |stable| | devel: | |devel| |
+------------------+------------------+

.. |stable| image:: http://scipion-test.cnb.csic.es:9980/badges/relion_prod.svg
.. |devel| image:: http://scipion-test.cnb.csic.es:9980/badges/relion_sdevel.svg


Installation
------------

You will need to use `3.0 <https://github.com/I2PC/scipion/releases/tag/V3.0.0>`_ version of Scipion to be able to run these protocols. To install the plugin, you have two options:

a) Stable version

   .. code-block::

      scipion installp -p scipion-em-relion

b) Developer's version

   * download repository

   .. code-block::

      git clone https://github.com/scipion-em/scipion-em-relion.git

   * install

   .. code-block::

      scipion installp -p path_to_scipion-em-relion --devel

RELION sources will be downloaded and compiled automatically with the plugin, but you can also link an existing installation. Default installation path assumed is ``software/em/relion-3.1``, if you want to change it, set *RELION_HOME* in ``scipion.conf`` file to the folder where the RELION is installed. If you need to use CUDA different from the one used during Scipion installation (defined by *CUDA_LIB*), you can add *RELION_CUDA_LIB* variable to the config file. Moreover, if you have to use a MPI for Relion different from Scipion MPI, you can set *RELION_MPI_BIN* and *RELION_MPI_LIB* variables in your shell environment - they will be recognized by Scipion.

To check the installation, simply run one of the tests. A complete list of tests can be displayed by executing ``scipion test --show --grep relion``

Supported versions
------------------

3.0, 3.1

In 2020 the plugin was updated to support the latest RELION 3.1. We discontinued any support for 2.x versions. We are still working towards stable plugin release, so some things might not work yet.

Protocols
---------

* 2D classification         
* 3D auto-refine            
* 3D classification         
* 3D initial model          
* 3D multi-body
* assign optics group
* auto-picking              
* auto-picking LoG          
* bayesian polishing        
* center averages           
* create 3d mask            
* ctf refinement            
* expand symmetry           
* export ctf                
* export particles          
* local resolution          
* motion correction
* particles extraction
* post-processing           
* preprocess particles      
* reconstruct               
* subtract projection
* symmetrize volume

References
----------

1. Scheres et al., JMB, 2012 
2. Scheres et al., JSB, 2012 
3. Kimanius et al., eLife, 2016 
4. Zivanov et al., eLife, 2018
