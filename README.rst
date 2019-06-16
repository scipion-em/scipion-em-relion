=============
Relion plugin
=============

This plugin provide wrappers around several programs of `RELION <https://www3.mrc-lmb.cam.ac.uk/relion/index.php/Main_Page>`_ software suite.

.. figure:: http://scipion-test.cnb.csic.es:9980/badges/relion_devel.svg
   :align: left
   :alt: build status

Installation
------------

You will need to use `2.0 <https://github.com/I2PC/scipion/releases/tag/V2.0.0>`_ version of Scipion to be able to run these protocols. To install the plugin, you have two options:

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

RELION sources will be downloaded and compiled automatically with the plugin, but you can also link an existing installation. Default installation path assumed is ``software/em/relion-3.0``, if you want to change it, set *RELION_HOME* in ``scipion.conf`` file to the folder where the RELION is installed. If you need to use CUDA different from the one used during Scipion installation (defined by *CUDA_LIB*), you can add *RELION_CUDA_LIB* variable to the config file. Moreover, if you have to use a MPI for Relion different from Scipion MPI, you can set *RELION_MPI_BIN* and *RELION_MPI_LIB* variables in your shell environment - they will be recognized by Scipion.

To check the installation, simply run one of the following Scipion tests:

.. code-block::

   scipion test pyworkflow.tests.em.workflows.test_workflow_streaming.TestRelionPickStreaming
   scipion test pyworkflow.tests.em.workflows.test_workflow_streaming.TestRelionExtractStreaming
   scipion test pyworkflow.tests.em.workflows.test_workflow_mixed_large.TestMixedRelionTutorial
   scipion test relion.tests.test_workflow_relion3.TestWorkflowRelion3Betagal
   scipion test relion.tests.test_convert_relion.TestReconstruct
   scipion test relion.tests.test_convert_relion.TestConvertBinaryFiles
   scipion test relion.tests.test_convert_relion.TestConversions
   scipion test relion.tests.test_convert_relion.TestAlignment
   scipion test relion.tests.test_protocols_relion.TestRelionSubtract
   scipion test relion.tests.test_protocols_relion.TestRelionSortParticles
   scipion test relion.tests.test_protocols_relion.TestRelionRefine
   scipion test relion.tests.test_protocols_relion.TestRelionPreprocess
   scipion test relion.tests.test_protocols_relion.TestRelionPostprocess
   scipion test relion.tests.test_protocols_relion.TestRelionLocalRes
   scipion test relion.tests.test_protocols_relion.TestRelionInitialModel
   scipion test relion.tests.test_protocols_relion.TestRelionExtractParticles
   scipion test relion.tests.test_protocols_relion.TestRelionExtractMovieParticles
   scipion test relion.tests.test_protocols_relion.TestRelionExportParticles
   scipion test relion.tests.test_protocols_relion.TestRelionExpandSymmetry
   scipion test relion.tests.test_protocols_relion.TestRelionCreate3dMask
   scipion test relion.tests.test_protocols_relion.TestRelionClassify3D
   scipion test relion.tests.test_protocols_relion.TestRelionClassify2D
   scipion test relion.tests.test_protocols_relion.TestRelionCenterAverages
   scipion test relion.tests.test_workflow_relion.TestWorkflowRelionPick
   scipion test relion.tests.test_workflow_relion.TestWorkflowRelionExtract

A complete list of tests can also be seen by executing ``scipion test --show --grep relion``

Supported versions
------------------

2.0.4, 2.1, 3.0

In 2018 the plugin was updated to support the latest (at that moment) RELION: 3.0. This required a lot of code refactoring and the support of old RELION versions (1.x) had to be discontinued. Many new protocols specific to new RELION release we added. You are welcome to check them out and give us your feedback. We are still fixing few issues related to RELION 3 support, please do let us know if something does not work as expected.

Protocols
---------

* 2D classification         
* 3D auto-refine            
* 3D classification         
* 3D initial model          
* 3D multi-body             
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
* motioncor                 
* movie particles extraction
* particle polishing        
* particles extraction      
* post-processing           
* preprocess particles      
* reconstruct               
* sort particles            
* subtract projection       
* symmetrize volume

References
----------

1. Scheres et al., JMB, 2012 
2. Scheres et al., JSB, 2012 
3. Kimanius et al., eLife, 2016 
4. Zivanov et al., eLife, 2018
