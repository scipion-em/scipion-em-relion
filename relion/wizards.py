# **************************************************************************
# *
# * Authors:     Jose Gutierrez (jose.gutierrez@cnb.csic.es) [1]
# *              J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [2]
# *
# * [1] Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# * [2] SciLifeLab, Stockholm University
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


import os
from collections import OrderedDict

import pyworkflow as pw
from pyworkflow.em import *
from pyworkflow.em.viewers import CoordinatesObjectView
from pyworkflow.em.wizard import *
import pyworkflow.em.metadata as md

from relion.constants import *
from relion.convert import writeSetOfMicrographs
from relion.protocols import (
    ProtRelionClassify3D, ProtRelionRefine3D, ProtRelionClassify2D,
    ProtRelionPreprocessParticles, ProtRelionAutopickLoG,
    ProtRelion2Autopick, ProtRelionCreateMask3D,
    ProtRelionSortParticles, ProtRelionInitialModel)

#===============================================================================
# MASKS
#===============================================================================

class RelionBackRadiusWizard(ParticleMaskRadiusWizard):
    _targets = [(ProtRelionPreprocessParticles, ['backRadius'])]
    _unit = UNIT_PIXEL
    
    def _getProtocolImages(self, protocol):
        return protocol.inputParticles
    
    def _getParameters(self, protocol):
        
        label, value = self._getInputProtocol(self._targets, protocol)
        
        protParams = {}
        protParams['input']= self._getProtocolImages(protocol)
        protParams['label']= label
        protParams['value']= value
        return protParams  
    
    def _getProvider(self, protocol):
        _objs = self._getParameters(protocol)['input']    
        return ParticleMaskRadiusWizard._getListProvider(self, _objs)
    
    def show(self, form):
        params = self._getParameters(form.protocol)
        _value = params['value']
        _label = params['label']
        ParticleMaskRadiusWizard.show(self, form, _value, _label, units=self._unit)


class RelionPartMaskDiameterWizard(RelionBackRadiusWizard):
    _targets = [(ProtRelionClassify2D, ['maskDiameterA']),
                (ProtRelionRefine3D, ['maskDiameterA']),
                (ProtRelionClassify3D, ['maskDiameterA']),
                (ProtRelionClassify2D, ['maskDiameterA']),
                (ProtRelionInitialModel, ['maskDiameterA'])]
    _unit = UNIT_ANGSTROM

    def _getParameters(self, protocol):
        protParams = RelionBackRadiusWizard._getParameters(self, protocol)
        # adjust to from diameter to radius
        protParams['value'] = protParams['value'] / 2

        return protParams

    def setVar(self, form, label, value):
        # adjust again from radius to diameter
        form.setVar(label, value * 2)

# We need this specific wizard for the sort protocol because
# this protocol have a particular way to grab the input images
class RelionSortMaskWizard(RelionPartMaskDiameterWizard):
    _targets = [(ProtRelionSortParticles, ['maskDiameterA'])]

    def _getProvider(self, protocol):
        if protocol.isInputClasses():
            images = [cls.clone()
                      for cls in protocol.inputSet.get().iterRepresentatives()]
        else:
            images = self._getParticles(protocol._allParticles(iterate=True))
        return ListTreeProvider(images)

    def _getProtocolImages(self, protocol):
        return None
        #return protocol._allParticles(iterate=True)


#===============================================================================
# FILTER
#===============================================================================

class RelionVolFilterWizard(FilterVolumesWizard):
    _targets = [(ProtRelionClassify3D, ['initialLowPassFilterA']),
                (ProtRelionRefine3D, ['initialLowPassFilterA']),
                (ProtRelionCreateMask3D, ['initialLowPassFilterA'])]
    
    def _getParameters(self, protocol):
        
        label, value = self._getInputProtocol(self._targets, protocol)
        protParams = {}

        if protocol.__class__.__name__ == 'ProtRelionCreateMask3D':
            protParams['input'] = protocol.inputVolume
        else:
            protParams['input'] = protocol.referenceVolume

        protParams['label'] = label
        protParams['value'] = value
        protParams['mode'] = FILTER_LOW_PASS_NO_DECAY
        return protParams  
    
    def _getProvider(self, protocol):
        _objs = self._getParameters(protocol)['input']    
        return FilterVolumesWizard._getListProvider(self, _objs)
    
    # def show(self, form):
    #     params = self._getParameters(form.protocol)
    #     # Value should be LowFreq=0, HighFreq and Decay for Low pass filter
    #     _value = params['value']
    #     _label = params['label']
    #     FilterVolumesWizard.show(self, form, _value, _label,
    #                              mode=FILTER_LOW_PASS,
    #                              unit=UNIT_ANGSTROM,
    #                              showDecay=False)

    def show(self, form):
        params = self._getParameters(form.protocol)
        protocol = form.protocol
        provider = self._getProvider(protocol)

        if provider is not None:
            args = {'mode': params['mode'],
                    'highFreq': params['value'],
                    'unit': UNIT_ANGSTROM
                    }

            args['showLowFreq'] = False
            args['showDecay'] = False

            d = BandPassFilterDialog(form.root, provider, **args)

            if d.resultYes():
                form.setVar('initialLowPassFilterA', d.samplingRate/d.getHighFreq())

        else:
            dialog.showWarning("Input volumes", "Select volumes first", form.root)
            

#===============================================================================
# PICKING
#===============================================================================

# class RelionPartDiameter(RelionPartMaskDiameterWizard):
#     _targets = [(ProtRelionAutopickFom, ['particleDiameter'])]
#
#     def _getProtocolImages(self, protocol):
#         return protocol.inputReferences


class Relion2AutopickParams(EmWizard):
    _targets = [(ProtRelion2Autopick, ['runType',
                                       'pickingThreshold',
                                       'interParticleDistance'])]

    def show(self, form):
        autopickProt = form.protocol

        if not autopickProt.hasAttribute('outputCoordinatesSubset'):
            form.showWarning("You should run the procotol in 'Optimize' mode "
                               "at least once before opening the wizard.")
            return

        project = autopickProt.getProject()
        micSet = autopickProt.outputMicrographsSubset
        micfn = micSet.getFileName()
        coordsDir = project.getTmpPath(micSet.getName())
        pw.utils.cleanPath(coordsDir)
        pw.utils.makePath(coordsDir)

        micStarFn = os.path.join(coordsDir, 'micrographs.xmd')

        # Set CTF information to the micrographs to be displayed in the
        # picking list
        autopickProt.micDict = OrderedDict()
        micDict, _ = autopickProt._loadInputList()

        def _preprocessMic(mic, micRow):
            mic.setCTF(micDict[mic.getMicName()].getCTF())

        writeSetOfMicrographs(micSet, micStarFn,
                              preprocessImageRow=_preprocessMic)

        # Create a folder in extra to backup the original autopick star files
        backupDir = autopickProt._getExtraPath('wizard-backup')
        pw.utils.cleanPath(backupDir)
        pw.utils.makePath(backupDir)
        pw.utils.copyPattern(autopickProt._getExtraPath("*autopick.star"),
                            backupDir)

        cmd = '%s relion_autopick ' % pw.getScipionScript()
        #cmd += '--i extra/%(micrographName).star '
        cmd += '--i input_micrographs.star '
        cmd += '--threshold %(threshold) --min_distance %(ipd) '
        cmd += ' --max_stddev_noise %(maxStddevNoise) '
        cmd += ' --read_fom_maps'
        cmd += autopickProt.getAutopickParams()

        convertCmd = pw.join('apps', 'pw_convert.py')
        convertCmd += ' --coordinates --from relion --to xmipp '
        convertCmd += ' --input %s' % micSet.getFileName()
        convertCmd += ' --output %s' % coordsDir
        convertCmd += ' --extra %s' % autopickProt._getExtraPath()

        args = {
            "threshold": autopickProt.pickingThreshold,
            'min_distance': autopickProt.interParticleDistance,
            'autopickCommand': cmd,
            'preprocessCommand': 'rm -rf %s/*.pos' % coordsDir,
            'convertCmd': convertCmd,
            'protDir': autopickProt.getWorkingDir(),
            'maxStddevNoise': autopickProt.maxStddevNoise
        }

        pickerProps = os.path.join(coordsDir, 'picker.conf')

        f = open(pickerProps, "w")
        f.write("""
        parameters = ipd,threshold,maxStddevNoise
        ipd.value = %(min_distance)s
        ipd.label = Inter-particles distance (A)
        ipd.help = Particles closer together than this distance will be consider to be a single cluster. From each cluster, only one particle will be picked.
        threshold.value =  %(threshold)s
        threshold.label = Threshold
        threshold.help = Use lower thresholds to pick more particles (and more junk probably).
        maxStddevNoise.value = %(maxStddevNoise)s
        maxStddevNoise.label = Max. stddev noise
        maxStddevNoise.help = Prevent picking in carbon areas, useful values probably between 1.0 and 1.2, use -1 to switch it off
        runDir = %(protDir)s
        preprocessCommand = %(preprocessCommand)s
        autopickCommand = %(autopickCommand)s
        convertCommand = %(convertCmd)s
        hasInitialCoordinates = true
        doPickAll = true
        """ % args)
        f.close()
        os.environ['XMIPP_EXTRA_ALIASES'] = 'micrograph=rlnMicrographName'
        process = CoordinatesObjectView(autopickProt.getProject(), micStarFn,
                                        coordsDir, autopickProt,
                                        mode=CoordinatesObjectView.MODE_AUTOMATIC,
                                        pickerProps=pickerProps).show()
        process.wait()
        myprops = pw.utils.readProperties(pickerProps)

        # Check if the wizard changes were accepted or just canceled
        if myprops.get('applyChanges', 'false') == 'true':
            form.setVar('pickingThreshold', myprops['threshold.value'])
            form.setVar('interParticleDistance', myprops['ipd.value'])
            form.setVar('maxStddevNoise', myprops['maxStddevNoise.value'])
            # Change the run type now to 'Compute' after using the wizard
            # and (supposedly) optimized parameters
            form.setVar('runType', RUN_COMPUTE)
            # Mark the wizard was used
            setattr(autopickProt, 'wizardExecuted', True)
        else:
            # If the wizard was not execute, we should restore the original
            # autopick star files in case their were modified by the wizard
            pw.utils.copyPattern(os.path.join(backupDir, "*autopick.star"),
                                autopickProt._getExtraPath())


class Relion2PartDiameter(RelionPartMaskDiameterWizard):
    _targets = [(ProtRelion2Autopick, ['particleDiameter'])]

    def _getProtocolImages(self, protocol):
        return protocol.inputReferences

    def show(self, form):
        prot = form.protocol
        if prot.useInputReferences():
            if prot.getInputReferences() is None:
                form.showWarning("Please select the input references first. ")
            else:
                RelionPartMaskDiameterWizard.show(self, form)
        else: # Gaussian blobs
            form.showWarning("This wizard only works when using input "
                             "references, not Gaussian blobs. ")


class RelionWizLogPickParams(EmWizard):
    _targets = [(ProtRelionAutopickLoG, ['minDiameter',
                                         'maxDiameter',
                                         'threshold'])]

    def show(self, form):
        autopickProt = form.protocol
        project = autopickProt.getProject()
        micSet = autopickProt.getInputMicrographs()
        micfn = micSet.getFileName()
        coordsDir = project.getTmpPath(micSet.getName())
        print("coordsDir: ", coordsDir)
        params, minDiameter, maxDiameter, threshold = autopickProt._getPickArgs()

        pw.utils.cleanPath(coordsDir)
        pw.utils.makePath(coordsDir, 'extra')
        pickerProps = os.path.join(coordsDir, 'picker.conf')
        micStarFn = os.path.join(coordsDir, 'input_micrographs.star')

        def _postprocessMic(mic, micRow):
            micFn = mic.getFileName()
            micBase = os.path.basename(micFn)
            pw.utils.createLink(micFn, os.path.join(coordsDir, micBase))
            micRow.setValue(md.RLN_MICROGRAPH_NAME, micBase)

        writeSetOfMicrographs(micSet, micStarFn, postprocessImageRow=_postprocessMic)

        f = open(pickerProps, "w")

        #params = params.replace('--odir ""', '--odir extra')
        autopickCmd = "%s relion_autopick " % pw.getScipionScript()
        autopickCmd += ' --i input_micrographs.star '
        autopickCmd += params
        autopickCmd += ' --LoG_diam_min %(mind) '
        autopickCmd += ' --LoG_diam_max %(maxd) '
        autopickCmd += ' --LoG_adjust_threshold %(threshold) '

        args = {
            "convert": pw.join('apps', 'pw_convert.py'),
            'coordsDir': coordsDir,
            'micsSqlite': micSet.getFileName(),
            "minDiameter": minDiameter,
            "maxDiameter": maxDiameter,
            "threshold": threshold,
            'projDir': project.getPath(), #autopickProt.getWorkingDir(),
            "autopickCmd": autopickCmd
        }

        f.write("""
        parameters = mind,maxd,threshold
        mind.value = %(minDiameter)s
        mind.label = Min. Diameter (A)
        mind.help = The smallest allowed diameter for the blob-detection algorithm. This should correspond to the smallest size of your particles in Angstroms.
        maxd.value = %(maxDiameter)s
        maxd.label = Max. Diameter (A)
        maxd.help = The largest allowed diameter for the blob-detection algorithm. This should correspond to the largest size of your particles in Angstroms.
        threshold.value =  %(threshold)s
        threshold.label = Threshold
        threshold.help = Lower threshold -> more particles
        runDir = %(coordsDir)s
        autopickCommand = %(autopickCmd)s
        convertCommand = %(convert)s --coordinates --from relion --to xmipp --input  %(micsSqlite)s --output %(coordsDir)s --extra %(coordsDir)s/
        hasInitialCoordinates = false
        doPickAll = true
        """ % args)
        f.close()
        process = CoordinatesObjectView(autopickProt.getProject(), micfn, coordsDir, autopickProt,
                                        pickerProps=pickerProps).show()
        process.wait()
        myprops = pw.utils.readProperties(pickerProps)

        if myprops['applyChanges'] == 'true':
            form.setVar('minDiameter', myprops['mind.value'])
            form.setVar('maxDiameter', myprops['maxd.value'])
            form.setVar('threshold', myprops['threshold.value'])

