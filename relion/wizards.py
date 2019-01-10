# **************************************************************************
# *
# * Authors:     Jose Gutierrez (jose.gutierrez@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

import pyworkflow.utils as pwutils
from pyworkflow.em import *
from pyworkflow.em.wizard import *

from relion.constants import *
from relion.convert import writeSetOfMicrographs
from relion.protocols import (
    ProtRelionClassify3D, ProtRelionRefine3D, ProtRelionClassify2D,
    ProtRelionPreprocessParticles,
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
        pwutils.cleanPath(coordsDir)
        pwutils.makePath(coordsDir)

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
        pwutils.cleanPath(backupDir)
        pwutils.makePath(backupDir)
        pwutils.copyPattern(autopickProt._getExtraPath("*autopick.star"),
                            backupDir)

        cmd = '%s relion_autopick ' % pw.getScipionScript()
        #cmd += '--i extra/%(micrographName).star '
        cmd += '--i input_micrographs.star '
        cmd += '--threshold %(threshold) --min_distance %(ipd) '
        cmd += ' --max_stddev_noise %(maxStddevNoise) '
        cmd += ' --read_fom_maps'
        cmd += autopickProt.getAutopickParams()

        convertCmd = pwutils.join('apps', 'pw_convert.py')
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
        ipd.help = Minimum distance (in Angstroms) between particles
        threshold.value =  %(threshold)s
        threshold.label = Threshold
        threshold.help = Picking threshold
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
        process = CoordinatesObjectView(autopickProt.getProject(), micStarFn,
                                        coordsDir, autopickProt,
                                        mode=CoordinatesObjectView.MODE_AUTOMATIC,
                                        pickerProps=pickerProps).show()
        process.wait()
        myprops = pwutils.readProperties(pickerProps)

        # Check if the wizard changes were accepted or just canceled
        if myprops.get('applyChanges', 'false') == 'true':
            form.setVar('pickingThreshold', myprops['threshold.value'])
            form.setVar('interParticleDistance', myprops['ipd.value'])
            form.setVar('maxStddevNoise', myprops['maxStddevNoise.value'])
            # Change the run type now to 'Compute' after using the wizard
            # and (supposedly) optimized parameters
            form.setVar('runType', ProtRelion2Autopick.RUN_COMPUTE)
            # Mark the wizard was used
            setattr(autopickProt, 'wizardExecuted', True)
        else:
            # If the wizard was not execute, we should restore the original
            # autopick star files in case their were modified by the wizard
            pwutils.copyPattern(os.path.join(backupDir, "*autopick.star"),
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
