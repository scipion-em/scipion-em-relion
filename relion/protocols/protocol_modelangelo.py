# **************************************************************************
# *
# * Authors:     Roberto marabini
# *
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
import os.path
from typing import List

from pwem.objects import Volume, AtomStruct, VolumeMask
from pyworkflow.protocol import params, LEVEL_ADVANCED
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol
from pyworkflow.protocol import GPU_LIST, USE_GPU
from pwem.convert.headers import Ccp4Header

from relion import Plugin


class ProtRelionModelAngelo(EMProtocol):
    """
    Relion protocol for continuous flexibility analysis.

    As of release 5.0, Relion comes with a machine-learning approach for
    the analysis of molecular motions and flexibility called DynaMight.
    DynaMight will fit molecular motions for each experimental particle image
    as a 3D deformation field, which is learnt using a variational auto-encoder.
    It also implements functionality to calculate a pseudo-inverse 3D deformation
    field that can then be used in a deformed weighted backprojection algorithm
    to obtain an improved 3D reconstruction of the consensus structure.

    """
    """
    ModelAngelo is an automatic atomic model building program for cryo-EM maps.
    With or without providing a sequence.
    """
    _label = 'modelangelo'

    OUTPUT_NAME = "pruned"
    OUTPUT_RAW_NAME = "raw"

    _possibleOutputs = {
        OUTPUT_NAME: AtomStruct,
        OUTPUT_RAW_NAME: AtomStruct
    }

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputVolume', params.PointerParam,
                      pointerClass=Volume,
                      label='Refined volume', important=True,
                      help='Refined cryo-em map.')

        form.addParam('inputSequenceS', params.MultiPointerParam,
                      pointerClass="Sequence", allowsNull=True, important=True,
                      label='Protein sequences',
                      help="Include here one or more sequences to be modeled\n"
                           "Leave empty to use the *model_no_seq* option.")

        form.addParam('inputMask', params.PointerParam,
                      pointerClass=VolumeMask,
                      label='Volume mask', allowsNull=True, important=True,
                      help='Search will be done inside the mask.\n'
                           'That is, voxels inside the mask should be NON zero.')

        form.addHidden(USE_GPU, params.BooleanParam, default=True,
                       label="Use GPU for execution",
                       help="This protocol has both CPU and GPU implementation. "
                            "Select the one you want to use.")

        form.addHidden(GPU_LIST, params.StringParam, default='0',
                       expertLevel=LEVEL_ADVANCED,
                       label="Choose GPU ID (single one)",
                       help="GPU device to be used")

        form.addParam('configFile', params.FileParam,
                      label="Configuration File",
                      default="",
                      expertLevel=LEVEL_ADVANCED,
                      help="""This option is only for VERY advanced users.\n
Below is an example of a config file:

{
  "standardize_mrc_args":
    {
      "target_voxel_size": 1.5,
      "crop_z": 0,
      "bfactor_to_apply": 0,
      "auto_mask": false
    },
  "ca_infer_args":
    {
      "model_checkpoint": "chkpt.torch",
      "bfactor": 0,
      "batch_size": 4,
      "stride": 16,
      "dont_mask_input": true,
      "threshold": 0.05,
      "save_real_coordinates": false,
      "save_cryo_em_grid": false,
      "do_nucleotides": false,
      "save_backbone_trace": false,
      "save_ca_grid": false,
      "crop": 6
    },
  "gnn_infer_args":
    {
      "num_rounds": 3,
      "crop_length": 200,
      "repeat_per_residue": 3,
      "esm_model": "esm1b_t33_650M_UR50S",
      "aggressive_pruning": false,
      "seq_attention_batch_size": 200
    }
}
""")

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.predictStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions ------------------------------
    def convertInputStep(self):
        """ convert 3D maps to MRC '.mrc' format
            with extension mrc, map extension is not handled by modelangelo
        """
        vol = self.inputVolume.get()
        inVolName = vol.getFileName()
        if inVolName.endswith(".mrc"):
            self.newFn = inVolName
        else:
            self.newFn = self._getExtraPath("inputVol.mrc")
            origin = vol.getOrigin(force=True).getShifts()
            sampling = vol.getSamplingRate()
            Ccp4Header.fixFile(inVolName, self.newFn, origin, sampling, Ccp4Header.START)  # ORIGIN

    def predictStep(self):
        seqs = self.inputSequenceS
        mask = self.inputMask.get()
        configFile = self.configFile.get()

        args = []
        if seqs:
            fasta = self.createInputFastaFile(seqs)
            args.extend(["build", "--fasta-path", fasta])
        else:
            args.append("build_no_seq")

        args.extend(["--volume-path", self.newFn,
                     "--output-dir", self._getExtraPath()])

        if mask:
            args.extend(["--mask-path", mask.getFileName()])

        # Gpu or cpu
        args.extend(["--device", ("%s" % self.getGpuList()[0]) if self.useGpu else "cpu"])

        if configFile:
            args.extend(["-c", configFile])

        try:
            # Call model angelo:
            self.runProgram(args)
            # self.runJob(Plugin.getModelAngeloCmd(), args)
        except Exception:
            # Modelangelo does not show error in the stdout, nor stderr we
            # should go and read the error information from a log file
            with open(self._getExtraPath("model_angelo.log")) as log:
                for line in log.read().splitlines():
                    self.error(line)
            self.info("ERROR: %s." % line)
            raise ChildProcessError("Model angelo has failed: %s. See error log "
                                    "for more details." % line) from None

    def createOutputStep(self):
        """Register atomic models, raw and pruned"""
        # check if files exists before registering
        # I think build_no_seq creates a single output file (no raw file)
        if os.path.exists(self._getExtraPath('extra_raw.cif')):
            self._registerAtomStruct(self.OUTPUT_RAW_NAME, self._getExtraPath('extra_raw.cif'))
        self._registerAtomStruct(self.OUTPUT_NAME, self._getExtraPath('extra.cif'))

    # --------------------------- INFO functions -----------------------------------
    def _validate(self):
        errors = []
        gpus = self.getGpuList()

        if len(gpus) > 1:
            errors.append('Only one GPU can be used.')

        return errors

    # -------------------------- UTILS functions ------------------------------
    def createInputFastaFile(self, seqs):
        """ Get sequence as string and create the corresponding fasta file. """

        fastaFileName = self._getExtraPath('sequence.fasta')

        with open(fastaFileName, "w") as f:
            for seq in seqs:
                s = seq.get()
                f.write(f"> {s.getId()}\n")
                f.write(f"{s.getSequence()}\n")

        return fastaFileName

    def _registerAtomStruct(self, name, path):
        if not os.path.exists(path):
            raise FileNotFoundError("Output %s not found." % path)

        output = AtomStruct(filename=path)
        self._defineOutputs(**{name: output})
        self._defineSourceRelation(self.inputVolume, output)

        seqs = self.inputSequenceS
        if seqs:
            for seq in seqs:
                self._defineSourceRelation(seq, output)

    # -------------------------- UTILS functions ------------------------------
    def runProgram(self, params: List[str]) -> None:
        program = "relion_python_modelangelo"
        self.runJob(f"{Plugin.getActivationCmd()} && {program}",
                    " ".join(params))

    def _getEnviron(self):
        env = Plugin.getEnviron()
        if 'LD_LIBRARY_PATH' in env:
            # this is required to avoid conflict btw DynaMight Qt5 libs
            # and system Qt libs
            del env['LD_LIBRARY_PATH']

        return env