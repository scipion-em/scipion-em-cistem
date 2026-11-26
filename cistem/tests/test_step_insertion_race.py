# **************************************************************************
# *
# * The generator inserts steps while the executor runs them: a step must
# * find its micrograph as soon as it exists.
# *
# **************************************************************************

import unittest
from collections import OrderedDict

from pwem.protocols.protocol import EMProtocol

from cistem.protocols.protocol_ctffind import CistemProtCTFFind
from cistem.protocols.protocol_picking import CistemProtFindParticles


class _Mic:
    def __init__(self, name):
        self._name = name

    def getMicName(self):
        return self._name


class _RacingExecutorHarness:
    """Runs each step the moment it is inserted, as the executor may."""
    _insertNewMics = EMProtocol._insertNewMics

    def __init__(self, batchSize):
        self.micDict = OrderedDict()
        self.initialIds = []
        self.streamClosed = True
        self._batchSize = batchSize
        self.missing = []
        self.ran = []

    def _getStreamingBatchSize(self):
        return self._batchSize

    def _insertFunctionStep(self, funcName, names, *args, **kwargs):
        names = names if isinstance(names, list) else [names]
        self.ran.extend(names)
        self.missing.extend(name for name in names if name not in self.micDict)
        return len(self.ran)


class _CtffindHarness(_RacingExecutorHarness):
    _insertNewMicsSteps = CistemProtCTFFind._insertNewMicsSteps
    _insertCtfStep = CistemProtCTFFind._insertCtfStep
    _insertCtfListStep = CistemProtCTFFind._insertCtfListStep

    def _getCtfArgs(self):
        return []


class _PickingHarness(_RacingExecutorHarness):
    _insertNewMicsSteps = CistemProtFindParticles._insertNewMicsSteps
    _insertPickMicrographStep = CistemProtFindParticles._insertPickMicrographStep
    _insertPickMicrographListStep = CistemProtFindParticles._insertPickMicrographListStep

    def _getPickArgs(self):
        return []


class TestStepsFindTheirMicrograph(unittest.TestCase):

    def _check(self, harnessClass, batchSize):
        harness = harnessClass(batchSize)
        mics = [_Mic('mic_%d.mrc' % i) for i in range(5)]
        harness._insertNewMicsSteps(mics)
        self.assertEqual(sorted(harness.ran), sorted(m.getMicName() for m in mics))
        self.assertEqual(harness.missing, [])
        self.assertEqual(list(harness.micDict), [m.getMicName() for m in mics])

    def test_CtffindOneByOne(self):
        self._check(_CtffindHarness, 1)

    def test_CtffindInBatches(self):
        self._check(_CtffindHarness, 2)

    def test_PickingOneByOne(self):
        self._check(_PickingHarness, 1)

    def test_PickingInBatches(self):
        self._check(_PickingHarness, 2)


if __name__ == '__main__':
    unittest.main()
