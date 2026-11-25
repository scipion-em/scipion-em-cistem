# **************************************************************************
# *
# * Resampling tilt series: one that fails is left out and resampled again
# * on Continue, and the run fails only when every one failed.
# *
# **************************************************************************

import os
import tempfile
import unittest
from unittest import mock

from pyworkflow.utils import weakImport

with weakImport('tomo'):
    from cistem.protocols.protocol_ts_resample import CistemProtTsResample


    class _Harness:
        runTsResample = CistemProtTsResample.runTsResample
        _getOutputFn = CistemProtTsResample._getOutputFn
        _stepResultPath = CistemProtTsResample._stepResultPath
        _stepResultFile = CistemProtTsResample._stepResultFile
        _prepareFailedForRetry = CistemProtTsResample._prepareFailedForRetry

        def __init__(self, root):
            self._root = root
            self.newXsize = 100
            self.newYsize = 100
            self.errors = []

        def _getExtraPath(self, *paths):
            return os.path.join(self._root, 'extra', *paths)

        def error(self, message):
            self.errors.append(message)


    class TestAFailedResamplingIsLeftOut(unittest.TestCase):
        def test_ItDoesNotFailTheStep(self):
            with tempfile.TemporaryDirectory() as tmp:
                prot = _Harness(tmp)
                prot.runJob = mock.Mock(side_effect=RuntimeError("resample crashed"))

                with mock.patch('cistem.protocols.protocol_ts_resample.Plugin.getProgram',
                                return_value='resample'):
                    result = prot.runTsResample('/data/TS_01.mrcs')

                self.assertTrue(os.path.exists(result))
            self.assertTrue(prot.errors)

        def test_ContinueResamplesAgainWhatHasNoResult(self):
            with tempfile.TemporaryDirectory() as tmp:
                prot = _Harness(tmp)
                done, failed = '/data/TS_01.mrcs', '/data/TS_02.mrcs'
                markers = {name: prot._stepResultFile(name) for name in (done, failed, 'output')}
                open(prot._getOutputFn(done), 'w').close()

                prot._prepareFailedForRetry([done, failed])

                remaining = {name for name, path in markers.items() if os.path.exists(path)}

            self.assertEqual({done}, remaining,
                             "TS_02 has no resampled stack: its step, and the output one, run again.")


if __name__ == '__main__':
    unittest.main()
