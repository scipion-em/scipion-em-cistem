# **************************************************************************
# *
# * ctffind on tilt series: a tilt series with no usable CTF is left out and
# * estimated again on Continue, and the run fails only when every one failed.
# *
# **************************************************************************

import os
import tempfile
import unittest
from types import SimpleNamespace
from unittest import mock

from pyworkflow.utils import weakImport

with weakImport('tomo'):
    from cistem.protocols.protocol_ts_ctffind import CistemProtTsCtffind, TS_STEPS


    class _Ts:
        def __init__(self, tsId):
            self._tsId = tsId

        def getTsId(self):
            return self._tsId


    def _md(tsId):
        return SimpleNamespace(ts=_Ts(tsId))


    class _CtfSet:
        def __init__(self, tsIds):
            self._tsIds = list(tsIds)

        def getUniqueValues(self, field):
            return list(self._tsIds)

        def __len__(self):
            return len(self._tsIds)


    class _Harness:
        _stepResultPath = CistemProtTsCtffind._stepResultPath
        _stepResultFile = CistemProtTsCtffind._stepResultFile
        _prepareFailedForRetry = CistemProtTsCtffind._prepareFailedForRetry
        _failedTsIds = CistemProtTsCtffind._failedTsIds
        _createOutput = CistemProtTsCtffind._createOutput
        closeStep = CistemProtTsCtffind.closeStep

        def __init__(self, root, published=()):
            self._root = root
            self.tsCtfMdList = [_md('TS_01'), _md('TS_02')]
            self._output = _CtfSet(published) if published else None
            self.errors = []
            self._closeOutputSet = mock.Mock()

        def _getExtraPath(self, *paths):
            return os.path.join(self._root, 'extra', *paths)

        def getOutputCtfTomoSet(self):
            return self._output

        def error(self, message):
            self.errors.append(message)


    class TestATiltSeriesWithoutCtfIsLeftOut(unittest.TestCase):
        def test_ItIsNotPublished(self):
            with tempfile.TemporaryDirectory() as tmp:
                prot = _Harness(tmp)
                prot._readTsCtfs = lambda mdObj: None

                prot._createOutput(_md('TS_01'))

            self.assertEqual({'TS_01'}, prot._failedTsIds())
            self.assertIsNone(prot._output, "Nothing was published for it.")
            self.assertTrue(prot.errors)

        def test_EveryTiltSeriesFailedFailsTheRunButClosesTheOutput(self):
            with tempfile.TemporaryDirectory() as tmp:
                prot = _Harness(tmp)
                prot._failedTsIds().update({'TS_01', 'TS_02'})

                with self.assertRaisesRegex(Exception, "failed for every one of the 2"):
                    prot.closeStep()

            prot._closeOutputSet.assert_called_once()

        def test_SomePublishedTiltSeriesAreNoFailure(self):
            with tempfile.TemporaryDirectory() as tmp:
                prot = _Harness(tmp, published=['TS_01'])
                prot._failedTsIds().add('TS_02')

                prot.closeStep()


    class TestContinueEstimatesAgainWhatFailed(unittest.TestCase):
        def test_OnlyTheStepsOfTheFailedTiltSeriesRunAgain(self):
            with tempfile.TemporaryDirectory() as tmp:
                prot = _Harness(tmp, published=['TS_01'])
                markers = {}
                for md in prot.tsCtfMdList:
                    for kind in TS_STEPS:
                        markers['%s_%s' % (kind, md.ts.getTsId())] = prot._stepResultFile(kind, md)
                markers['close'] = prot._stepResultFile('close')

                prot._prepareFailedForRetry()

                remaining = {name for name, path in markers.items() if os.path.exists(path)}

            self.assertEqual({'%s_TS_01' % kind for kind in TS_STEPS}, remaining,
                             "TS_02 is not published: its steps, and the closing one, run again.")


if __name__ == '__main__':
    unittest.main()
