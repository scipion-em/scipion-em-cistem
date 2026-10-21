# **************************************************************************
# *
# * Two micrographs whose files share a basename (in different folders)
# * must never share an artefact: one would overwrite the other's results.
# *
# **************************************************************************

import os
import tempfile
import unittest

from cistem.protocols.protocol_ctffind import CistemProtCTFFind
from cistem.protocols.protocol_picking import CistemProtFindParticles
from cistem.protocols.protocol_streaming_base import CistemStreamingBase


class _Mic:
    def __init__(self, objId, fileName, yDim=100):
        self._objId, self._fileName, self._yDim = objId, fileName, yDim

    def getObjId(self):
        return self._objId

    def getFileName(self):
        return self._fileName

    def getMicName(self):
        return os.path.basename(self._fileName)

    def getYDim(self):
        return self._yDim


class _Folders(CistemStreamingBase):
    def __init__(self, root):
        self._root = root

    def _getExtraPath(self, *paths):
        return os.path.join(self._root, 'extra', *paths)

    def _getTmpPath(self, *paths):
        return os.path.join(self._root, 'tmp', *paths)


class _CtffindFolders(_Folders):
    _getMicExtra = CistemProtCTFFind._getMicExtra
    _getPsdPath = CistemProtCTFFind._getPsdPath
    _getCtfOutPath = CistemProtCTFFind._getCtfOutPath
    _getCtfAvrotPath = CistemProtCTFFind._getCtfAvrotPath


class _PickingFolders(_Folders):
    _getPltFn = CistemProtFindParticles._getPltFn
    _getLogFn = CistemProtFindParticles._getLogFn
    _getStackFn = CistemProtFindParticles._getStackFn


class TestEveryMicrographHasItsOwnArtefacts(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.root = self._tmp.name
        os.makedirs(os.path.join(self.root, 'extra'))
        self.mics = [_Mic(1, '/data/a/mic.mrc'), _Mic(2, '/data/b/mic.mrc')]

    def tearDown(self):
        self._tmp.cleanup()

    def test_CtffindResultsOfSameNamedMicrographsDoNotCollide(self):
        protocol = _CtffindFolders(self.root)

        for getPath in (protocol._getPsdPath, protocol._getCtfOutPath,
                        protocol._getCtfAvrotPath):
            self.assertNotEqual(getPath(self.mics[0]), getPath(self.mics[1]))

    def test_PickingFilesOfSameNamedMicrographsDoNotCollide(self):
        protocol = _PickingFolders(self.root)

        for getPath in (protocol._getPltFn, protocol._getLogFn):
            self.assertNotEqual(getPath(self.mics[0]), getPath(self.mics[1]))

    def test_AProjectWrittenBeforeKeepsItsResults(self):
        # Only the old, unscoped file exists: Continue must still find it.
        protocol = _CtffindFolders(self.root)
        legacy = os.path.join(self.root, 'extra', 'mic_ctf.txt')
        open(legacy, 'w').close()

        self.assertEqual(legacy, protocol._getCtfOutPath(self.mics[0]))

    def test_PickingReadsEachMicrographFromItsOwnFile(self):
        protocol = _PickingFolders(self.root)
        for mic, (x, y) in zip(self.mics, ((10, 20), (30, 40))):
            with open(protocol._getPltFn(mic), 'w') as fh:
                fh.write('%d %d\n' % (y, x))

        coords = _CoordSet()
        CistemProtFindParticles.readCoordsFromMics(protocol, None, self.mics,
                                                   coords)

        self.assertEqual([(1, 10.0), (2, 30.0)],
                         [(c.getMicId(), c.getX()) for c in coords.items])


class _CoordSet:
    def __init__(self):
        self.items = []

    def getBoxSize(self):
        return 64

    def append(self, coord):
        self.items.append(coord)


if __name__ == '__main__':
    unittest.main()
