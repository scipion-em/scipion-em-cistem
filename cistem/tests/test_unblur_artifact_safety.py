# **************************************************************************
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
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
"""Artefacts Unblur owns: the failure report and the per-movie folder.

The failure report is a payload the viewer reads back, not completion
state, so dropping it breaks a user-visible feature without any of the
backend independence it was traded for. The movie folder is removed
through a shell, which turns a space in a project path into the deletion
of something else entirely.
"""
import os
import shutil
import tempfile
import unittest

from cistem.protocols.protocol_picking import CistemProtFindParticles
from cistem.protocols.protocol_unblur import CistemProtUnblur


class _Movie:
    def __init__(self, objId):
        self._objId = objId

    def getObjId(self):
        return self._objId


class _FailureReportHarness:
    """Exercises the real read/write pair against a temporary extra dir."""

    _writeFailedList = CistemProtUnblur._writeFailedList
    _readFailedList = CistemProtUnblur._readFailedList
    _getAllFailed = CistemProtUnblur._getAllFailed

    def __init__(self, extraDir):
        self._extraDir = extraDir

    def _getExtraPath(self, *paths):
        return os.path.join(self._extraDir, *paths)


class TestFailedMovieReportSurvives(unittest.TestCase):

    def setUp(self):
        self.extraDir = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.extraDir, True)
        self.harness = _FailureReportHarness(self.extraDir)

    def testAFailedMovieCanBeReadBack(self):
        self.harness._writeFailedList([_Movie(7), _Movie(9)])

        self.assertEqual(
            self.harness._readFailedList(),
            [7, 9],
            "The viewer offers 'show failed movies' by reading this back. "
            "Dropping the write makes it answer 'no failed movies' for a "
            "run where movies did fail.",
        )

    def testFailuresAccumulateAcrossBatches(self):
        self.harness._writeFailedList([_Movie(1)])
        self.harness._writeFailedList([_Movie(4)])

        self.assertEqual(self.harness._readFailedList(), [1, 4])

    def testNoFailuresReadsBackEmpty(self):
        self.assertEqual(self.harness._readFailedList(), [])

    def testTheReportLivesUnderTheProtocolExtraDir(self):
        self.harness._writeFailedList([_Movie(3)])

        written = os.listdir(self.extraDir)

        self.assertEqual(
            len(written),
            1,
            "The report is this protocol's own payload file and belongs in "
            "its own working directory, nowhere else.",
        )


class _PickingFailureHarness:
    """FindParticles keeps the same report, written when a pick fails."""

    _writeFailedList = CistemProtFindParticles._writeFailedList
    _readFailedList = CistemProtFindParticles._readFailedList
    _getAllFailed = CistemProtFindParticles._getAllFailed

    def __init__(self, extraDir):
        self._extraDir = extraDir

    def _getExtraPath(self, *paths):
        return os.path.join(self._extraDir, *paths)


class TestFailedMicrographReportSurvives(unittest.TestCase):

    def setUp(self):
        self.extraDir = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.extraDir, True)
        self.harness = _PickingFailureHarness(self.extraDir)

    def testAFailedMicrographCanBeReadBack(self):
        self.harness._writeFailedList([_Movie(12)])

        self.assertEqual(
            self.harness._readFailedList(),
            [12],
            "A micrograph whose picking failed produces no coordinates, "
            "which on its own is indistinguishable from one that simply "
            "had no particles. The report is the only durable record.",
        )

    def testNoFailuresReadsBackEmpty(self):
        self.assertEqual(self.harness._readFailedList(), [])


class _CleanupHarness:
    """Exercises the real folder cleanup against a temporary workspace."""

    _cleanMovieFolder = CistemProtUnblur._cleanMovieFolder

    def __init__(self, tmpDir):
        self._tmpDir = tmpDir
        self.messages = []

    def _getTmpPath(self, *paths):
        return os.path.join(self._tmpDir, *paths)

    def info(self, message):
        self.messages.append(message)

    def warning(self, message):
        self.messages.append(message)


class TestMovieFolderCleanupStaysInside(unittest.TestCase):

    def setUp(self):
        self.root = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.root, True)

    def _workspace(self, name):
        path = os.path.join(self.root, name)
        os.makedirs(path)
        return path

    def testAMovieFolderIsRemoved(self):
        tmpDir = self._workspace('tmp')
        harness = _CleanupHarness(tmpDir)
        movieFolder = os.path.join(tmpDir, 'movie_000001')
        os.makedirs(movieFolder)

        harness._cleanMovieFolder(movieFolder)

        self.assertFalse(os.path.exists(movieFolder))

    def testASpaceInThePathDoesNotDeleteASibling(self):
        """The folder name is handed to a shell, which splits on spaces."""
        tmpDir = self._workspace('my project/tmp')
        harness = _CleanupHarness(tmpDir)
        movieFolder = os.path.join(tmpDir, 'movie_000001')
        os.makedirs(movieFolder)

        # What the shell would take as the second argument.
        sibling = os.path.join(self.root, 'my')
        os.makedirs(os.path.join(sibling, 'keep-me'), exist_ok=True)

        harness._cleanMovieFolder(movieFolder)

        self.assertTrue(
            os.path.exists(os.path.join(sibling, 'keep-me')),
            "A space in the project path turned the cleanup into the "
            "deletion of an unrelated directory.",
        )
        self.assertFalse(os.path.exists(movieFolder))

    def testAFolderOutsideTheWorkspaceIsRefused(self):
        tmpDir = self._workspace('tmp')
        harness = _CleanupHarness(tmpDir)
        outsider = self._workspace('not-mine')

        harness._cleanMovieFolder(outsider)

        self.assertTrue(
            os.path.exists(outsider),
            "Cleanup must only ever touch paths this protocol controls.",
        )

    def testTheWorkspaceItselfIsRefused(self):
        tmpDir = self._workspace('tmp')
        harness = _CleanupHarness(tmpDir)

        harness._cleanMovieFolder(tmpDir)

        self.assertTrue(
            os.path.exists(tmpDir),
            "Removing the working directory itself would take every other "
            "movie's in-flight data with it.",
        )


if __name__ == '__main__':
    unittest.main()
