###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################


class GPlaceException(Exception):
    """ Base exception for all gPlace exceptions thrown in this project. """

    def __init__(self, message=''):
        Exception.__init__(self, message)


class GPlaceExit(Exception):
    """Raised when gPlace is to quietly exit."""

    def __init__(self, message=''):
        Exception.__init__(self, message)


class EngineNotBuilt(GPlaceException):
    """ Thrown when the compiled placement engine cannot be found. """

    def __init__(self, message=''):
        GPlaceException.__init__(self, message)


class EngineNotExact(GPlaceException):
    """ Thrown when the compiled engine fails a bit-exactness self-test. """

    def __init__(self, message=''):
        GPlaceException.__init__(self, message)


class RefpkgMalformed(GPlaceException):
    """ Thrown when a reference package is missing files or uses an unsupported model. """

    def __init__(self, message=''):
        GPlaceException.__init__(self, message)


class AlignmentMalformed(GPlaceException):
    """ Thrown when a query alignment has the wrong length or unknown characters. """

    def __init__(self, message=''):
        GPlaceException.__init__(self, message)


class BatchfileMalformed(GPlaceException):
    """ Thrown when the format of the batchfile is malformed. """

    def __init__(self, message=''):
        GPlaceException.__init__(self, message)


class FileNotFound(GPlaceException):
    """ Thrown when a file is not found. """

    def __init__(self, message=''):
        GPlaceException.__init__(self, message)
