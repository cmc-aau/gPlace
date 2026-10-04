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

"""Whole-file reads and writes that close their file (no ResourceWarning)."""

import json


def read_text(path):
    with open(path) as fh:
        return fh.read()


def read_json(path):
    with open(path) as fh:
        return json.load(fh)


def write_text(path, text):
    with open(path, 'w') as fh:
        fh.write(text)


def write_json(path, obj, **kw):
    with open(path, 'w') as fh:
        json.dump(obj, fh, **kw)
