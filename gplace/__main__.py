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

import logging
import sys
import traceback

from gplace import __author__, __copyright__, __version__
from gplace.biolib_lite.exceptions import BioLibError
from gplace.biolib_lite.logger import logger_setup
from gplace.cli import get_main_parser
from gplace.exceptions import GPlaceExit, GPlaceException
from gplace.main import OptionsParser


def print_help():
    print('''\

              ...::: gPlace v%s :::...

  Placement:
    place         -> Place aligned queries on a reference package (same result as pplacer).
    place_batch   -> Place many alignments, pooling queries per reference package.
    compare       -> Compare two jplace files placement by placement.

  Planning:
    benchmark     -> Find the fastest --cpus, and the memory and time a job needs (e.g. for SLURM).

  Testing:
    test          -> Place the bundled genomes and check the output is identical to pplacer.
    check_install -> Check the engine is installed and bit-exact.

  Use: gplace <command> -h for command specific help
    ''' % __version__)


def main():
    # -------------------------------------------------
    # get and check options
    args = None
    if len(sys.argv) == 1:
        print_help()
        sys.exit(0)
    elif sys.argv[1] in {'-v', '--v', '-version', '--version'}:
        print(f"gplace: version {__version__} {__copyright__} {__author__}")
        sys.exit(0)
    elif sys.argv[1] in {'-h', '--h', '-help', '--help'}:
        print_help()
        sys.exit(0)
    else:
        args = get_main_parser().parse_args()

    # setup logger
    log_dir = None
    if hasattr(args, 'out_dir') and args.out_dir:
        log_dir = args.out_dir
    elif hasattr(args, 'out_file') and args.out_file:
        import os
        log_dir = os.path.dirname(os.path.abspath(args.out_file))
    logger_setup(log_dir, "gplace.log", "gPlace", __version__, False,
                 hasattr(args, 'debug') and args.debug)
    logger = logging.getLogger('timestamp')

    # -------------------------------------------------
    # do what we came here to do
    try:
        gp_parser = OptionsParser(__version__, log_dir)
        gp_parser.parse_options(args)
    except SystemExit:
        logger.error('Controlled exit resulting from early termination.')
        sys.exit(1)
    except KeyboardInterrupt:
        logger.error('Controlled exit resulting from interrupt signal.')
        sys.exit(1)
    except GPlaceExit as e:
        if len(str(e)) > 0:
            logger.error('{}'.format(e))
        logger.error('Controlled exit resulting from an unrecoverable error or warning.')
        sys.exit(1)
    except (GPlaceException, BioLibError) as e:
        msg = 'Controlled exit resulting from an unrecoverable error or warning.\n\n'
        msg += '=' * 80 + '\n'
        msg += 'EXCEPTION: {}\n'.format(type(e).__name__)
        msg += '  MESSAGE: {}\n'.format(e)
        msg += '_' * 80 + '\n\n'
        msg += traceback.format_exc()
        msg += '=' * 80
        logger.error(msg)
        sys.exit(1)
    except Exception as e:
        msg = 'Uncontrolled exit resulting from an unexpected error.\n\n'
        msg += '=' * 80 + '\n'
        msg += 'EXCEPTION: {}\n'.format(type(e).__name__)
        msg += '  MESSAGE: {}\n'.format(e)
        msg += '_' * 80 + '\n\n'
        msg += traceback.format_exc()
        msg += '=' * 80
        logger.error(msg)
        sys.exit(1)


if __name__ == '__main__':
    main()
