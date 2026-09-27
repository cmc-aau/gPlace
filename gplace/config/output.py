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

from os.path import join

PATH_LOG = 'gplace.log'
PATH_WARNINGS = 'gplace.warnings.log'

# gplace test
DIR_TEST = 'test'
PATH_TEST_JPLACE = join(DIR_TEST, 'gplace_test.jplace')
PATH_TEST_REPORT = join(DIR_TEST, 'gplace_test_report.tsv')

# gplace place_batch
PATH_BATCH_SUMMARY = 'gplace.batch_summary.tsv'

JPLACE_FIELDS = ['distal_length', 'edge_num', 'like_weight_ratio', 'likelihood', 'pendant_length']
