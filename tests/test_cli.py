import gzip
import io
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stderr

import gplace
from gplace.cli import get_main_parser
from gplace.engine import testsets
from gplace.engine.jplace import compare_jplace
from gplace.main import OptionsParser


class TestCli(unittest.TestCase):

    def test_subcommands(self):
        p = get_main_parser()
        a = p.parse_args(['place', '--msa_file', 'q.fa', '--refpkg', 'r', '--out_file', 'o.json', '--cpus', '4'])
        self.assertEqual((a.subparser_name, a.cpus), ('place', '4'))
        self.assertEqual(p.parse_args(['test', '--cpus', 'auto']).cpus, 'auto')
        a = p.parse_args(['place_batch', '--batchfile', 'b.tsv', '--out_dir', 'out'])
        self.assertEqual(a.subparser_name, 'place_batch')
        for cmd in (['compare', '--jplace_a', 'a', '--jplace_b', 'b'], ['test'], ['check_install'], ['build_engine'],
                    ['benchmark', '--refpkg', 'r', '--queries', '100']):
            self.assertEqual(p.parse_args(cmd).subparser_name, cmd[0])

    def test_model_options(self):
        p = get_main_parser()
        place = ['place', '--msa_file', 'q.fa', '--refpkg', 'r', '--out_file', 'o.json']
        self.assertFalse(p.parse_args(place).fix_gamma_underflow)
        self.assertTrue(p.parse_args(place + ['--fix_gamma_underflow']).fix_gamma_underflow)
        self.assertTrue(p.parse_args(['place_batch', '--batchfile', 'b', '--out_dir', 'o',
                                      '--fix_gamma_underflow']).fix_gamma_underflow)
        self.assertTrue(p.parse_args(['test', '--models']).models)

    def test_tiling_options(self):
        p = get_main_parser()
        place = ['place', '--msa_file', 'q.fa', '--refpkg', 'r', '--out_file', 'o.json']
        a = p.parse_args(place)
        self.assertEqual((a.tile_edges, a.max_memory), (None, None))
        self.assertEqual(p.parse_args(place + ['--tile_edges', '10000']).tile_edges, 10000)
        self.assertEqual(p.parse_args(place + ['--max_memory', '64G']).max_memory, '64G')
        for cmd in (['place_batch', '--batchfile', 'b.tsv', '--out_dir', 'out'], ['test'],
                    ['benchmark', '--refpkg', 'r', '--queries', '100']):
            self.assertEqual(p.parse_args(cmd + ['--max_memory', '32G']).max_memory, '32G')
        with self.assertRaises(SystemExit), redirect_stderr(io.StringIO()):
            p.parse_args(place + ['--tile_edges', '10000', '--max_memory', '64G'])


class TestPplacerCommandLine(unittest.TestCase):
    """`gplace place` takes pplacer's command line (pplacer v1.1.alpha19's pplacer.ml / prefs.ml)."""

    def test_pplacer_names(self):
        p = get_main_parser()
        a = p.parse_args(['place', '-c', 'r.refpkg', '-o', 'o.jplace', '-j', '4', '--keep-at-most', '3', 'q.fasta'])
        self.assertEqual((a.refpkg, a.out_file, a.cpus, a.keep_at_most, a.alignment, a.msa_file),
                         ('r.refpkg', 'o.jplace', '4', 3, ['q.fasta'], None))
        a = p.parse_args(['place', '-c', 'r', '--out-dir', 'od', 'q1.fasta', 'q2.fasta'])
        self.assertEqual((a.out_dir, a.alignment, a.out_file), ('od', ['q1.fasta', 'q2.fasta'], None))
        a = p.parse_args(['place', '-c', 'r', '--mmap-file', 'm', '--verbosity', '0', '--timing', 'q.fa'])
        self.assertEqual((a.mmap_file, a.verbosity, a.timing, a.pplacer_unsupported), ('m', 0, True, []))

    def test_unsupported_options(self):
        p = get_main_parser()
        a = p.parse_args(['place', '-c', 'r', '--map-identity', '--groups', '3', '--fig-tree', 'f', 'q.fa'])
        self.assertEqual(a.pplacer_unsupported, ['--map-identity', '--groups', '--fig-tree'])
        a = p.parse_args(['place', '-c', 'r', '--groups', '1', '--fig-cutoff', '0', '--fantasy', '0',
                          '--fantasy-frac', '0.2', '--map-mrca-min', '0.9', 'q.fa'])
        self.assertEqual(a.pplacer_unsupported, [])          # pplacer's defaults: no effect

    def test_jplace_names(self):
        for path, name in (('q.fasta', 'q.jplace'), ('dir/q.fasta.gz', 'q.jplace'), ('q.v2.fasta', 'q.v2.jplace'),
                           ('q.sto', 'q.jplace'), ('noext', 'noext.jplace')):
            self.assertEqual(OptionsParser.pplacer_jplace_name(path), name)


class TestPplacerCommandLineRun(unittest.TestCase):
    """pplacer's command line, run: output names and content against pplacer's output."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        with open(os.path.join(testsets.TEST_DATA, 'queries.fasta'), 'rb') as src, \
                gzip.open(os.path.join(self.tmp, 'q.fasta.gz'), 'wb') as dst:
            dst.write(src.read())
        shutil.copy(os.path.join(testsets.TEST_DATA, 'queries.fasta'), os.path.join(self.tmp, 'queries.fasta'))

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def gplace(self, *args):
        env = dict(os.environ, PYTHONPATH=os.path.dirname(os.path.dirname(os.path.abspath(gplace.__file__))))
        return subprocess.run([sys.executable, '-m', 'gplace', 'place', *args], cwd=self.tmp, env=env,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

    def assertAsPplacer(self, path):
        r = compare_jplace(path, os.path.join(testsets.TEST_DATA, 'expected_pplacer.jplace'))
        self.assertTrue(r['same_tree'] and r['same_records'])
        self.assertEqual(r['identical'], r['genomes'])

    def test_files_named_as_pplacer(self):
        p = self.gplace('-j', '2', '-c', testsets.MINI, 'queries.fasta', '--out-dir', 'od', 'q.fasta.gz')
        self.assertEqual(p.returncode, 0, p.stdout)
        for name in ('queries.jplace', 'q.jplace'):
            self.assertAsPplacer(os.path.join(self.tmp, 'od', name))

    def test_o_and_errors(self):
        p = self.gplace('-c', testsets.MINI, '-j', '2', '-o', 'x.json', 'queries.fasta')
        self.assertEqual(p.returncode, 0, p.stdout)
        self.assertAsPplacer(os.path.join(self.tmp, 'x.json'))
        p = self.gplace('-c', testsets.MINI, '-o', 'y.json', 'queries.fasta', 'q.fasta.gz')
        self.assertNotEqual(p.returncode, 0)
        self.assertIn('several alignments', p.stdout)
        p = self.gplace('-c', testsets.MINI, '--map-identity', 'queries.fasta')
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("pplacer's --map-identity: not supported", p.stdout)
        p = self.gplace('-c', testsets.MINI)
        self.assertEqual(p.returncode, 0, p.stdout)          # pplacer: a warning, exit status 0
        self.assertIn('No query alignment', p.stdout)
        self.assertFalse([f for f in os.listdir(self.tmp) if f.endswith(('.jplace', '.json')) and f != 'x.json'])


if __name__ == '__main__':
    unittest.main()
