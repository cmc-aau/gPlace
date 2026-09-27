import unittest

from gplace.cli import get_main_parser


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


if __name__ == '__main__':
    unittest.main()
