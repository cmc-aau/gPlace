# gdb Python script: run with  gdb -q -batch -x <this file> --args pplacer ...  (see README.md)
# gdb Python: log every three-taxon evaluation of stock pplacer (bit-exact, hex floats).
#  O            start of Three_tax.optimize
#  P <x>        set_pend_bl argument
#  D <x>        set_dist_bl argument
#  L <x>        mat_log_like3_c return value (before exponent term)
import gdb, os
out = open(os.environ["TRACE_OUT"], "w", buffering=1)

def boxed(reg):
    a = int(gdb.parse_and_eval(reg)) & 0xFFFFFFFFFFFFFFFF
    return float(gdb.parse_and_eval("*(double*)%d" % a)).hex()

class Arg(gdb.Breakpoint):
    def __init__(self, spec, tag):
        super().__init__(spec); self.tag = tag
    def stop(self):
        out.write(self.tag if self.tag == "O" else "%s %s" % (self.tag, boxed("$rbx")))
        out.write("\n"); return False

class RetLL(gdb.FinishBreakpoint):
    def stop(self):
        out.write("L %s\n" % boxed("$rax")); return False
    def out_of_scope(self):
        pass

class LL(gdb.Breakpoint):
    def stop(self):
        RetLL(gdb.newest_frame(), internal=True); return False

Arg("camlThree_tax__optimize_10631", "O")
Arg("camlThree_tax__set_pend_bl_10608", "P")
Arg("camlThree_tax__set_dist_bl_10611", "D")
LL("mat_log_like3_c")
gdb.execute("set pagination off")
gdb.execute("set follow-fork-mode child")
gdb.execute("set detach-on-fork on")
gdb.execute("set confirm off")
gdb.execute("run")
gdb.execute("quit")
