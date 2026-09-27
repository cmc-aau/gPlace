# gdb Python script: run with  gdb -q -batch -x <this file> --args pplacer ...  (see README.md)
import gdb, os, struct
D = os.environ["DUMP_DIR"]; inf = gdb.selected_inferior
def gmat(addr):   # gsl_matrix {size1,size2,tda,data,block,owner}
    s1, s2, tda, data = struct.unpack("<qqqq", inf().read_memory(addr, 32).tobytes())
    return inf().read_memory(data, s1 * tda * 8).tobytes(), (s1, s2, tda)
def gvec(addr):   # gsl_vector {size,stride,data,...}
    n, stride, data = struct.unpack("<qqq", inf().read_memory(addr, 24).tobytes())
    return inf().read_memory(data, n * stride * 8).tobytes()
class Fin(gdb.FinishBreakpoint):
    def __init__(self, fr, ev, evec): super().__init__(fr, internal=True); self.ev, self.evec = ev, evec
    def stop(self):
        open(f"{D}/symmv_eval.bin", "wb").write(gvec(self.ev))
        open(f"{D}/symmv_evec.bin", "wb").write(gmat(self.evec)[0])
        gdb.post_event(lambda: gdb.execute("kill")); return True
class Symmv(gdb.Breakpoint):
    def stop(self):
        a = int(gdb.parse_and_eval("$rdi")) & (2**64-1)
        b, dims = gmat(a)
        if dims[0] != 20: return False
        open(f"{D}/symmv_in.bin", "wb").write(b); open(f"{D}/symmv_dims.txt", "w").write(str(dims))
        Fin(gdb.newest_frame(), int(gdb.parse_and_eval("$rsi")) & (2**64-1), int(gdb.parse_and_eval("$rdx")) & (2**64-1))
        return False
Symmv("gsl_eigen_symmv")
for c in ("set pagination off", "set confirm off", "run"): gdb.execute(c)
