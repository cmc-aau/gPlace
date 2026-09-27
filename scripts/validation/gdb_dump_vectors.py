# gdb Python script: run with  gdb -q -batch -x <this file> --args pplacer ...  (see README.md)
# gdb Python: dump bit-exact arrays from stock pplacer.
#  dediag_{k}_{u,lambda,uit}.bin : first 3 calls of dediagonalize (model eigvecs + exp(rate*bl*l))
#  ll3_{statd,x,y,z}.bin         : first mat_log_like3_c call (prox_v, dist_v, pend_v evolved glvs)
import gdb, os, struct
D = os.environ["DUMP_DIR"]; inf = lambda: gdb.selected_inferior()
def ba(reg):
    v = int(gdb.parse_and_eval(reg)) & 0xFFFFFFFFFFFFFFFF
    mem = inf().read_memory(v + 8, 48).tobytes()
    data, nd, flags, proxy, d0, d1 = struct.unpack("<qqqqqq", mem)
    n = d0 * (d1 if nd >= 2 else 1)
    return inf().read_memory(data, n * 8).tobytes(), (nd, d0, d1)
class Dediag(gdb.Breakpoint):
    k = 0
    def stop(self):
        if Dediag.k < 3:
            for nm, r in (("u", "$rsi"), ("lambda", "$rdx"), ("uit", "$rcx")):
                b, dims = ba(r); open(f"{D}/dediag_{Dediag.k}_{nm}.bin", "wb").write(b)
            Dediag.k += 1
        return False
class LL3(gdb.Breakpoint):
    done = False
    def stop(self):
        if not LL3.done:
            for nm, r in (("statd", "$rdi"), ("x", "$rsi"), ("y", "$rdx"), ("z", "$rcx")):
                b, dims = ba(r); open(f"{D}/ll3_{nm}.bin", "wb").write(b)
                open(f"{D}/ll3_dims.txt", "a").write(f"{nm} {dims}\n")
            LL3.done = True
        return False
Dediag("dediagonalize"); LL3("mat_log_like3_c")
for c in ("set pagination off", "set follow-fork-mode child", "set detach-on-fork on", "set confirm off", "run", "quit"):
    gdb.execute(c)
