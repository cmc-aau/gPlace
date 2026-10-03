# gdb Python script: run with  gdb -q -batch -x <this file> --args pplacer --version  (see README.md)
# log1p evaluated inside GTDB-Tk's pplacer binary: every argument GSL 1.16's gsl_cdf_gamma_Pinv passes
# for discrete gamma rates with 21..300 categories (log1p(-P), P = i/n > 0.95), plus random and special
# inputs. Writes log1p_inputs.txt and log1p_stock_vals.txt (IEEE bits, hex) in the current directory.
import gdb, random
random.seed(19)
xs = []
for n in range(21, 301):
    xs += [-(i / n) for i in range(1, n) if i / n > 0.95]
for _ in range(600):
    xs.append(-(10 ** random.uniform(-20, 0)) * random.random())        # (-1, 0)
for _ in range(600):
    xs.append(10 ** random.uniform(-20, 10) * random.random())          # (0, 1e10)
for _ in range(200):
    xs.append(-1 + 10 ** random.uniform(-16, -1))                       # near -1
xs += [0.0, -0.0, 1e-300, -1e-300, 5e-324, 2.0 ** -29, -2.0 ** -29, 2.0 ** -54, 0.41422, -0.2929, 1.0, 1e300,
       -1.0, -2.0, float('inf'), float('-inf'), float('nan')]
# arguments and results are passed as raw bits, so signed zeros and NaN signs survive
import struct
gdb.execute("set pagination off"); gdb.execute("set confirm off")
gdb.execute("break main"); gdb.execute("run")
dbl = gdb.lookup_type("double")
fn = gdb.parse_and_eval("(double(*)(double))log1p").dereference()
bits = lambda v: struct.unpack("<Q", struct.pack("<d", v))[0]
out = open("log1p_stock_vals.txt", "w"); inp = open("log1p_inputs.txt", "w")
for x in xs:
    arg = gdb.Value(struct.pack("<d", x), dbl)
    assert bits(float(arg)) == bits(x)
    v = float(fn(arg))
    inp.write("%016x\n" % bits(x)); out.write("%016x\n" % bits(v))
out.close(); inp.close(); gdb.execute("kill"); gdb.execute("quit")
