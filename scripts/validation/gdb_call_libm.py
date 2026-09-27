# gdb Python script: run with  gdb -q -batch -x <this file> --args pplacer ...  (see README.md)
import gdb, random
random.seed(7)
cases = []
for _ in range(1500):
    a = random.uniform(-3, 3) * 10 ** random.uniform(-6, 1); b = random.uniform(-3, 3) * 10 ** random.uniform(-6, 1)
    cases.append(("h", a, b))
for _ in range(3000): cases.append(("e", -random.uniform(0, 60) * random.random() ** 2, 0.0))
for _ in range(3000): cases.append(("l", 10 ** random.uniform(-300, 1) * random.random(), 0.0))
gdb.execute("set pagination off"); gdb.execute("set confirm off")
gdb.execute("break main"); gdb.execute("run")
out = open("stock_vals.txt", "w"); inp = open("inputs.txt", "w")
for op, a, b in cases:
    if op == "h": e = "((double(*)(double,double))hypot)(%r,%r)" % (a, b)
    elif op == "e": e = "((double(*)(double))exp)(%r)" % a
    else: e = "((double(*)(double))log)(%r)" % a
    v = float(gdb.parse_and_eval(e))
    inp.write("%s %s %s\n" % (op, float(a).hex(), float(b).hex())); out.write(v.hex() + "\n")
out.close(); inp.close(); gdb.execute("kill"); gdb.execute("quit")
