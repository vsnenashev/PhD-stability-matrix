#ins = (inv(k_global).dot(g_global)).toarray()
#print(la.eigvals(ins))
#insa = (inv(g_global).dot(k_global)).toarray()
#print(la.eigvals(insa))
#insi = k_global.dot(inv(g_global)).toarray()
#print(la.eigvals(insi))
#print(lam_i)
import pandas as pd
import ezdxf as ez
from ezdxf.addons import iterdxf

doc = ez.new()
doc1 = iterdxf.opendxf('big.dxf')
msp = doc.modelspace()
for entity in doc.modelspace():
    doc.write(entity)


doc.saveas("Cut.dxf")
doc1.close()
