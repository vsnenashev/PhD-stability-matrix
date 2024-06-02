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
