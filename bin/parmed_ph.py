import parmed as pmd
import sys

top=sys.argv[1]
gro=sys.argv[2]
name_top=sys.argv[3]
name_rst=sys.argv[4]


gmx_top = pmd.load_file(top, xyz=gro)

gmx_top.save(name_rst, format='rst7')
gmx_top.save(name_top, format='amber')


