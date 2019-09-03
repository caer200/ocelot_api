import fbf
import sys
fn = sys.argv[1]
latmat = fbf.get_pbc(fn)
print(latmat)
fbf.write_obj(latmat, 'latmat.pkl')
