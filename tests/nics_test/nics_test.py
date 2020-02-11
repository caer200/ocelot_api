from ocelot.schema.conformer import MolConformer
from ocelot.task.nics import NICSjob
"""
this will wirte sigma-0.gjf and total-0.gjf
then print zz values from sigma-0.log and total-0.log
"""

def ws(s, fn):
    with open(fn, 'w') as f:
        f.write(s)

def rf(fn):
    with open(fn, 'r') as f:
        s = f.read()
        return s


omol = MolConformer.from_file('CNBr_s0_optgeo.gjf')

nicsjob = NICSjob(omol, nmrscheme='GIAO')

nrings = len(omol.backbone.rings)

functional = 'LC-wpbe'
basis_set='6-311++G(d,p)'
route_parameters={
    'iop':{'3/108':'0156600000', '3/107':'0156600000'},
    'NMR': 'GIAO',
    'SCF': {'qc':'', 'maxcycle':'512'}
                  }
link0_parameters={
    '%rwf': 'Br_s0-0-nicsxy-backbone-sigma.rwf',
    '%nosave': '',
    '%chk': 'Br_s0-0-nicsxy-backbone-sigma.chk',
}

# C-wPBE/6-311++g(d,p) IOP(3/107=0156600000) IOP(3/108=0156600000) NMR=GIAO SCF=(qc,MAXCYCLE=512)


sigma_outs_list, total_outs_list, xticks, xnumbers, pt_idx, pts = nicsjob.gen_input(
    step_size=0.4, nrings=nrings, maxnbq=40, height=1.7, normaldirection=0, charge=0, spin_multiplicity=1,
    title=None, functional=functional, basis_set=basis_set, route_parameters=route_parameters,
    input_parameters={},
    link0_parameters=link0_parameters, dieze_tag='#P', gen_basis=None
)

isigma = 0
itotal = 0
nsigma = len(sigma_outs_list)

for i in range(nsigma):
    sigma_out = sigma_outs_list[i]
    total_out = total_outs_list[i]
    ws(sigma_out, 'sigma-{}.gjf'.format(isigma))
    ws(total_out, 'total-{}.gjf'.format(itotal))
    isigma += 1
    itotal += 1

sigma_logs = rf('sigma-0.log')
total_logs = rf('total-0.log')
sigma_tensor_zzs = nicsjob.get_zzmstensor(sigma_logs)
total_tensor_zzs = nicsjob.get_zzmstensor(total_logs)
print(sigma_tensor_zzs)
print(total_tensor_zzs)

#
#
