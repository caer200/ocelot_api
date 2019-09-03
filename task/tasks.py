from api.routines.pbc import CIFparser, PBCparser
from fireworks.core.firework import FWAction, FiretaskBase
import os
from collections import OrderedDict
import abc
import warnings
from api.schema.config import Config
from api.schema.dimer import DimerCollection


class OCELOT_task:

    def __init__(self, name):
        self.name = name
        self.result = OrderedDict()
        self.result['task_name'] = name
        self.run()

    @abc.abstractmethod
    def run(self, *args, **kwargs):
        return


class task_cifparse(OCELOT_task):
    """
    read cif string to config, clean up disorder
    """
    def run(self, cifstring=''):
        if cifstring == '':
            raise TypeError('cifstring is empty')
        try:
            cp = CIFparser.from_cifstring(cifstring)
        except:
            warnings.warn('cannot parse cifstring {} in task {}'.format(cifstring, self.name))
            return None
        self.result['cif_parser'] = cp.as_dict()


class task_configparse(OCELOT_task):
    """
    read config cif (no disorder) into molecule, bone, hbone, sidechain and dimers
    """
    def run(self, configcifstr):
        config = Config.from_cifstring(configcifstr)
        self.result['config'] = config.as_dict()

        for i in range(config.z):
            dimers_ref_i = config.dimers[i].flatten()
            dimers_ref_i = sorted(dimers_ref_i, key=lambda x: x.vslipnorm)[:27*config.z]
            dcstr = DimerCollection([dimer for dimer in dimers_ref_i if not dimer.is_identical]).get_xyz_string()





        hbone_config, hbone_config_pstructure, terminated_backbones= config.get_bone_config()

        bc_dimers, bc_dimers_transv_fcs = self.boneconfig.get_dimers()
        report = dict()










# class Task_CIFparse(FiretaskBase):
#     _fw_name = "CIFparse"
#
#     def run_task(self, fw_spec):
#         input_array = fw_spec['input_array']
#         m_sum = sum(input_array)
#
#         print("The sum of {} is: {}".format(input_array, m_sum))
#
#         return FWAction(stored_data={'sum': m_sum},
#                         mod_spec=[{'_push': {'input_array': m_sum}}])
