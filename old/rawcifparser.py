from pymatgen.io.cif import CifFile


class RawCifParser:
    """
    clean/parse the raw cif file, write files:
        config0.cif     cif with no disorder
        config1.cif ...
    """
    def __init__(self, cifname):
        self.raw_cifname = cifname
        self.raw_d = list(CifFile.from_file(self.raw_cifname).data.items())[0][1].data
        self.configds = self.config_split()
        self.config_number = len(self.configds)
        self.config_files = self.write_all()

    def config_split(self):
        raw_d = self.raw_d
        config_ds = []
        if '_atom_site_disorder_group' in raw_d.keys():
            dg_labels = sorted([line for line in list(set([s.strip('-') for s in raw_d['_atom_site_disorder_group']]))
                                if line != u'.'])
            if len(dg_labels) > 0:
                for l in dg_labels:
                    config_d = dict(
                        name='data_config' + l,
                        a=raw_d['_cell_length_a'],
                        b=raw_d['_cell_length_b'],
                        c=raw_d['_cell_length_c'],
                        A=raw_d['_cell_angle_alpha'],
                        B=raw_d['_cell_angle_beta'],
                        C=raw_d['_cell_angle_gamma'],
                    )
                    if '_space_group_symop_operation_xyz' in raw_d.keys():
                        config_d['symop'] = raw_d['_space_group_symop_operation_xyz']
                    elif '_symmetry_equiv_pos_as_xyz' in raw_d.keys():
                        config_d['symop'] = raw_d['_symmetry_equiv_pos_as_xyz']
                    sites = []
                    for j in range(len(raw_d['_atom_site_label'])):
                        if raw_d['_atom_site_disorder_group'][j].strip('-') == l \
                                or raw_d['_atom_site_disorder_group'][j] == u'.':
                            site = [
                                raw_d['_atom_site_label'][j],
                                raw_d['_atom_site_type_symbol'][j],
                                raw_d['_atom_site_fract_x'][j],
                                raw_d['_atom_site_fract_y'][j],
                                raw_d['_atom_site_fract_z'][j],
                                '1.0',
                                raw_d['_atom_site_disorder_group'][j]
                            ]
                            sites.append(site)
                    config_d['msites'] = sites
                    config_ds.append(config_d)
            else:
                config_d = dict(
                    name='data_config0',
                    a=raw_d['_cell_length_a'],
                    b=raw_d['_cell_length_b'],
                    c=raw_d['_cell_length_c'],
                    A=raw_d['_cell_angle_alpha'],
                    B=raw_d['_cell_angle_beta'],
                    C=raw_d['_cell_angle_gamma'],
                )
                if '_space_group_symop_operation_xyz' in raw_d.keys():
                    config_d['symop'] = raw_d['_space_group_symop_operation_xyz']
                elif '_symmetry_equiv_pos_as_xyz' in raw_d.keys():
                    config_d['symop'] = raw_d['_symmetry_equiv_pos_as_xyz']
                sites = []
                for j in range(len(raw_d['_atom_site_label'])):
                    site = [
                        raw_d['_atom_site_label'][j],
                        raw_d['_atom_site_type_symbol'][j],
                        raw_d['_atom_site_fract_x'][j],
                        raw_d['_atom_site_fract_y'][j],
                        raw_d['_atom_site_fract_z'][j],
                        '1.0',
                        '.'
                    ]
                    sites.append(site)
                config_d['msites'] = sites
                config_ds.append(config_d)

        else:
            config_d = dict(
                name='data_config0',
                a=raw_d['_cell_length_a'],
                b=raw_d['_cell_length_b'],
                c=raw_d['_cell_length_c'],
                A=raw_d['_cell_angle_alpha'],
                B=raw_d['_cell_angle_beta'],
                C=raw_d['_cell_angle_gamma'],
            )
            if '_space_group_symop_operation_xyz' in raw_d.keys():
                config_d['symop'] = raw_d['_space_group_symop_operation_xyz']
            elif '_symmetry_equiv_pos_as_xyz' in raw_d.keys():
                config_d['symop'] = raw_d['_symmetry_equiv_pos_as_xyz']
            sites = []
            for j in range(len(raw_d['_atom_site_label'])):
                site = [
                    raw_d['_atom_site_label'][j],
                    raw_d['_atom_site_type_symbol'][j],
                    raw_d['_atom_site_fract_x'][j],
                    raw_d['_atom_site_fract_y'][j],
                    raw_d['_atom_site_fract_z'][j],
                    '1.0',
                    '.'
                ]
                sites.append(site)
            config_d['msites'] = sites
            config_ds.append(config_d)
        return config_ds

    def write_all(self):
        fns = []
        # fns.append('clean.cif')
        for i in range(self.config_number):
            self.write_config_dict(self.configds[i], 'config' + str(i) + '.cif')
            fns.append('config' + str(i) + '.cif')
        return fns

    @staticmethod
    def write_config_dict(d, cifn):
        """
        write config_dict into minimum readable cif file
        :return:
        """
        with open(cifn, 'w') as f:
            f.write(d['name'] + '\n')
            f.write(
                'loop_\n _symmetry_equiv_pos_as_xyz\n' + '\n'.join(["'" + xyz + "'" for xyz in d['symop']]) + '\n\n')
            f.write('_cell_length_a\t' + d['a'] + '\n')
            f.write('_cell_length_b\t' + d['b'] + '\n')
            f.write('_cell_length_c\t' + d['c'] + '\n')
            f.write('_cell_angle_alpha\t' + d['A'] + '\n')
            f.write('_cell_angle_beta\t' + d['B'] + '\n')
            f.write('_cell_angle_gamma\t' + d['C'] + '\n')
            f.write('loop_\n _atom_site_label\n _atom_site_type_symbol\n _atom_site_fract_x\n _atom_site_fract_y\n '
                    '_atom_site_fract_z\n _atom_site_occupancy\n _atom_site_disorder_group\n')
            for s in d['msites']:
                line = '\t'.join(s)
                f.write(line + '\n')
