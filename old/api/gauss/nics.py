

def nics_sigma_structure(self, normaldirection=0):

    lgfr = self.largest_fr()
    geocs = [r.geoc for r in lgfr]  # nx3
    dmat = squareform(pdist(np.array(geocs)))
    maxi, maxj = np.unravel_index(dmat.argmax(), dmat.shape)
    fr = sorted(lgfr, key=lambda r: np.linalg.norm(r.geoc - lgfr[maxi].geoc))
    available_normals = fr[0].normals
    normals = [r.normal_along(available_normals[normaldirection]) for r in fr]

    frsites = []
    for r in lgfr:
        for s in r.sites:
            if s not in frsites:
                frsites.append(s)

    sigma_sites = deepcopy(self.sites)

    terminated_sites_id = []
    added_hydrogens = []

    for i in range(len(lgfr)):
        for s in lgfr[i].sites:
            if s.insaturation == 1 and s.id not in terminated_sites_id:
                terminated_sites_id.append(s.id)
                hs = Site('H', 1.0*normals[i]+s.coords)
                added_hydrogens.append(hs)
            for nb in s.nbs:
                if nb.id not in terminated_sites_id and nb.insaturation == 1 and nb not in frsites and all([nb not in ring.sites for ring in self.rings]) :
                    terminated_sites_id.append(nb.id)
                    hs = Site('H', 1.0*normals[i]+nb.coords)
                    added_hydrogens.append(hs)
    return Sitelist(sigma_sites + added_hydrogens)

def nics_line_scan_path(self, step_size, nrings, height=1.7, normaldirection=0):
    lgfr = self.largest_fr()[:nrings]
    geocs = [r.geoc for r in lgfr]  # nx3
    dmat = squareform(pdist(np.array(geocs)))
    maxi, maxj = np.unravel_index(dmat.argmax(), dmat.shape)
    fr = sorted(lgfr, key=lambda r: np.linalg.norm(r.geoc - lgfr[maxi].geoc))
    available_normals = fr[0].normals
    normals = [r.normal_along(available_normals[normaldirection]) for r in fr]

    pts = []
    ring_idx = []
    xnumbers = []
    xticks = [0]
    cross_point = np.zeros(3)
    for i in range(len(fr)-1):
        ring = fr[i]
        nb_ring = fr[i+1]
        segend1 = ring.geoc
        segend2 = nb_ring.geoc
        # nb_ring center -- cross_pt -- ring center -- prev_cross_pt -- prev ring center
        if i == 0:
            bond_centers = sorted([b.center for b in ring.bonds], key=lambda bc: angle_btw(segend2-segend1, bc-segend1))
            cross_point, start_point = sorted([bond_centers[0], bond_centers[-1]], key=lambda c: np.linalg.norm(c-0.5*(segend2+segend1)))
            subpath = genlinepts(start_point, ring.geoc, step_size)
            xnumbers += [np.linalg.norm(p-subpath[0])+xticks[-1] for p in subpath]
            xticks += [xticks[-1] + np.linalg.norm(ring.geoc - start_point)]
            subpath = [p + normals[i]*height for p in subpath]
            pts += subpath
            ring_idx += [i]*len(subpath)


            subpath = genlinepts(ring.geoc, cross_point, step_size)
            xnumbers += [np.linalg.norm(p-subpath[0])+xticks[-1] for p in subpath]
            xticks += [xticks[-1] + np.linalg.norm(ring.geoc - cross_point)]
            subpath = [p + normals[i]*height for p in subpath]
            pts += subpath
            ring_idx += [i]*len(subpath)
        elif i != 0 and i != len(fr)-2:
            prev_cross_point = cross_point
            bond_centers = sorted([b.center for b in ring.bonds], key=lambda bc: angle_btw(segend2-segend1, bc-segend1))
            cross_point, start_point = sorted([bond_centers[0], bond_centers[-1]], key=lambda c: np.linalg.norm(c-0.5*(segend2+segend1)))
            subpath = genlinepts(prev_cross_point, ring.geoc, step_size)
            xnumbers += [np.linalg.norm(p-subpath[0])+xticks[-1] for p in subpath]
            xticks += [xticks[-1] + np.linalg.norm(prev_cross_point - ring.geoc)]
            subpath = [p + normals[i]*height for p in subpath]
            pts += subpath
            ring_idx += [i]*len(subpath)

            subpath = genlinepts(ring.geoc, cross_point, step_size)
            xnumbers += [np.linalg.norm(p-subpath[0])+xticks[-1] for p in subpath]
            xticks += [xticks[-1] + np.linalg.norm(cross_point - ring.geoc)]
            subpath = [p + normals[i]*height for p in subpath]
            pts += subpath
            ring_idx += [i]*len(subpath)
        elif i == len(fr)-2:
            prev_cross_point = cross_point
            bond_centers = sorted([b.center for b in ring.bonds], key=lambda bc: angle_btw(segend2-segend1, bc-segend1))
            cross_point, start_point = sorted([bond_centers[0], bond_centers[-1]], key=lambda c: np.linalg.norm(c-0.5*(segend2+segend1)))
            subpath = genlinepts(prev_cross_point, ring.geoc, step_size)
            xnumbers += [np.linalg.norm(p-subpath[0])+xticks[-1] for p in subpath]
            xticks += [xticks[-1] + np.linalg.norm(prev_cross_point - ring.geoc)]
            subpath = [p + normals[i]*height for p in subpath]
            pts += subpath
            ring_idx += [i]*len(subpath)

            subpath = genlinepts(ring.geoc, cross_point, step_size)
            xnumbers += [np.linalg.norm(p-subpath[0])+xticks[-1] for p in subpath]
            xticks += [xticks[-1] + np.linalg.norm(cross_point - ring.geoc)]
            subpath = [p + normals[i]*height for p in subpath]
            pts += subpath
            ring_idx += [i]*len(subpath)

            subpath = genlinepts(cross_point, nb_ring.geoc, step_size)
            xnumbers += [np.linalg.norm(p-subpath[0])+xticks[-1] for p in subpath]
            xticks += [xticks[-1] + np.linalg.norm(cross_point - nb_ring.geoc)]
            subpath = [p + normals[i+1]*height for p in subpath]
            pts += subpath
            ring_idx += [i+1]*len(subpath)

            bond_centers = sorted([b.center for b in nb_ring.bonds], key=lambda bc: angle_btw(segend1-segend2, bc-segend2))
            cross_point, end_point = sorted([bond_centers[0], bond_centers[-1]], key=lambda c: np.linalg.norm(c-0.5*(segend2+segend1)))
            subpath = genlinepts(nb_ring.geoc, end_point, step_size)
            xnumbers += [np.linalg.norm(p-subpath[0])+xticks[-1] for p in subpath]
            xticks += [xticks[-1] + np.linalg.norm(end_point - nb_ring.geoc)]
            subpath = [p + normals[i+1]*height for p in subpath]
            pts += subpath
            ring_idx += [i+1]*len(subpath)
    return pts, ring_idx, xnumbers, xticks
