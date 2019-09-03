import fbf

lat = fbf.read_obj('latmat.pkl')
print(lat)
corejunction_node = fbf.read_obj('corejunction_node.pkl')
cjs1, cjs2, outv_cjs1, outv_cjs2, labeled_mol = fbf.analyze_bone('bone1.xyz')
cnode = fbf.vtree2bltree(corejunction_node)
cnode1 = fbf.align_decompose_bltree2vtree(cnode, outv_cjs1)
cnode1 = fbf.append_sidechain_number(cnode1, '1')
cnode1.symbol = 'C_cjunction1'
cnode2 = fbf.align_decompose_bltree2vtree(cnode, outv_cjs2)
cnode2 = fbf.append_sidechain_number(cnode2, '2')
cnode2.symbol = 'C_cjunction2'

# import anytree
# #print(cnode1.children)
# #print(cnode2.children)
# print(anytree.RenderTree(cnode1, style=anytree.AsciiStyle()))
# print(anytree.RenderTree(cnode2, style=anytree.AsciiStyle()))



from timeit import default_timer as timer
start = timer()
survival = fbf.grow([labeled_mol], 30, cnode1, cnode2, lat, 'inversion')
for i in range(len(survival)):
    survival[i].to('xyz', 'sur_' + str(i) + '.xyz')
    fbf.pack_mol2pbc(survival[i], lat).to('cif', 'sur_'+str(i)+'.cif')
end = timer()
print(end - start)
