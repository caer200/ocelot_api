#  This file is (c) 2012, Noel O'Boyle
#  All rights reserved.
#
#  This file is available under the BSD license
#  available from http://www.opensource.org/licenses/bsd-license.php

def getCanonicalLabels(inchi, aux):
    # In the InChI, use the connection layer from /R if available
    recon_start = inchi.find("/r")
    if recon_start == -1:
        split_inchi = inchi.split("/")
    else:
        split_inchi = [""] + inchi[recon_start+1:].split("/")
        recon_start_b = aux.find("/R:")
        assert recon_start_b >= 0
        aux = aux[recon_start_b+1:]

    split_aux = aux.split("/")
    assert split_aux[2][0] == "N"
    
    # In the Auxiliary Info, use the reconnected canonical labels if present
    # - otherwise use the normal labels.
    # Either way, then adjust using any following /f section.
        
    canlabels = [map(int, x.split(",")) for x in split_aux[2][2:].split(";")]

    fixedH_start = aux.find("/F")
    if fixedH_start >= 0:
        broken = aux[fixedH_start+3:].split("/")[0]
        mols = broken.split(";")
        new_canlabels = []
        tot = 0
        for mol in mols:
            if mol == "m":
                mol = "1m"
            if mol.endswith("m"):
                mult = int(mol[:-1])
                new_canlabels += canlabels[tot:tot+mult]
                tot += mult
            else:
                new_canlabels.append(map(int, mol.split(",")))
                tot += 1
        canlabels = new_canlabels

    # Flatten canlabels
    canlabels = [item for sublist in canlabels for item in sublist]

    return canlabels

if __name__ == "__main__":
    # The following examples are from Table 1 in the paper
    print "Table One"
    examplesTableOne = """InChI=1/C2H6O/c1-2-3/h3H,2H2,1H3	AuxInfo=1/0/N:3,2,1/rA:3OCC/rB:s1;s2;/rC:;;;
InChI=1/C2H2BrClO/c3-2(5)1-4/h1H2	AuxInfo=1/0/N:2,3,5,1,4/rA:5ClCCOBr/rB:s1;s2;d3;s3;/rC:;;;;;
InChI=1/2CH4/h2*1H4	AuxInfo=1/0/N:1;2/rA:2CC/rB:;/rC:;;
InChI=1/C2H6.CH4/c1-2;/h1-2H3;1H4	AuxInfo=1/0/N:2,3;1/E:(1,2);/rA:3CCC/rB:;s2;/rC:;;;
InChI=1/C2H6/c1-2/h1-2H3	AuxInfo=1/0/N:1,3/E:(1,2)/rA:3CHC/rB:s1;s1;/rC:;;;
InChI=1S/CH4O/c1-2/h2H,1H3/i1D	AuxInfo=1/0/N:1,3/rA:3CH.i2O/rB:s1;s1;/rC:;;;"""

    for line in examplesTableOne.split("\n"):
        inchi, aux = line.split()
        print getCanonicalLabels(inchi, aux)

    # The following examples are from Table 2 in the paper
    print "\nTable Two"
    examplesTableTwo = """InChI=1/C2H2O4/c3-1(4)2(5)6/h(H,3,4)(H,5,6)/p-1	AuxInfo=1/1/N:1,4,2,3,5,6/E:(1,2)(3,4,5,6)/gE:(1,2)/rA:6COO-COO/rB:d1;s1;s1;d4;s4;/rC:;;;;;;
InChI=1/C2H2O4/c3-1(4)2(5)6/h(H,3,4)(H,5,6)/p-1	AuxInfo=1/1/N:1,4,2,3,5,6/E:(1,2)(3,4,5,6)/gE:(1,2)/rA:6COOCOO-/rB:d1;s1;s1;d4;s4;/rC:;;;;;;
InChI=1/C2H2O4/c3-1(4)2(5)6/h(H,3,4)(H,5,6)/p-1/fC2HO4/h3H/q-1	AuxInfo=1/1/N:1,4,2,3,5,6/E:(1,2)(3,4,5,6)/gE:(1,2)/F:4,1,6,5,2,3/E:(5,6)/rA:6COO-COO/rB:d1;s1;s1;d4;s4;/rC:;;;;;;
InChI=1/C2H2O4/c3-1(4)2(5)6/h(H,3,4)(H,5,6)/p-1/fC2HO4/h3H/q-1	AuxInfo=1/1/N:1,4,2,3,5,6/E:(1,2)(3,4,5,6)/gE:(1,2)/F:1,4,3,2,5,6/E:(5,6)/rA:6COOCOO-/rB:d1;s1;s1;d4;s4;/rC:;;;;;;"""
    for line in examplesTableTwo.split("\n"):
        inchi, aux = line.split()
        print getCanonicalLabels(inchi, aux)

    # Here's a reconnected metal example
    print "\nReconnected Metal"
    line = """InChI=1/CH2O2.ClH.2H2N.Pt/c2-1-3;;;;/h1H,(H,2,3);1H;2*1H2;/q;;2*-1;+4/p-2/fCHO2.Cl.2H2N.Pt/h;1h;;;/q2*-1;3m/rCH5ClN2O2Pt/c2-7(3,4)6-1-5/h1H,3-4H2
AuxInfo=1/1/N:3,2,4;7;5;6;1/E:(2,3);;;;/F:5m/E:m;;;;/CRV:;;2*1-1;/rA:7PtOCONNCl/rB:s1;s2;d3;s1;s1;s1;/rC:;;;;;;;/R:/0/N:3,7,5,6,4,2,1/E:(3,4)"""
    inchi, aux = line.split()
    print getCanonicalLabels(inchi, aux)
