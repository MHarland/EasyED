from tetrahedron import tetrahedron

tetrahedron.calcOccupation()
for key, val in tetrahedron.occupation.items():
    print key,': ',val
print 'N_tot =', tetrahedron.getTotalOccupation()
