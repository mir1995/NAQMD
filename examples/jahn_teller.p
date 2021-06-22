set term png  
set output "jahn_teller.png"

eplus(x,y, gamma) = gamma*(x**2 + y**2) + sqrt(x**2 + y**2)
eminus(x,y, gamma) = gamma*(x**2 + y**2) - sqrt(x**2 + y**2)

set hidden3d
set isosamples 30
#set view map
#set contour surface
set view 80,30

gamma = 0
splot [-0.2:0.2][-0.2:0.2] eplus(x, y, gamma),  eminus(x, y, gamma)

set output
