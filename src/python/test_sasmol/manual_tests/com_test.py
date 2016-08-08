import os
import sasmol.sasmol as sasmol
import sasmol.calculate as calculate

pdbDataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','pdb_common')+os.path.sep

m = sasmol.SasMol(0)

try:
    m.read_pdb("hiv1_gag.pdb")
except:
    m.read_pdb(pdbDataPath+"hiv1_gag.pdb")

frame = 0

print 'com = ',m.calculate_center_of_mass(frame)

m.setCom = m.calculate_center_of_mass(frame)

print 'ugly = ',calculate.Calculate.calculate_center_of_mass(m, frame)

calc_com = calculate.Calculate.calculate_center_of_mass

print 'prettier  = ',calc_com(m, frame)



