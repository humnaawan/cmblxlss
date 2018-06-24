import readcol
import pandas as pd
import time

path= "/global/cscratch1/sd/awan/takahashi_data/skyhalo_nres12r000.halo"
tag= 'attempt1'

outDir= "/global/homes/a/awan/SO/CMBL-LSS-CrossCorr/data"

keep= ['ID', 'Mvir', 'z_halo', 'theta_i', 'phi_i']
Mmin, Mmax= 5e12, 8e13

# ------------------------------------------------------------------------------------------------------
startTime= time.time()
readme= ''

# read in the data
colNames, cols= readcol.readcol(path, twod= False, names= True)
readme+= 'Columns read in: %s\n'%(colNames)
#readme+= 'Data shape: %s\n'%(np.shape(cols),)

colNames= colNames[1:]
# put into dataframe; easier to handle
data= {}
for ith, col in enumerate(colNames):
    data[col]= cols[ith]
data= pd.DataFrame(data)

readme+= 'Columns in the dataframe: %s\n'%(data.keys(),)

# ----------------------------------------
# apply cuts
# ----------------------------------------
# find all the galaxies with parentID= -1 (central gals).
subData= data[data.parent_ID == -1]
subData.reset_index(drop=True)
nKeep= len(subData)

readme+= 'Dropped %s rows with parentID!= -1\n'%(len(data)-nKeep)
readme+= 'Final nRows: %s\n'%nKeep

# ----------------------------------------
# implement mass cuts
subData= subData[(subData.Mvir> Mmin) & (subData.Mvir< Mmax)]
subData.reset_index(drop=True)

readme+= 'Dropped %s rows with Mvir not between %s to %s\n'%(nKeep-len(subData), Mmin, Mmax)
readme+= 'Final nRows: %s\n'%len(subData)

# ----------------------------------------
# save file
filename= '%s/haloCat_readcol_%sgals.csv'%(outDir, len(subData))
subData.to_csv(filename, columns= keep, index=False)

readme+= '\n## Saved %s in %s\n'%(filename, outDir)
readme+= '\nTime taken: %s min'%((time.time()-startTime)/60.)

readme_file= open('%s/readme_%s.txt'%(outDir, tag), 'a')
readme_file.write(readme)
readme_file.close()


