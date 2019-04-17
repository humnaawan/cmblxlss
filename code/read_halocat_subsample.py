import numpy as np
import pandas as pd
import os

home= '/Users/humnaawan/cmb_related_data/'

halocat_path= '%s/skyhalo_nres12r000.halo'%home
outDir= '%s/interm'%home

debug= False
keep= ['ID', 'Mvir', 'z_halo', 'theta_i', 'phi_i']
Mmin, Mmax= 5e12, 8e13

# ------------------------------------------------------------------------------------------------------
if debug:
        printProgress= True
        nSample= 500
        chunkSize= 10000
        outDir= '%s_debug'%(outDir)
else:
        printProgress= False
        nSample= 500000
        chunkSize= 100000
        outDir= '%s_%sgals'%(outDir, nSample)
        
outDir= '%s/'%(outDir)
if not os.path.exists(outDir): os.makedirs(outDir)
    
print('Outdir: %s'%outDir)

# read in the header and choose the columns that are readable.
header= pd.read_table(halocat_path, nrows=4, header=None, )
colNames= [(f) for f in header[0][0].split(' ') if (f!='' and not f.__contains__('#') and \
                                                not f.__contains__('rock'))]

    
########################################################################################################
def read_Takahashi_HaloCat(path, colNames, nrows= 30, skiprows= 0, printProgress= True):
        # skiprows input should not include the header
        nCol= len(colNames)
        data= pd.read_table(path, nrows= nrows, header=None, skiprows= 4+skiprows)
        for i, row in enumerate(data[0][:]):
                if (i==0): 
                        dataArr= np.array([float(f) for f in row.split(' ') \
                                           if (f!='' and not f.__contains__('list'))])
                else:
                        dataArr= np.vstack([dataArr, np.array([float(f) for f in row.split(' ') \
                                                               if (f!='' and not f.__contains__('list'))])])
        if printProgress:
                print('%s rows read in.'%len(dataArr))
                print('%s columns read in:\n%s\n'%(nCol, colNames))
                
        return pd.DataFrame(dataArr, columns= colNames
                           )

########################################################################################################                
goodN= 0
i= 0
while (goodN<nSample):
        data= read_Takahashi_HaloCat(path= halocat_path, colNames= colNames,
                                     nrows= chunkSize,
                                     skiprows= i*chunkSize,
                                     printProgress= printProgress)
    
        # find all the galaxies with parentID= -1 (central gals).
        # sample cuts: Mvir between 5e12 to 8e13 for now.
        subData= data[data.parent_ID == -1]
        subData.reset_index(drop=True)
        nKeep= len(subData)
        if printProgress:
                print('Dropped %s rows with parentID!= -1'%(len(data)-nKeep) )
                print('Final nRows: %s\n'%nKeep)

        # implement mass cuts
        subData= subData[(subData.Mvir> Mmin) & (subData.Mvir< Mmax)]
        subData.reset_index(drop=True)

        if printProgress:
                print('Dropped %s rows with Mvir not between 5e12 to 8e13'%(nKeep-len(subData)))
                print('Final nRows: %s'%len(subData))

        goodN+= len(subData)

        filename= '%s/output_%s.csv'%(outDir, i)
        subData.to_csv(filename, columns= keep)

        if printProgress: print('\n## Saved %s\n'%filename)
        i+=1

        if printProgress:
                print('%s galaxies; iteration # %s\n'%(goodN, i))
                print('------------------------------')

# ------------------------------------------------------------------------------------------------------
# combine the files now.
os.chdir(outDir)
filenames= [f for f in os.listdir(outDir) if not f.__contains__('gals')]

result= pd.read_csv(filenames[0])

for filename in filenames[1:]:
    if printProgress: print('Reading %s'%filename)
    result = pd.concat([result, pd.read_csv(filename)], ignore_index=True)

filename= '%s/finalOutput_%sgals.csv'%(outDir, len(result))
result.to_csv(filename, columns= keep, index=False)

print('\n## Saved %s\n'%filename)

