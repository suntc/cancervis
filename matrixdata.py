'''
Created on 3016-8-38

@author: sun
'''

import scipy.io
import numpy as np
import h5py

def mat_to_json(path, curvepath, ch=1):
	try:
		mat = scipy.io.loadmat(path)
	except NotImplementedError as e:
		mat = h5py.File(path)
	curvemat = scipy.io.loadmat(curvepath)
	arr_dict = {}
	for key in mat:
		if key not in ['__header__', '__globals__', '__version__']:
			arr_dict[key] = mat.get(key)
			arr_dict[key] = np.array(arr_dict[key])
			arr_dict[key] = np.where(np.isnan(arr_dict[key]), -1, arr_dict[key])
	## use channel 1 as an example
	yr, xr = np.nonzero(arr_dict['lifet_ch'+str(ch)])
	to_export = []
	minx = min(xr)
	miny = min(yr)
	shape = arr_dict['tumor'].shape
	for i in range(min(yr), max(yr)):## row
		for j in range(min(xr), max(xr)):
			to_export.append({'lifetime':arr_dict['lifet_ch'+str(ch)][i][j], 'intensity':arr_dict['int_ch'+str(ch)][i][j], 'tumor':arr_dict['tumor'][i][j], 'x':j-minx, 'y':i-miny})
	
	output = open("./matData_ch"+str(ch)+".js","wb")
	print >> output, "var matData_ch" + str(ch) + " = "
	print >> output, to_export
	print >> output, "var maxx_ch" + str(ch) + " = "
	print >> output, max(xr)-minx
	print >> output, "var maxy_ch" + str(ch) + " = "
	print >> output, max(yr)-miny
	output.close()
	#output = open("./curveData_ch"+str(ch)+".js","wb")
	## build dict
	curveDict = {}
	for i, row in enumerate(yr):
			col = xr[i]
			print row, col
			curveDict[str(row-miny) + ',' + str(col-minx)] = curvemat['C'][row,col,0:625].tolist()
	#output = open("./curveData_ch"+str(ch)+".json","wb")
	#print >> output, '{"curves":', curvemat['C'][miny:max(yr)+1,minx:max(xr)+1,0:625].tolist(), "}"
	#print >> output, curvemat['c1'][300:303,300:303,0:635].tolist()
	output = open("./curveDataDict_ch"+str(ch)+".js","wb")
	print >> output, ("var curveDataDict_ch"+str(ch)+" ="), curveDict
	output.close()

def fat_mat(path, fpath):
	fmat = scipy.io.loadmat(fpath)
	mat = scipy.io.loadmat(path)

	arr_dict = {}
	for key in mat:
		if key not in ['__header__', '__globals__', '__version__']:
			arr_dict[key] = mat.get(key)
			arr_dict[key] = np.array(arr_dict[key])
			arr_dict[key] = np.where(np.isnan(arr_dict[key]), -1, arr_dict[key])
	## use channel 1 as an example
	yr, xr = np.nonzero(arr_dict['lifet_ch1'])
	minx = min(xr)
	miny = min(yr)
	maxx = max(xr)
	maxy = max(yr)
	fat_mat = fmat['fat']
	output = open("./fatMat.js","wb")
	print >> output, ("var fatMat ="), fat_mat[miny:maxy+1,minx:maxx+1].tolist()
	output.close()
	
def save_curvemat(curvepath):
	curvemat = scipy.io.loadmat(curvepath)
	output = open('curvemat_ch1.pkl', 'wb')
	np.save(output, curvemat['c1'][:,:,0:635])
	output.close()
	
def curve_eudist(path, curvepath, ch=1):
	## load data
	mat = scipy.io.loadmat(path)
	curvemat = scipy.io.loadmat(curvepath)
	arr_dict = {}
	for key in mat:
		if key not in ['__header__', '__globals__', '__version__']:
			arr_dict[key] = mat.get(key)
			arr_dict[key] = np.array(arr_dict[key])
			arr_dict[key] = np.where(np.isnan(arr_dict[key]), -1, arr_dict[key])
	row, col = np.nonzero(arr_dict['lifet_ch'+str(ch)])
	mincol = min(col)
	minrow = min(row)
	maxcol = max(col)
	maxrow = max(row)
	## compute pairwise eudist
	from scipy.spatial.distance import cdist
	curve_in_use = curvemat['C'][minrow:maxrow+1,mincol:maxcol+1,0:635] ## need plus 1 because maxrow should be included
	trow, tcol = np.nonzero(arr_dict['tumor'])
	tcol = tcol - mincol ## convert to selected range
	trow = trow - minrow
	tumor_coords = np.dstack((trow,tcol))[0]
	tumor_curves = curve_in_use[trow,tcol]
	nonzerocol = col - mincol
	nonzerorow = row - minrow
	nonzero_coords = np.dstack((nonzerorow,nonzerocol))[0]
	nonzero_coords = np.array(nonzero_coords, 'i')
	nonzero_curves = curve_in_use[nonzerorow,nonzerocol]
	tumor_curves_eudist = cdist(tumor_curves, nonzero_curves)
	nearest_index = np.argsort(tumor_curves_eudist)
	nearest_dict = {}
	top = 100
	for i, coord in enumerate(tumor_coords):
		nearest_dict[str(coord.tolist()[0])+','+str(coord.tolist()[1])] = nonzero_coords[ nearest_index[i][0:top] ].tolist()
	output = open("./curveDist"+str(ch)+".js","wb")
	print >> output, "var curveDist_ch" + str(ch) + " = "
	print >> output, nearest_dict
	output.close()
	'''
	## too slow
	curve_in_use = np.reshape(curve_in_use, [curve_in_use.shape[0]*curve_in_use.shape[1],curve_in_use.shape[3]])
	curve_eudist = cdist(curve_in_use,curve_in_use)
	output = open("./curveDist.js", "wb")
	print >> output, "var curveDist = "
	print >> output, curve_eudist.tolist()
	output.close()
	'''

def curve_eudist_all(path, curvepath):
  ## load data
	mat = scipy.io.loadmat(path)
	curvemat = scipy.io.loadmat(curvepath)
	arr_dict = {}
	for key in mat:
		if key not in ['__header__', '__globals__', '__version__']:
			arr_dict[key] = mat.get(key)
			arr_dict[key] = np.array(arr_dict[key])
			arr_dict[key] = np.where(np.isnan(arr_dict[key]), -1, arr_dict[key])
	row, col = np.nonzero(arr_dict['lifet_ch1'])
	mincol = min(col)
	minrow = min(row)
	maxcol = max(col)
	maxrow = max(row)
	## compute pairwise eudist
	from scipy.spatial.distance import cdist
	curve_in_use = curvemat['c1'][minrow:maxrow+1,mincol:maxcol+1,0:635] ## need plus 1 because maxrow should be included
	nonzerocol = col - mincol
	nonzerorow = row - minrow
	nonzero_coords = np.dstack((nonzerorow,nonzerocol))[0]
	nonzero_coords = np.array(nonzero_coords, 'i')
	nonzero_curves = curve_in_use[nonzerorow,nonzerocol]
	curves_eudist = cdist(nonzero_curves, nonzero_curves)
	import pickle
	print "dump"
	output = open("./curveDist_all.pkl","wb")
	pickle.dump(curves_eudist, output)
	print "end dump"
	output.close()


def neighbor_contrast(path, curvepath, ch=1):
	
	def get_bound_sum(coord, d):
		dist = 0
		count = 0
		r = coord[0]
		c = coord[1]
		cind = nonzero_coords_dict[tuple(coord)]
		maxr = indexmat.shape[0]
		maxc = indexmat.shape[1]
		tocon = [indexmat[max(r-d+1,0):r+d,max(c-d,0)], indexmat[max(r-d+1,0):r+d,min(c+d,maxc-1)]]
		if r-d >= 0:
			tocon.append(indexmat[r-d,max(c-d, 0):c+d+1])
		if r+d < maxr:
			tocon.append(indexmat[r+d,max(c-d, 0):c+d+1])
		#indexmat[max(r-d, 0),max(c-d, 0):c+d+1], indexmat[min(r+d, maxr),max(c-d, 0):c+d+1], indexmat[max(r-d+1,0):r+d,max(c-d,0)], indexmat[max(r-d+1,0):r+d,min(c+d,maxc)]
		bound_coord = np.concatenate(tocon)
		for item in bound_coord:
			key = tuple(item)
			if key in nonzero_coords_dict:
				count += 1
				dist += curves_eudist[cind, nonzero_coords_dict[key]]
		return dist, count
	
	mat = scipy.io.loadmat(path)
	curvemat = scipy.io.loadmat(curvepath)
	arr_dict = {}
	for key in mat:
		if key not in ['__header__', '__globals__', '__version__']:
			arr_dict[key] = mat.get(key)
			arr_dict[key] = np.array(arr_dict[key])
			arr_dict[key] = np.where(np.isnan(arr_dict[key]), -1, arr_dict[key])
	row, col = np.nonzero(arr_dict['lifet_ch'+str(ch)])
	mincol = min(col)
	minrow = min(row)
	maxcol = max(col)
	maxrow = max(row)
	## compute pairwise eudist
	from scipy.spatial.distance import cdist
	curve_in_use = curvemat['C'][minrow:maxrow+1,mincol:maxcol+1,0:635] ## need plus 1 because maxrow should be included
	nonzerocol = col - mincol
	nonzerorow = row - minrow
	nonzero_coords = np.dstack((nonzerorow,nonzerocol))[0]
	nonzero_coords = np.array(nonzero_coords, 'i')
	nonzero_coords_dict = {tuple(item):i for i, item in enumerate(nonzero_coords)}
	#inp = open("./curveDist_all.pkl","rb")
	#curves_eudist.load(inp)
	#inp.close()
	#curves_eudist = np.ones( [nonzero_coords.shape[0], nonzero_coords.shape[0]] )
	nonzero_curves = curve_in_use[nonzerorow,nonzerocol]
	curves_eudist = cdist(nonzero_curves, nonzero_curves)
	## delete arrays
	del curve_in_use
	del curvemat
	del mat
	####
	neighbor_contrast_dict = {}
	top = 50
	#boolmat = arr_dict['lifet_ch1'][minrow:maxrow+1,mincol:maxcol+1]
	ri, ci = np.meshgrid(range(0, maxrow+1-minrow), range(0, maxcol+1-mincol))
	indexmat = np.dstack( (ci, ri) )
	for i, c in enumerate(nonzero_coords):
		print i, c
		key = str(c[0])+','+str(c[1])
		neighbor_contrast_dict[key] = []
		dist_sum = 0
		count_sum = 0
		#p_row, p_col = np.nonzero(boolmat[c[0]-top:c[0]+top,c[1]-top,c[1]+top)
		#p_row = p_row + np.max(0, c[0]-top)
		#p_col = p_col + np.max(0, c[1]-top)
		#nonzero_neighbor_coords = np.dstack((p_row,p_col))[0]
		for d in range(1, top):
			#print i, ':', c, d
			dist, count = get_bound_sum(c, d)
			dist_sum += dist
			count_sum += count
			if count_sum == 0:
				neighbor_contrast_dict[key].append(-1)
			else:
				neighbor_contrast_dict[key].append(dist_sum/count_sum)
			
	
	print "dump"
	output = open("./neignborDist_ch"+str(ch)+".js","wb")
	print >> output, "var neignborDist_ch" + str(ch) + " = "
	print >> output, neighbor_contrast_dict
	print "end dump"
	output.close()
	


if __name__ == "__main__":
	import sys
	if	len(sys.argv) > 1:
		num = sys.argv[1]
		#mat_to_json("case3_ch"+str(num)+".mat", "./curve_ch"+str(num)+".mat", int(num))
	fat_mat("case3_ch1.mat","case3_fat")
	#print "finish mat_to_json"
	#curve_eudist("case3_ch1.mat", "./curve_ch1.mat", 1)
	#print "finish curve_eudist"
	#curve_eudist_all("case3_ch1.mat", "./curve_ch1.mat")
	#neighbor_contrast("case3_ch1.mat", "./curve_ch1.mat", 1)