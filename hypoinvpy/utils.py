#This file contains utility functions for the package.
import pandas as pd
import os
from importlib import resources as impresources
from hypoinvpy import templates
#####
def get_filelist(dir=None,extension=None,pattern=None,sort=True):
    """
    Get list of files with absolute path, by specifying the format extension.
    Modified from SeisGo.utils.get_filelist()

    ===========PARAMETERS=============
    dir: directory containing the files.
    extension: file extension (the ending format tag), for example "h5" for asdf file.
    pattern: pattern to use in searching. Wildcards are NOT considered here.
    sort: (optional) to sort the list, default is True.

    ============RETURN=============
    flist: the list of file names with paths.
    """
    if dir is None:
        dir="."
    if extension is None:
        flist=[os.path.join(dir,f) for f in os.listdir(dir)]
    else:
        flist=[os.path.join(dir,f) for f in os.listdir(dir) if f[-len(extension):].lower()==extension.lower()]
    if pattern is not None:
        flist2=[]
        for f in flist:
            if f.find(pattern)>=0: flist2.append(f)
        flist=flist2
    if sort:
        return  sorted(flist)
    else:
        return  flist
##    
def stainfo_json2csv(infile,outfile=None):
    """
    Convert seismic station information from json database format to CSV format (with the option of return a Pandas.DataFrame object).
    
    ====== PARAMETERS =======
    infile: input file name of the json database. [required]
    outfile: output file name. Default is None (return the dataframe object). If specified (not None), the dataframe will be saved as a csv file.

    ====RETURN====
    outdata: [optional] return the object when outfile is not specified (None)

    """
    indata=pd.read_json(infile,orient='index')
    indata['station']=indata.index
    indata=indata.reset_index(drop=True)

    # allocate arrays/list for output columns
    net_all=[]
    sta_all=[]
    chan_all=[]
    lat_all=[]
    lon_all=[]
    ele_all=[]
    for i in range(len(indata)):
        # [lat[i], lon[i],ele[i]]=coords_all[i]
        for j in range(len(indata.channels[i])):
            net_all.append(indata.network[i])
            sta_all.append(indata.station[i])
            chan_all.append(indata.channels[i][j])
            lat_all.append(indata.coords[i][0])
            lon_all.append(indata.coords[i][1])
            ele_all.append(indata.coords[i][2])
    outdata=pd.DataFrame(list(zip(net_all,sta_all,chan_all,lat_all,lon_all,ele_all)),columns=['network','station','channel',
                                                                             'latitude','longitude','elevation'])
    if outfile is not None: #save to file.
        fhead=os.path.split(outfile)[0]
        if len(fhead)>0:
            if not os.path.isdir(fhead):os.makedirs(fhead)
        outdata.to_csv(outfile, index=False)
    else: #return the pandas.DataFrame object.
        return outdata

def get_hypinv_template_list():
    hypinvtemplatedir=impresources.files(templates)
    templatelist=get_filelist(hypinvtemplatedir,pattern='hypinv_template_')
    templatetail=[]
    for tf in templatelist:
        ftail=os.path.split(tf)[1]
        templatetail.append(ftail)
    #
    return(templatetail)
#
def load_hypinv_template(template_name=None):
    if template_name is None:
        raise ValueError('template_name NOT specified. run get_hypinv_template_list() to get a list of available templates.')
    else:
        if os.path.isfile(template_name):
            inp_file=template_name
        else:
            inp_file = (impresources.files(templates) / template_name)
        f=open(inp_file)
        lines=f.readlines()
        f.close()
        
        return lines
# format output
def write_csv(fout, line, evid, lat_code, lon_code):
    """
    Write hypoinverse output lines into CSV format with earthquake source parameters.

    ==========PARAMETERS==============
    fout: output CSV file name.
    line: hypoinverse output earthquake parameter line.
    evid: evid for record. 
    lat_code and lon_code: the character labeling the latitude and longitude values. E.g., lat_code='N', lon_code='W' for 
        geographical locations in the northern hemisphere with west longitude.
    """
    #
    grd_ele = 0 #the original code from Hypo_Interface_Py added the correction for grid elevation. Force to 0 now. 
    #Not clear why this was added. Will work on this. Be aware of this point. Noted by Xiaotao for HypoInvPy
    
    mag = 0.0 #force magnitude to 0 for now. Will add magnitude information later.
    codes = line.split()
    date, hrmn, sec = codes[0:3]
    dtime = date + hrmn + sec.zfill(5)
    lat_deg = float(line[20:22])
    lat_min = float(line[23:28])
    lat = lat_deg + lat_min/60 if lat_code=='N' else -lat_deg - lat_min/60
    lon_deg = float(line[29:32])
    lon_min = float(line[33:38])
    lon = lon_deg + lon_min/60 if lon_code=='E' else -lon_deg - lon_min/60
    dep = float(line[38:44])
    
    fout.write('{},{:.4f},{:.4f},{:.1f},{:.2f},{}\n'.format(dtime, lat, lon, dep+grd_ele, mag, evid))

