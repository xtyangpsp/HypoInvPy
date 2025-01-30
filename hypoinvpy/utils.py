#This file contains utility functions for the package.
import pandas as pd
import os,json
from importlib import resources as impresources
from hypoinvpy import templates
import warnings
from tqdm import tqdm
from datetime import datetime
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
def stainfo_json2csv(infile,outfile=None,informat=None):
    """
    Convert seismic station information from json database format to CSV format (with the option of return a Pandas.DataFrame object).
    
    ====== PARAMETERS =======
    infile: input file name of the json database. [required]
    outfile: output file name. Default is None (return the dataframe object). If specified (not None), the dataframe will be saved as a csv file.
    informat: input json format, either from EQTransformer (eqt) or GAMMA.
    ====RETURN====
    outdata: [optional] return the object when outfile is not specified (None)

    """
    if informat is None:
        informat='EQT'
        warnings.warn('default input station format (eqt or gamma) not set. will use: '+informat)
    # allocate arrays/list for output columns
    net_all=[]
    sta_all=[]
    chan_all=[]
    lat_all=[]
    lon_all=[]
    ele_all=[]
    if informat.lower() == 'eqt':
        indata=pd.read_json(infile,orient='index')
        indata['station']=indata.index
        indata=indata.reset_index(drop=True)
        for i in range(len(indata)):
            # [lat[i], lon[i],ele[i]]=coords_all[i]
            for j in range(len(indata.channels[i])):
                net_all.append(indata.network[i])
                sta_all.append(indata.station[i])
                chan_all.append(indata.channels[i][j])
                lat_all.append(indata.coords[i][0])
                lon_all.append(indata.coords[i][1])
                ele_all.append(indata.coords[i][2])
    elif informat.lower() == 'gamma':
        with open(infile, "r") as file:
            indata = json.load(file)
        #
        for station,details in indata.items():
            net,sta,loc,ch=station.split('.')
            # print(net,sta,loc,ch)
            component=details['component']
            for i in range(len(component)):
                chan=ch+component[i]
                net_all.append(net)
                sta_all.append(sta)
                chan_all.append(chan)
                lon_all.append(details['longitude'])
                lat_all.append(details['latitude'])
                ele_all.append(details['elevation(m)'])
    ##
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
def write_csv(fout, line, evid, lat_code, lon_code, mag_dict=None):
    """
    Write hypoinverse output lines into CSV format with earthquake source parameters.

    ==========PARAMETERS==============
    fout: output CSV file name.
    line: hypoinverse output earthquake parameter line.
    evid: evid for record. 
    lat_code and lon_code: the character labeling the latitude and longitude values. E.g., lat_code='N', lon_code='W' for 
        geographical locations in the northern hemisphere with west longitude.
    mag_dict: a dictionary containing the magnitude of each event id in the formm of {'id',magnitude}. Default None.
    """
    #
    grd_ele = 0 #the original code from Hypo_Interface_Py added the correction for grid elevation. Force to 0 now. 
    #Not clear why this was added. Will work on this. Be aware of this point. Noted by Xiaotao for HypoInvPy
    
    # mag = 0.0 #force magnitude to 0 for now. Will add magnitude information later.
    codes = line.split()
    date, hrmn, sec = codes[0:3]
    dtime = date + hrmn + sec.zfill(5)
    lat_deg = float(line[20:22])
    lat_min = float(line[23:28])
    lat = lat_deg + lat_min/60 if lat_code=='N' else -lat_deg - lat_min/60
    lon_deg = float(line[29:32])
    lon_min = float(line[33:38])
    lon = lon_deg + lon_min/60 if lon_code=='E' else -lon_deg - lon_min/60
    dep = float(line[38:45])
    mag = float(line[45:52])
    if mag_dict is not None:
        mag = mag_dict[str(evid)]
    
    fout.write('{},{:.4f},{:.4f},{:.1f},{:.2f},{}\n'.format(dtime, lat, lon, dep+grd_ele, mag, evid))

def conv_gamma(eventfile, pickfile, outfile='phase_output.phs',outformat='Y2000',
               qc=False, separator=',',default_component=None,v=False):
    """
    Function to convert picks from GaMMA associator to HypoInvPy format, 
    based on gamma2hypoinverse.py in Weiqiang Zhu's QuakeFlow. (github.com/AI4EPS/QuakeFlow/HypoDD)
    
    Parameters:
    eventfile: name of file which contains the events from GaMMA. 
                The output file of event catalog from the GaMMA implementation from QuakeFlow above.
    pickfile: name of file which contains the picks from GaMMA. 
                The output file for phase picks name from the GaMMA implementation from QuakeFlow above.
    outfile: name of the output phase file. Default is 'phase_output.phs'.
    outformat: output phase data format, based on formats for Hypoinverse. Default is 'Y2000'.
    qc: If true, this converter will check to ensure that every event has both 
                P picks and S picks (events with only one or the other will cause HypoInverse to 
                fail). Setting to False omits this check. Defaults to False.
    separator: line separator, defaults to ','.
    default_component: force to use the specified component if the comp_code in the pick file is empty.
                Only one letter is allowed. Default is None.
    """
    if default_component is not None and len(default_component)>1:
        raise ValueError('default_component can only one one letter. Current %s has a length of %d.'
                         %(default_component,len(default_component)))
    picks=pd.read_csv(pickfile, sep=separator)
    events=pd.read_csv(eventfile, sep=separator)

    # Remove duplicate columns
    picks = picks.loc[:, ~picks.columns.duplicated()]
    # Rename 'event_idx' to 'event_index' to be consistent with the catalog
    picks.rename(columns={"event_idx": "event_index"}, inplace=True)

    picks_eventwise = picks.groupby("event_index").groups

    outfile_handle=open(outfile,'w')
    # print(events)
    for i in tqdm(range(len(events))):
        temporarybuffer = [] #don't add lines directly to file, catch here and check first
        has_p = False
        has_s = False

        event = events.iloc[i]
        # print(event)
        event_time = datetime.strptime(event["time"], "%Y-%m-%dT%H:%M:%S.%f").strftime("%Y%m%d%H%M%S%f")[:-4]

        lat_degree = int(event["latitude"])
        lat_minute = (event["latitude"] - lat_degree) * 60 * 100
        south = "S" if lat_degree <=0 else " "

        lon_degree = int(event["longitude"])
        lon_minute = (event["longitude"] - lon_degree) * 60 * 100
        east = "E" if lat_degree >=0 else " "
        event_mag = event["magnitude"]
        depth = event["depth(m)"] / 1000
        event_line = f"{event_time}{abs(lat_degree):2d}{south}{abs(lat_minute):4.0f}{abs(lon_degree):3d}{east}{abs(lon_minute):4.0f}{depth:5.0f}" 
        event_final = event_line + f"{event_mag:3.1f}{' ':97}{event['event_index']:10}\n"
        temporarybuffer.append(event_final)
        if v:
            print(event_line)

        picks_idx = picks_eventwise[event["event_index"]]
        for j in picks_idx:
            pick = picks.iloc[j]

            network_code, station_code, comp_code, channel_code = pick['id'].split('.')
            phase = pick['type']
            phase_weight = min(max(int((1-pick['prob']) / (1 - 0.3) * 4) - 1, 0), 3)
            picktime = datetime.strptime(pick["timestamp"], "%Y-%m-%dT%H:%M:%S.%f")
            pickmin = picktime.strftime("%Y%m%d%H%M")
            picksec = picktime.strftime("%S%f")[:-4]
            phase_amplitude=pick['phase_amplitude']

            if default_component is not None and len(comp_code)<1:
                comp_code = default_component

            templine = f"{station_code:<5}{network_code:<2}  {channel_code:<2}{comp_code:<1}"
            
            if phase.upper() == 'P':
                has_p = True
                pick_line = f"{templine:<13} P {phase_weight:<1d}{pickmin} {picksec}                    {phase_amplitude:<7}"
            elif phase.upper() == "S":
                has_s = True
                pick_line = f"{templine:<13}   4{pickmin} {'':<12}{picksec} S {phase_weight:<1d}    {phase_amplitude:7f}"
            else:
                warnings.warn("Phase Type Error: "+phase.upper())
                continue
            if v:
                print(pick_line)
            temporarybuffer.append(pick_line + "\n")
        
        if len(temporarybuffer)>1:
            #add lines to file if there is both a P and S pick, or if filter is off.
            has_pick=False
            count=0
            for line in temporarybuffer:
                if qc:
                    if has_p and has_s:
                        outfile_handle.write(line)
                        count = count + 1
                        has_pick=True
                else:
                    outfile_handle.write(line)
                    count = count + 1
                    has_pick=True
            #

            if has_pick:
                outfile_handle.write("\n")
                if v:
                    print('-> Saved %d phase picks. '%(count))
            else:
                if v:
                    print('-> No pick to save. Picks might be dropped after QC.')
        else:
            continue
    
    outfile_handle.close()
