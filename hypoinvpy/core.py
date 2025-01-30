#This module contains core functions for running HypoInvPy interface.
#Import needed packages first.
import pandas as pd
import os
import numpy as np
from hypoinvpy import utils
import subprocess
#
#
def reformat_stainfo(infile,outfile, informat='csv',channel=True,channel_default='HHZ',
                     force_channel_type=None,lat_code='N', lon_code='W',rename_component_dict=None,
                     ignore_component=False):
    """
    Reformat station information file for hypoinver run. Input CSV needs to have at least the following columns:
    network,station,channel,latitude,longitude,elevation.

    Additional columns will be ignored.

    ====== PARAMETERS =======
    infile: input file name (including path).
    outfile: output file name (including path).
    informat: input file format, choose from 'csv', 'json', 'json-eqt','json-gamma'. Default is 'csv'.
    channel: is channel column available. Default is True. Otherwise, channel_default is used.
    channel_default: when channel is not available (channel=False), this value will be used. Default='HHZ'
    force_channel_type: treat all channel types (i.e., EH?, BH?, or HH?) as the same. Default=None (use true channel).
    rename_component_dict: dictionary used to rename component label, e.g., HH1 to HHN. Default is {'1':'N','2':'E'}
    lat_code: code for latitudes. Default is 'N'.
    lon_code: code for longitudes. Default is 'W'. lat_code and lon_code are used to format the coordinates.
    ignore_component: if True, only save the channel type information, such as BH instead of BHZ. Default False.
    """
    if force_channel_type is not None:
        if len(force_channel_type) != 2:
            raise ValueError('force_channel_type has to be two characters. Wrong value: '+force_channel_type)
    if informat.lower() == 'csv':
        if infile[-4:].lower() =='json':
            raise ValueError(infile+' may be a json file. change informat argument to json.')
        indata=pd.read_csv(infile)
    elif informat.lower() == 'json' or informat.lower() == 'json-eqt':
        if infile[-3:].lower() =='csv':
            raise ValueError(infile+' may be a csv file. change informat argument to csv.')
        indata=utils.stainfo_json2csv(infile)
    elif informat.lower() == 'json-gamma':
        if infile[-3:].lower() =='csv':
            raise ValueError(infile+' may be a csv file. change informat argument to csv.')
        indata=utils.stainfo_json2csv(infile,informat='gamma')
    else:
        raise ValueError('input file format of %s not recoganized. Use "csv" or "json".'%(informat))
    #
    if rename_component_dict is not None: 
        rename_keys=list(rename_component_dict.keys())

    fhead=os.path.split(outfile)[0]
    if len(fhead)>0:
        if not os.path.isdir(fhead):os.makedirs(fhead)
    fout = open(outfile,'w')
    for i in range(len(indata)):
        net,sta,lat,lon,ele=[indata['network'][i],indata['station'][i],indata['latitude'][i],indata['longitude'][i],indata['elevation'][i]]
        lat, lon, ele = abs(lat), abs(lon), int(ele)
        lat_deg = int(lat)
        lat_min = 60*(lat-int(lat))
        lon_deg = int(lon)
        lon_min = 60*(lon-int(lon))
        lat = '{:2} {:7.4f}{}'.format(lat_deg, lat_min, lat_code)
        lon = '{:3} {:7.4f}{}'.format(lon_deg, lon_min, lon_code)
        if not channel:
            chan=channel_default
        else:
            chan=indata['channel'][i]
        if force_channel_type is not None:
            chan=force_channel_type+chan[-1]
        if rename_component_dict is not None: #force to rename component label.
            if chan[-1] in rename_keys:
                chan=chan[0:2]+str(rename_component_dict[chan[-1]])
        if ignore_component:
            chan=chan[0:2] #drop the component information.
        # hypoinverse format 2
        fout.write("{:<5} {}  {}  {}{}{:4}\n".format(sta, net, chan,lat, lon, ele))
    fout.close()

#
#
def generate_parfile(config,pardir='input',outdir='output',template=None,magline=None):
    """
    Generate parameter files for running hypoinverse, with given HypoInvConfig object.

    ========PARAMETERS=========
    config: A HypoInvConfig object containing all controling parameters.
    pardir: output directory to store the parameter files. Default: input.
    outdir: hypoinverse relocation output directory. Default: output.
    template: provide template or use the built-in template file for parameters.
    magline: the line in the parameter file to compute magnitude. Defalt is None.

    ========RETURNS===========
    filelist: list of the parameter files, including path.
    """
    if not os.path.isdir(pardir):
        os.makedirs(pardir)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    filelist=[]
    for ztr in config.ztrlist:
        # set control file
        fhyp = os.path.join(pardir,'%s-%s.hyp'%(config.run_tag, ztr))
        filelist.append(fhyp)

        #save parameters by modifying the template parameters.
        fout=open(fhyp,'w')
        if template is None:
            lines=utils.load_hypinv_template(config.template_parfile)
        else:
            lines=utils.load_hypinv_template(template)
        for line in lines:
            # loc params
            if line[0:3]=='ZTR': line = "ZTR %s F \n"%ztr
            if line[0:3]=='RMS': line = "RMS %s \n"%config.rms_weight
            if line[0:3]=='DI1': line = "DI1 %s \n"%config.dist_initial
            if line[0:3]=='DIS': line = "DIS %s \n"%config.dist_weight
            if line[0:3]=='WET': line = "WET %s \n"%config.weight_code
            # i/o paths
            if line[0:3]=='STA': line = "STA '%s' \n"%config.station_file
            if line[0:3]=='PHS': line = "PHS '%s' \n"%config.phase_file
            if line[0:5]=='CRE 1': line = "CRE 1 '%s' %s T \n"%(config.pmodel, config.ref_ele)
            if line[0:5]=='CRE 2': line = "CRE 2 '%s' %s T \n"%(config.smodel, config.ref_ele)
            if line[0:3]=='POS': line = "POS %s \n"%(config.pos)
            if line[0:3]=='SUM': line = "SUM '%s/%s-%s.sum' \n"%(outdir,config.run_tag, ztr)
            if line[0:3]=='MIN': line = "MIN %d \n"%(config.min_nsta)
            if line[0:3]=='PRT': 
                line = "PRT '%s/%s-%s.ptr' \n"%(outdir,config.run_tag, ztr) if config.get_prt else ''
            if line[0:3]=='ARC': 
                line = "ARC '%s/%s-%s.arc' \n"%(outdir,config.run_tag, ztr) if config.get_arc else ''
            #if line[0:3]=='H71': line = "H71 1 1 3" #use hypoinverse summary output format (first integer)
            if line[0:3]=='STO': 
                continue
            fout.write(line)
        #get magnitude
        if magline is not None:
            # line="MAG 1 T 1 1\n"
            line = "MFL '%s/%s-%s.mag' \n"%(outdir,config.run_tag, ztr)
            fout.write(line)
            fout.write(magline)
        line="STO" #stop the program.
        fout.write(line)

        fout.close()
    #
    return filelist
####
def merge_summary(filelist,file_good,file_bad,lat_code,lon_code,mag_dict=None):
    """
    Extract earthquake parameters based on final quality after merging all summary files.

    ======PARAMETERS======
    filelist: summary file list from hypoinverse.
    file_good,file_bad: file names to store good and bad earthquakes.
    lat_code,lon_code: latitude and longitude codes.
    mag_dict: magnitude dictionary in the form of {'id',mag} for all events. Default None (use hypoinverse output).
            mag_dict could also be specified as a catalog file in csv format. E.g., the catalog from GAMMA.
    """
    if mag_dict is not None:
        if isinstance(mag_dict,str):
            events=pd.read_csv(mag_dict)
            mag_dict_use=dict()
            for i in range(len(events.time)):
                event = events.iloc[i]
                # print(event['event_index'].astype(str))
                mag_dict_use[event['event_index'].astype(str)] = event['magnitude']
            #
        elif isinstance(mag_dict,dict):
            mag_dict_use = mag_dict
        else:
            raise ValueError('mag_dict is wrong in type. CSV catalog or dictionary.')
    else:
        mag_dict_use = mag_dict
    fout_bad = open(file_bad,'w')
    fout_good = open(file_good,'w')
    
    # read sum files
    sum_dict = {}
    for fsum in filelist:
        f=open(fsum); sum_lines=f.readlines(); f.close()
        for sum_line in sum_lines:
            evid = sum_line.split()[-1]
            if evid not in sum_dict: sum_dict[evid] = [sum_line]
            else: sum_dict[evid].append(sum_line)
    
    # merge sum lines
    for evid, sum_lines in sum_dict.items():
        sum_list = []
        dtype = [('line','O'),('is_loc','O'),('qua','O'),('azm','O'),('npha','O'),('rms','O')]
        for sum_line in sum_lines:
            codes = sum_line.split()
            is_loc = 1 # whether loc reliable
            if '-' in codes or '#' in codes: is_loc = 0
            qua = sum_line[80:81]
            npha = 1 / float(sum_line[52:55])
            azm  = float(sum_line[56:59])
            rms  = float(sum_line[64:69])
            sum_list.append((sum_line, is_loc, qua, azm, npha, rms))
        sum_list = np.array(sum_list, dtype=dtype)
        sum_list = np.sort(sum_list, order=['qua','azm','npha','rms'])
        sum_list_loc = sum_list[sum_list['is_loc']==1]
        num_loc = len(sum_list_loc)
        # if no reliable loc
        if num_loc==0:
            sum_list_loc = sum_list
            utils.write_csv(fout_bad, sum_list_loc[0]['line'], evid,lat_code,lon_code,mag_dict=mag_dict_use)
        else:
            utils.write_csv(fout_good, sum_list_loc[0]['line'], evid,lat_code,lon_code,mag_dict=mag_dict_use)
    
    fout_bad.close()
    fout_good.close()

    print('Earthquakes are saved in: '+file_good+' and '+file_bad+' for good and bad sources.')
#####
#####
class HypoInvConfig(object):
    """
    Container class to store key configuration parameters for running HypoInverse.
    """
    def __init__(self,phase_file=None,station_file=None,pmodel=None,smodel=None,poisson=1.73,
               run_tag='hyp',hypoinv_bin='hyp1.40',get_prt=False,get_arc=False,
               lat_code='N',lon_code='W',ref_ele=0.0,grd_ele=0.0,
               ztrlist = np.arange(0,20,1),rms_weight='4 0.3 1 3',dist_initial = '1 50 1 2',
               dist_weight = '4 20 1 3',weight_code='1 0.6 0.3 0.2',min_nsta=4):
        self.run_tag = run_tag
        self.hypoinv_bin=hypoinv_bin
        # i/o paths
        self.station_file = station_file
        self.phase_file = phase_file
        self.get_prt = get_prt
        self.get_arc = get_arc
        # geo ref
        self.lat_code = lat_code
        self.lon_code = lon_code
        self.ref_ele = ref_ele # ref ele for CRE mod (max sta ele)
        self.grd_ele = grd_ele # typical station elevation
        # loc params
        self.ztrlist = ztrlist #initial depth for the relocation run.
        self.p_weight = 0 # weight code index
        self.s_weight = 1
        self.min_nsta=min_nsta
        self.rms_weight = rms_weight
        self.dist_initial = dist_initial
        self.dist_weight = dist_weight
        self.weight_code = weight_code
        self.template_parfile = utils.get_hypinv_template_list()[1]
        self.pmodel = pmodel #'input/velo_p_eg.cre'
        self.smodel = smodel #[None, 'input/velo_s_eg.cre'][1]
        self.poisson = poisson #1.73 # provide smod or pos
      
def run_hypoinv(parfilelist,hyp_bin='hyp1.40'):
    for fhyp in parfilelist:
        # 2. run hypoinverse
        p = subprocess.Popen([hyp_bin], stdin=subprocess.PIPE,encoding='utf-8')
        s = "@{}".format(fhyp) + '\n'
        p.communicate(s)
