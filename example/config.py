""" Configure file for HypoInverse interface
    Download hypoInverse at https://www.usgs.gov/software/hypoinverse-earthquake-location
"""
import numpy as np

class Config(object):
  def __init__(self):

    self.ctlg_code = 'eg_hyp'
    # i/o paths
    self.fsta = 'input/station_eg.csv'
    self.fsta_out = 'output/%s.sta'%self.ctlg_code
    self.fpha = 'input/eg.pha'
    self.fpha_out = 'output/%s_hypo.phs'%self.ctlg_code
    self.out_ctlg = 'output/%s.ctlg'%self.ctlg_code
    self.out_pha = 'output/%s.pha'%self.ctlg_code
    self.out_pha_full = 'output/%s_full.pha'%self.ctlg_code
    self.out_sum = 'output/%s.sum'%self.ctlg_code
    self.out_bad = 'output/%s_bad.csv'%self.ctlg_code
    self.out_good = 'output/%s_good.csv'%self.ctlg_code
    self.get_prt = False
    self.get_arc = False
    self.keep_fsums = False
    # geo ref
    self.lat_code = 'N'
    self.lon_code = 'W'
    self.ref_ele = 2.5 # ref ele for CRE mod (max sta ele)
    self.grd_ele = 1.5 # typical station elevation
    self.fsums = 'output/%s-*.sum'%self.ctlg_code
    # loc params
    self.num_workers = 5
    self.ztr_rng = np.arange(0,20,1)
    self.p_wht = 0 # weight code index
    self.s_wht = 1
    self.rms_wht = '4 0.3 1 3'
    self.dist_init = '1 50 1 2'
    self.dist_wht = '4 20 1 3'
    self.wht_code = '1 0.6 0.3 0.2'
    self.fhyp_temp = 'temp_hyp/temp_vp-vs.hyp'
    self.pmod = 'input/velo_p_eg.cre'
    self.smod = [None, 'input/velo_s_eg.cre'][1]
    self.pos = 1.73 # provide smod or pos
