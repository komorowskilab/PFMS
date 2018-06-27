'''
@attention: Excuting a set of peakfinders in parallel or sequentially
@author: Husen Umer
@organization: LCB centre for BioInformatics at Uppsala University
@since: 7 mar 2011
@version: 1.0

'''

from src import Erange
from src import MACS
from src import FindPeaks
from src import HPeak
from src import SISSR
from src import CISGenome
from src import SeqSite

class ExecutePeakFinders:

    def parallel(self, max_cpu_use, peakfinders, pf_path, pf_options, pf_info, Results):
        from multiprocessing import Pool
        Results = []

        pool = Pool(max_cpu_use)
        if "sissr" in peakfinders:
            pool.apply_async(func=SISSR.SISSR_call, args=(pf_info['label'], pf_info['data_file'], pf_info['pfms_path']+pf_path['sissr'], pf_info['wiggle'], pf_info['control_file'], pf_options['sissr']),callback=Results.append)
        if "macs" in peakfinders:
            if pf_path['macs'] == "macs":
                pool.apply_async(func=MACS.MACS_call, args=(pf_info['label'], pf_info['data_file'], pf_path['macs'], pf_info['wiggle'], pf_info['control_file'], pf_options['macs']), callback=Results.append)
            else:
                pool.apply_async(func=MACS.MACS_call, args=(pf_info['label'], pf_info['data_file'], pf_info['pfms_path']+pf_path['macs'], pf_info['wiggle'], pf_info['control_file'], pf_options['macs']), callback=Results.append)
        if "erange" in peakfinders:
            pool.apply_async(func=Erange.Erange_call, args=(pf_info['label'], pf_info['data_file'], pf_info['pfms_path']+pf_path['erange'], pf_info['wiggle'], pf_info['control_file'], pf_options['erange']), callback=Results.append)
        if "cisgenome" in peakfinders:
            pool.apply_async(func=CISGenome.CISGenome_call, args=(pf_info['label'], pf_info['data_file'], pf_info['pfms_path']+pf_path['cisgenome'], pf_info['wiggle'], pf_info['control_file'], pf_options['cisgenome']), callback=Results.append)
        if "findpeaks" in peakfinders:
            pool.apply_async(func=FindPeaks.FindPeaks_call, args=(pf_info['label'], pf_info['data_file'], pf_info['pfms_path']+pf_path['findpeaks'], pf_info['wiggle'], pf_info['control_file'], pf_options['findpeaks'], pf_info['java_max_memory_usage']), callback=Results.append)
        if "seqsite" in peakfinders:
            pool.apply_async(func=SeqSite.SeqSite_call, args=(pf_info['label'], pf_info['data_file'], pf_info['pfms_path']+pf_path['seqsite'], pf_info['wiggle'], pf_info['control_file'], pf_options['seqsite']), callback=Results.append)
        if "hpeak" in peakfinders:
            pool.apply_async(func=HPeak.HPeak_call, args=(pf_info['label'], pf_info['data_file'], pf_info['pfms_path']+pf_path['hpeak'], pf_info['wiggle'], pf_info['control_file'], pf_options['hpeak']), callback=Results.append)
        
        pool.close()
        pool.join()
        return Results


    def sequential(self, peakfinders, pf_path, pf_options, pf_info, Results):
        Results = []
        
        if "sissr" in peakfinders:
            Results.append(SISSR.SISSR_call(pf_info['label'], pf_info['data_file'], pf_info['pfms_path']+pf_path['sissr'], pf_info['wiggle'], pf_info['control_file'], pf_options['sissr']))
        if "macs" in peakfinders:
            if pf_path['macs'] == "macs":
                Results.append(MACS.MACS_call(pf_info['label'], pf_info['data_file'], pf_path['macs'], pf_info['wiggle'], pf_info['control_file'], pf_options['macs']))
            else:
                Results.append(MACS.MACS_call(pf_info['label'], pf_info['data_file'], pf_info['pfms_path']+pf_path['macs'], pf_info['wiggle'], pf_info['control_file'], pf_options['macs']))
        if "erange" in peakfinders:
            Results.append(Erange.Erange_call(pf_info['label'], pf_info['data_file'], pf_info['pfms_path']+pf_path['erange'], pf_info['wiggle'], pf_info['control_file'], pf_options['erange']))
        if "cisgenome" in peakfinders:
            Results.append(CISGenome.CISGenome_call(pf_info['label'], pf_info['data_file'], pf_info['pfms_path']+pf_path['cisgenome'], pf_info['wiggle'], pf_info['control_file'], pf_options['cisgenome']))
        if "findpeaks" in peakfinders:
            Results.append(FindPeaks.FindPeaks_call(pf_info['label'], pf_info['data_file'], pf_info['pfms_path']+pf_path['findpeaks'], pf_info['wiggle'], pf_info['control_file'], pf_options['findpeaks'], pf_info['java_max_memory_usage']))
        if "seqsite" in peakfinders:
            Results.append(SeqSite.SeqSite_call(pf_info['label'], pf_info['data_file'], pf_info['pfms_path']+pf_path['seqsite'], pf_info['wiggle'], pf_info['control_file'], pf_options['seqsite']))
        if "hpeak" in peakfinders:
            Results.append(HPeak.HPeak_call(pf_info['label'], pf_info['data_file'], pf_info['pfms_path']+pf_path['hpeak'], pf_info['wiggle'], pf_info['control_file'], pf_options['hpeak']))
        

        return Results



