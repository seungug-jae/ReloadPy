#!/usr/bin/env python3
import sys
import os
import itertools
from typing import OrderedDict
import numpy as np
from text_io import skip_lines, skip_to, read_next_line, convert_str2num

_reload_options = ['null_refueling', 'num_reload_comp', 'objective_func', 'reload_refpwr', 'update_refpwr', \
    'dif3d_check', 'dif3d_check_all', 'use_persent', 'num_persent_case', \
    'normalize_power', 'fixed_peak_power', 'peak_power_limit', 'assembly_power_limit', \
    'consider_burnup', 'min_burnup', 'burnup_data', \
    'consider_keff', 'min_deltak', 'max_deltak', 'consider_worth', 'max_worth', \
    'adjust_cycle_length', 'min_cycle_length', 'target_eoc_keff', 'deprate_ratio', 'eps_cycle', \
    'num_fuel_region', 'fuel_regions', 'candidate_range_size', 'candidate_range', 'fixed_sequence', 'margin_power_limit', \
    'eps_power', 'wgt_pekfac']
_numeric_options = ['num_cycles', 'num_past_cycles','num_reload_comp','num_persent_case','peak_power_limit',\
    'assembly_power_limit','min_burnup','min_deltak','max_deltak','max_worth', 'min_cycle_length', 'target_eoc_keff', \
    'disperse_mindist', 'disperse_numpast', 'initial_sequence', 'candidate_range', 'candidate_range_size', \
    'eps_cycle', 'eps_keoc', 'deprate_ratio','ufrac_mean','ufrac_std','salt_dens', 'salt_frac','niso_salt', 'initial_cycle', 'margin_power_limit',\
    'eps_power','wgt_pekfac']
_defined_objectives = ['assembly_power','peak_power','burnup','reactivity','power_shape','power_deep','reactivity_slow']

def reload_input_scan(file):
    ''' read input file '''
    with open(file,'r') as f:      # read fresh fuel composition from reload.in
        # read general control options, may define a class to handle this but unnecessary for the simple case
        # default options -----------------------------------------------------------------------------------------------------------
        options = {'case_name':'test', 'rebus_exe':'rebus.x', 'dif3d_exe':'dif3d.x', 'persent_exe':'persent.x', \
            'reload_utility':'reload.x', 'type14_utility':'mk1415.x', 'normalize_power':'T'}
        options['objective_func']      = 'power_shape'
        options['reload_refpwr']       = 'bol'  # 'bol', 'flat', 'user', relevant only if objective_func = power_shape / power_deep
        options['margin_power_limit']  = 0.0
        options['eps_power']           = 1.0e-3
            # default reference power distribution is from BOL, 'flat' means all fuel assemblies have the core-averaged power
            # 'user' means a 'reload.refpwr' file needs to be provided by users
        options['update_refpwr']       = 'F'   # whether the reference power distribution is updated, needed only if objective_func = power_shape / power_deep
        options['num_reload_comp']     = 1
        options['fixed_peak_power']    = 'T'
        options['peak_power_limit']    = -1.0
        options['assembly_power_limit']= -1.0
        options['consider_burnup']     = 'T'
        options['min_burnup']          = -1.0
        options['burnup_data']         = 'reload.burnup'
        options['consider_keff']       = 'T'
        options['min_deltak']          = -1.0
        options['max_deltak']          = -1.0
        options['consider_worth']      = 'T'
        options['max_worth']           = 0.002
        options['use_persent']         = 'F'
        options['num_persent_case']    = 1
        options['dif3d_check']         = 'F'
        options['dif3d_check_all']     = 'F'
        options['adjust_cycle_length'] = 'T'
        options['min_cycle_length']    = 1.0   # one day
        options['target_eoc_keff']     = 1.002
        options['deprate_ratio']       = None  # depletion rate correction factor
        options['eps_cycle']           = 0.1   # 10%, relative tolerance for cycle length comparison
        options['eps_keoc']            = -1.0  # absolute tolerance for keff comparison during cycle length iteration, user can specify a positive value as fixed tolerance
            # default is -1.0, which means the tolerance will be determined as: [eps_cycle * min_cycle_length * (boc_keff-eoc_keff) / cycle_length]
        options['disperse_reload']     = 'F'
        options['disperse_mindist']    = -1.0   # meant to exclude previous positions directly
        options['disperse_numpast']    = 0      # (maximum) number of previously reloaded assemblies for which the distance is checked, default value 0 means all previous positions
        options['initial_sequence']    = None
        options['initial_cycle']    = None
        options['null_refueling']      = 'F'
        options['candidate_range']     = None
        options['candidate_range_size']= 0
        options['num_past_cycles']     = 0
        options['rebus_input']         = None
        options['random_enri']         = False
        options['fixed_sequence']      = False
        options['wgt_pekfac']          = 1.0
        # end default options -----------------------------------------------------------------------------------------------------------
        
        # flags not specified through input
        options['edit_region_fluence'] = True  # whether REBUS edited both peak discharge burnup and fast fluence
        
        # start to read input
        skip_to(f,'&control')
        flag = None
        while 1:
            line = read_next_line(f,'!')
            if line:
                if line[0] == '/': break
                key_in_line = '=' in line
                line = line.split('=')
                if len(line) == 1:
                    if key_in_line: # first line of certain parameter with values on continued lines
                        flag = line[0].strip().lower()
                        continue
                    elif flag:    # continued line
                        data = line[0].split()
                    else:
                        raise ValueError(f'no option name specified for input line {line[0]}')
                else:
                    flag = line[0].strip().lower()
                    data = line[1].split()
                if flag in _numeric_options:
                    data = convert_str2num(data)
                if len(line) == 1:    # continued line of values
                    if not isinstance(options[flag], list): options[flag] = [options[flag]]
                    options[flag].extend(data)
                else:
                    options[flag] = data[0] if len(data) == 1 else data
            else:
                raise ValueError('input section &control is not ended with "/"')
        if 'num_cycles' not in options:
            raise ValueError(f'total number of depletion cycles (num_cycles) must be specified in {file}')
        options['reload_refpwr'] = options['reload_refpwr'].lower()
        
        # read reloaded compositions
        reload_composition = {}
        skip_to(f,'&material')
        materials = {}
        name = None
        while 1:
            line = read_next_line(f,'!')
            if not line: raise ValueError('input section &material is not ended with "/"')
            if line[0] == '/': break    # end of &material block
            if 'composition_name' in line:
                line = line.split('=')
                index = int(line[0].split('(')[1].split(')')[0])
                name = line[1].strip()
                materials[index] = name
                reload_composition[name] = []      # [[isotope_name, atom_density], ...]
            elif 'composition' in line:     # this line should not have atom density input
                index = int(line.split(',')[1].split(')')[0])
                if index in materials:
                    name = materials[index]
                else:
                    raise ValueError('input error: composition name should be given before atom densities')
            else:
                if name is None:
                    raise ValueError('input error: atom densities given before composition name')
                data = line.split()
                ndata = len(data)
                if ndata%2 != 0: raise ValueError('atomic composition should be given as pairs of isotope and atom density')
                for i in range(0,len(data),2):
                    reload_composition[name].append([data[i], float(data[i+1])])
        options['num_reload_comp'] = len(reload_composition)

        # read random composition data
        if (options['random_enri']=='T'):
            skip_to(f, '&random_comp')
            flag = None
            while 1:
                line = read_next_line(f,'!')
                if line:
                    if line[0] == '/': break
                    key_in_line = '=' in line
                    line = line.split('=')
                    if len(line) == 1:
                        if key_in_line: # first line of certain parameter with values on continued lines
                            flag = line[0].strip().lower()
                            continue
                        elif flag:    # continued line
                            data = line[0].split()
                        else:
                            raise ValueError(f'no option name specified for input line {line[0]}')
                    else:
                        flag = line[0].strip().lower()
                        data = line[1].split()
                    if flag in _numeric_options:
                        data = convert_str2num(data)
                    if len(line) == 1:    # continued line of values
                        if not isinstance(options[flag], list): options[flag] = [options[flag]]
                        options[flag].extend(convert_str2num(data))
                    else:
                        options[flag] = data[0] if len(data) == 1 else data
                else:
                    raise ValueError('input section &control is not ended with "/"')
        
    return (options, reload_composition)

# ----------------------------------------------------------------------
def find_duplicate_numbers(data, index=False):
    from collections import Counter, defaultdict
    if index:
        D = defaultdict(list)
        for i, item in enumerate(data):
            D[item].append(i)
        return {k:v for k,v in D.items() if len(v)>1}
    else:
        return [k for k,v in Counter(data).items() if v>1]    

# ----------------------------------------------------------------------
def reload_input_check(input_file, options, output):
    output.write(' *** control options from input file: '+input_file+' (options might be modified by code)\n')
    output.write(f'     working directory: run_{options["case_name"]}\n\n')
    # check and modify input options 
    if options['null_refueling'] == 'T':
        options['use_persent'] = 'T'
        options['adjust_cycle_length'] = 'F'
        options['deprate_ratio'] = None
        options['disperse_reload'] = 'F'
        options['initial_sequence'] = None
    else:
        if options['consider_worth'] == 'T':
            options['use_persent'] = 'T'
        if options['adjust_cycle_length'] == 'T':
            options['consider_keff'] = 'T'
            options['min_deltak'] = 0.000001
            if options['min_burnup'] < 0.0:
                options['min_burnup'] = 1.0e-8 # relax this constraint
        if options['objective_func'] not in _defined_objectives:
            raise ValueError(f'invalid input for objective_func, valid options are: {_defined_objectives}')
        if options['objective_func'] == 'burnup':  
            options['consider_burnup'] = 'T'   # burnup evaluation is required
        elif options['objective_func'] == 'reactivity':
            options['use_persent'] = 'T'       # reactivity evaluation is required
        if options['objective_func'] in ['power_shape','power_deep']:
            if options['reload_refpwr'] not in ['bol','flat','user']:
                raise ValueError('invalid input for reload_refpwr, valid options are: bol, flat, or user')
            if options['reload_refpwr'] == 'flat': options['update_refpwr'] = 'F'
        else:
            options['reload_refpwr'] = None
            options['update_refpwr'] = None
        if options['consider_keff'] == 'T' and options['use_persent'] == 'F':
            options['dif3d_check'] = 'T'
        if options['initial_sequence'] and len(options['initial_sequence']) % 2 != 0:
            raise ValueError('manual reloading sequence for a few initial cycles should be specified as pairs of ring and position numbers')
        if options['num_persent_case'] < 0: 
            raise ValueError('invalid input for num_persent_case, which cannot be negative')

    if options['candidate_range']:
        if len(options['candidate_range']) % 2 != 0:
            raise ValueError('prescribed candidate positions should be specified as pairs of ring and position numbers')
        else:
            options['candidate_range_size'] = len(options['candidate_range']) // 2
            RPs = np.reshape(options['candidate_range'],(options['candidate_range_size'],2))
            rr = find_duplicate_numbers(RPs[:,0],index=True)
            dup = []
            for r in rr:
                pp = find_duplicate_numbers([RPs[i,1] for i in rr[r]])
                dup.extend([f'({r},{p})' for p in pp])
            if len(dup) > 0:
                raise ValueError('duplicated candidate positions in the user input:\n'+' '.join(dup))

    if options['consider_worth'] == 'T':
        options['num_reload_comp'] -= 1
        if options['num_reload_comp'] < 1:
            output.write('Input Error: at least one fresh fuel and one empty assembly composition should be provided in &material section')
            raise ValueError('input error: no empty assembly composition provided in &material section')

    # check executable paths
    for exe in [options['rebus_exe'], options['reload_utility'], options['type14_utility']]:
        if not os.path.exists(exe):
            raise ValueError(f'required executable {exe} is not found')
    dif3d_exe, persent_exe = options['dif3d_exe'], options['persent_exe']
    if options['consider_keff'] == 'T':
        if not os.path.exists(dif3d_exe):
            raise ValueError(f'DIF3D executable {dif3d_exe} is not found')
        if options['use_persent'] == 'T':
            if not os.path.exists(persent_exe):
                raise ValueError(
                    f'PERSENT executable "{persent_exe}" is not found')
    elif options['dif3d_check'] == 'T':
        if not os.path.exists(dif3d_exe):
            raise ValueError(f'DIF3D executable {dif3d_exe} is not found')

    # print a copy of control options
    for opt in options:
        if options[opt] is None:
            continue
        elif isinstance(options[opt], list):
            if len(options[opt]) == 0: continue
            output.write(f' {opt:25s} = ' + ' '.join([f'{val}' for val in options[opt]]) + '\n')
        else:
            output.write(f' {opt:25s} = {options[opt]}\n')

    # check rebus input and make sure flux (RTFLUX/NHFLX0) and power density (PWDINT) files are saved
    rebus_input = options['rebus_input']
    if rebus_input is None and options['num_past_cycles'] > 0:
        rebus_input = f'run_{options["case_name"]}/rebus.in_{options["num_past_cycles"]+1}'
    elif rebus_input is None:
        raise ValueError('Input Error: REBUS input file not specified')
    with open(rebus_input,'r') as fin:
        line = skip_to(fin,'=A.STP027')
        freeform = 'UNFORM' in line
        while True:
            line = fin.readline()
            if line == '': break
            if 'UNFORM=' in line or 'DATASET=' in line: break
            if line[:3] == '03':
                if freeform:
                    flag = line[2:].split()
                    if len(flag) < 7 or int(flag[6]) == 0: options['edit_region_fluence'] = False
                    if len(flag) < 8 or int(flag[7]) == 0: options['edit_region_fluence'] = False
                else:
                    if int(line[42:48]) == 0 or int(line[48:54]) == 0: options['edit_region_fluence'] = False
            if line[:2] == '06':
                if freeform:
                    flag = int(line[2:].split()[0])
                else:
                    flag = int(line[11:18].strip())
                if flag != 0:
                    raise ValueError('Input Error: flux/power data files should be saved at each time step,'\
                        ' check type 06 card of REBUS input A.STP027')

    # check VARIANT option
    with open(rebus_input,'r') as fin:
        line = skip_to(fin,'=A.DIF3D')
        freeform = 'UNFORM' in line
        if not freeform and options['use_persent']=='T':
            if output: output.write('Input Error: when PERSENT is invoked, A.DIF3D must be provided in free format\n')
            raise ValueError('when PERSENT is invoked, A.DIF3D must be provided in free format')
        enabled_VARIANT = False
        while 1:
            line = fin.readline()
            if line[:2] == '  ': continue
            if 'UNFORM=' in line or 'DATASET=' in line: break
            if line[:2] == '12':
                enabled_VARIANT = True
                break
        options['use_VARIANT'] = enabled_VARIANT
        if enabled_VARIANT:
            order = abs(int(line.split()[2]))
            if(order > 9999): order = order%10000
            options['VARIANT_Pn_order'] = order//100 if order > 99 else order//10
        else:
            options['VARIANT_Pn_order'] = None
            output.write('Warning: if VARIANT option is not enabled in A.DIF3D, PERSENT will use default angular approximation (usually P3)\n')

# ----------------------------------------------------------------------
def prepare_utility_input(code,options):
    from rebus_files import fetch_anip3_info
    info = fetch_anip3_info('ANIP3')
    if code == 'reload':
        with open('reload.in','w') as f:
            line = '&control\n'            
            indent = 23
            for opt in _reload_options:
                if options[opt] is None: continue
                if isinstance(options[opt], list):
                    if len(options[opt]) == 0: continue
                    val = options[opt][0]
                    form = '{:13.5E} ' if isinstance(val,float) else '{:} '
                    val = form.format(val)
                    j = 1
                    for i in range(1, len(options[opt])):
                        val += form.format(options[opt][i])
                        j += 1
                        if j == 10:
                            val += '\n' + ' '*indent
                            j = 0
                    if j != 0: val += '\n'
                else:
                    val = options[opt]
                line += f' {opt:19s} = {val}\n'
            if info['nsect'] > 0:
                line += f' num_hex_sector = {info["nsect"]}\n'
            f.write(line+'/\n')

            f.write('\n&material\n')
            with open('reload.comp','r') as f0:
                while 1:
                    line = read_next_line(f0, comment_tag='!')
                    if line.strip() == '&material': break 
                while 1:
                    line = read_next_line(f0, comment_tag='!')
                    if line.strip() == '/': break 
                    f.write(line+'\n')
            f.write('/\n')
    else:
        raise ValueError(f'unrecorganized code {code}')

# ----------------------------------------------------------------------
def read_estimated_powers(file):
    # read estimated peak assembly powers with reloaded assemblies
    line = skip_to(file,'*** Estimating assembly powers (Watt) with reloaded composition')
    skip_lines(file,4)
    reload_pekass = {}      # {reload_pos:[peak_ass,peak_assembly_power]}
    reload_asspwr = {}      # {reload_pos: [eoc power of reloaded assembly, reloaded assembly power]}
    while 1:
        line = file.readline().split()
        if len(line) == 0: break
        reload_pekass[int(line[0])] = [int(line[-1]), float(line[-2])]
        reload_asspwr[int(line[0])] = [float(line[4]), float(line[5])]
    # read estimated peak power densities with reloaded assemblies
    line = skip_to(file,'*** Estimating peak power density (W/cc) with reloaded composition')
    skip_lines(file,4)
    reload_peak = {}        # {reload_pos:[peak_ass,peak_node,peak_density]}
    while 1:
        line = file.readline().split()
        if len(line) == 0: break
        reload_peak[int(line[0])] = [int(line[6]), int(line[7]), float(line[5])]
    return (reload_pekass, reload_peak, reload_asspwr)

# ----------------------------------------------------------------------
def read_estimated_worths(file):
    line = skip_to(file,'*** Estimating reactivity change due to reloaded composition')
    skip_lines(file,5)
    data = {}
    while 1:
        line = file.readline().split()
        if len(line) == 0: break
        data[int(line[0])] = [float(line[4]), float(line[5])]
    return data

# ----------------------------------------------------------------------
def reload_output_scan(file):
    ''' read reload.out and aggregate power, reactivity, burnup, etc. for each cycle
    returned data:
        dict :: results
            'null_refueling'        : True/False
            'peak_power_limit'      : peak power limit
            'assembly_power_limit'  : assembly-integerated power limit
            'boc_peak'              : [IASS, IRING, IPOS, INODE, peak_power_density], peak power density at BOC
            'eoc_peak'              : [IASS, IRING, IPOS, INODE, peak_power_density], peak power density at EOC
            'boc_pekass'            : [IASS, IRING, IPOS, peak_assembly_power], peak assembly power at BOC
            'eoc_pekass'            : [IASS, IRING, IPOS, peak_assembly_power], peak assembly power at EOC
            'this_cycle_length'     : cycle length for this cycle
            'reload_comp'           : name of chosen reload composition (useful when multiple reload compositions are tested)
                        
            'reload_results'        : {keyword:result, ...} results corresponding to each reload composition,
                'reload_pos'        : [IASS, IRING, IPOS], reloading assembly position
                'reload_pwr'        : [power at EOC, power after reloading], power of reloaded assembly
                'reload_peak'       : [IASS, IRING, IPOS, INODE, peak_power_density],...}, estimated peak power density after reloading
                'reload_pekass'     : [IASS, IRING, IPOS, peak_assembly_power], estimated peak assembly power after reloading
                'burnup'            : burnup of discharged assembly in this cycle
                'worth'             : [list of string lines of reactivity worth results in reload.out]
                'reload_map'        : [regions to be reloaded with comp_name]
                'search_success'    : True/False, whether a reload position is found
                'next_target_keff'  : estimated new tarted EOC keff for the next cycle
                'next_cycle_length' : estimated cycle length for the next cycle

                The following are only for null refueling case
                'ass2rp'            : {IASS:[IRING,IPOS], ...}  mapping IASS to (IRING, IPOS) for candidate positions
                'null'              : {comp_name:{sub_dict}, ...}
                    each sub_dict contains the following items
                    'peak_power'    : {IASS: peak_assembly, peak_assembly_power}  reload position and corresponding peak assembly power
                    'peak_pwrden'   : {IASS: peak_assembly, peak_node, peak_power_density}   reload position and corresponding peak power density
                    'worth'         : {IASS: [fresh_worth, added_reactivity]} reload position and corresponding reactivity worth
                
        Or empty dict {} when reload search not completed successfully
    '''
    try:
        with open(file,'r') as f:
            line = skip_to(f,'Constraints for reloading position:')
            consider_cycle = 'cycle length' in line
            consider_keff  = consider_cycle or 'reactivity' in line
            null_refueling = 'null_refueling' in line
            results = {}
            results['null_refueling'] = null_refueling
            line = skip_to(f,'number of fresh fuel compositions =')
            num_reload_comp = int(line.split('=')[1])
            reload_comp_names = f.readline().split()
            
            # read peak assembly power and power density at BOC and EOC
            line = skip_to(f,'Peak assembly power at BOC occurs').split(',')
            bocpwr_ass = int(line[0].split()[-1])
            bocpwr_val = float(line[1].split('=')[1])
            line = skip_to(f,'Peak node-averaged power density at BOC occurs').split(',')
            boc_peak_ass = int(line[0].split()[-1])
            boc_peak_nod = int(line[1].split()[-1])
            boc_peak_val = float(line[2].split('=')[1])
            ass2rp = {}
            line = skip_to(f,'*** Assembly powers at the end of cycle (EOC) ***')
            skip_lines(f,2)
            while 1:
                line = f.readline().split()
                if len(line) == 0: break
                iass = int(line[0])
                ass2rp[iass] = [int(line[1]), int(line[2])]
            line = skip_to(f,'Peak assembly power at EOC occurs').split(',')
            eocpwr_ass = int(line[0].split()[-1])
            eocpwr_val = float(line[1].split('=')[1])
            line = skip_to(f,'Peak node-averaged power density at EOC occurs').split(',')
            eoc_peak_ass = int(line[0].split()[-1])
            eoc_peak_nod = int(line[1].split()[-1])
            eoc_peak_val = float(line[2].split('=')[1])
            results['boc_pekass'] = [bocpwr_ass] + ass2rp[bocpwr_ass] + [bocpwr_val]
            results['eoc_pekass'] = [eocpwr_ass] + ass2rp[eocpwr_ass] + [eocpwr_val]
            results['boc_peak'] = [boc_peak_ass] + ass2rp[boc_peak_ass] + [boc_peak_nod, boc_peak_val]
            results['eoc_peak'] = [eoc_peak_ass] + ass2rp[eoc_peak_ass] + [eoc_peak_nod, eoc_peak_val]
            
            # cycle length
            if consider_cycle:
                line = skip_to(f,'cycle length (days) =')
                results['this_cycle_length'] = float(line.split('=')[1])
            else:
                results['this_cycle_length'] = -1.0
            f.seek(0)

            if null_refueling: # read estimated peak powers
                results['null'] = OrderedDict()
                results['ass2rp'] = ass2rp
                for iload in range(num_reload_comp):
                    comp_name = reload_comp_names[iload]
                    results['null'][comp_name] = {}
                    pekpwr, pwrden, dummy = read_estimated_powers(f)
                    results['null'][comp_name]['peak_power'] = pekpwr
                    results['null'][comp_name]['peak_pwrden'] = pwrden                    
                    results['null'][comp_name]['worth'] = read_estimated_worths(f)
                
            else: # normal refueling
                # peak assembly/node averaged powers
                reload_pekass, reload_peak, reload_asspwr = read_estimated_powers(f)
                # peak power limits
                line = skip_to(f,'limit of peak power density (W/cc) =')
                results['peak_power_limit'] = float(line.split('=')[1])
                line = skip_to(f,'limit of assembly-integrated power (W) =')
                results['assembly_power_limit'] = float(line.split('=')[1])
                # assembly worth
                if consider_keff: results['worth'] = read_estimated_worths(f)

            # final result of reload search
            if null_refueling:
                results['search_success'] = True
            else:
                line = skip_to(f,'search for reloading position:')
                results['search_success'] = 'SUCCESS' in line
                if results['search_success']:
                    line = skip_to(f,'reloading position (assembly index):')
                    reload_pos = int(line.split(':')[1])
                    pekpwr_ass, pekpwr_val = reload_pekass[reload_pos]      # estimated peak assembly power after reloading
                    peak_ass, peak_nod, peak_val = reload_peak[reload_pos]  # estimated peak power density after reloading
                    line = skip_to(f,'discharge burnup (MWD/MT) of reloaded assembly')
                    results['burnup'] = float(line.split(':')[1])
                    if consider_cycle:
                        line = skip_to(f,'targeted keff at EOC of next cycle')
                        results['next_target_keff'] = float(line.split(':')[1])
                        line = skip_to(f,'estimated cycle length (days) for next cycle')
                        results['next_cycle_length'] = float(line.split(':')[1])
                    results['reload_pos'] = [reload_pos] + ass2rp[reload_pos]
                    results['reload_pekass'] = [pekpwr_ass] + ass2rp[pekpwr_ass] + [pekpwr_val]
                    results['reload_peak'] = [peak_ass] + ass2rp[peak_ass] + [peak_nod, peak_val]    # estimated peak node position and power density 
                    results['reload_pwr'] = reload_asspwr[reload_pos]
                    # read regions to be modified
                    reload_map = {}
                    skip_to(f,'*** Summary of reloaded assembly')
                    skip_lines(f,3)
                    line = f.readline().split()
                    name = line[1]  # reloaded composition
                    results['reload_comp'] = name
                    reload_map[name] = []
                    for reg in line[5:]:
                        if reg not in reload_map[name]:
                            reload_map[name].extend([reg])
                    results['reload_map'] = reload_map
        return results
    except Exception:
        return {}

# ----------------------------------------------------------------------
def set_manual_reload(sequence,filename):
    ''' put a pair of ring and position numbers for refueling in file filename to do manual reloading
    '''
    if sequence is None:
        if os.path.exists(filename): os.remove(filename)
        return None
    else:
        with open(filename,'w') as f:
            f.write('{}  {}'.format(*sequence[:2]))
        return None if len(sequence)==2 else sequence[2:]

# ----------------------------------------------------------------------
def import_history_burnup(file):
    burnup = {}
    with open(file,'r') as f:
        lines = f.read().splitlines()
        for line in lines[4:]:
            data = line.split()
            burnup[data[0]] = [float(data[3]), float(data[2]), int(data[1])]
    return burnup

# ----------------------------------------------------------------------
def export_burnup_for_reload(file, burnup):
    with open(file,'w') as f:
        line = f'{len(burnup)}\n'
        for reg in burnup:
            #!   region_name     heavy_metal_mass    burnup      burned_cycles
            line += f'  {reg:8s}  {burnup[reg][0]:13.5E}  {burnup[reg][1]:13.5E}  {burnup[reg][2]:3d}\n'
        f.write(line)

# ----------------------------------------------------------------------
def print_reload_summary(output, cycles=None, cycle_length=None, results=None, print_discharge=False):
    from rebus_files import scan_rebus_output
    ''' print summary of reload sequence, keff history, and peak powers in each cycle
    cycles - list of reload cycle indices (starting from 1)
    results - returned data of function reload_output_scan 
    print_discharge - whether to edit nuclide densities of discharged fuel composition
    '''
    if cycles is None:
        cycles = sorted([int(f.split('_')[1]) for f in os.listdir() if 'reload.out_' == f[:11]])
    num_cycle = len(cycles)
    null_refueling = False
    if results is None:  # post-processing
        results = {}
        for icycle in cycles:
            file = f'reload.out_{icycle}'
            print(f'reading reload output {file}')
            results[icycle] = reload_output_scan(file)
            if results[icycle]['null_refueling']: 
                print_discharge = False
                null_refueling = True
            
    if cycle_length is None:
        cycle_length = []
        for icycle in cycles:
            cycle_length.append(results[icycle]['this_cycle_length'])

    keff = []
    peakf= []
    assden = None     # will store assembly averaged atom density of active isotopes
    for i in range(num_cycle):
        icycle = cycles[i]
        print(f'reading rebus output rebus.out_{icycle}')
        if print_discharge:
            discharged = list(results[icycle]['reload_map'].values())[0]        # discharged regions
            rebusout = scan_rebus_output(f'rebus.out_{icycle}', keff=True, peakf=True, discharged=discharged)
            discharged = rebusout['discharged']         # composition of discharged regions
            if assden is None:
                isolist = discharged['active_isotopes']
                niso = len(isolist)
                assden = np.ndarray((num_cycle, niso), dtype=float)
            atdn = np.asarray([0.0]*niso, dtype=float)
            totvol = 0.0
            regvol = discharged['region_vol']   # normally not needed. just in case region volume changes
            for reg in regvol:
                atdn += np.asarray(discharged[reg], dtype=float) * regvol[reg]
                totvol += regvol[reg]
            assden[i] = 1.0/totvol * atdn
        else:
            rebusout = scan_rebus_output(f'rebus.out_{icycle}', keff=True, peakf=True)
        keff.append(rebusout['keff'])
        peakf.append(rebusout['peakf'])

    # ------------------ null refueling summary
    if null_refueling:
        fresh_fuels = [name for name in results[1]['null']]
        candidates = [i for i in results[1]['null'][fresh_fuels[0]]['peak_power']]
        ass2rp = results[1]['ass2rp']
        form1 = '{:6d} {:6d} {:6d}  '
        form2 = '{:13.5f}'
        form3 = '{:13.5e}'
        header = f'\n  candidate position     burned cycles\n'\
        '  IASS  IRING   IPOS '+''.join([f'{i:13d}' for i in cycles])+'\n\n'
        output.write('\n ====== Estimated fresh assembly worth and peak power quantities ======\n')
        for comp in fresh_fuels:
            output.write(f'\n*** {comp}: fresh assembly worth\n')
            output.write(header)
            for ipos in candidates:
                line = form1.format(ipos,*ass2rp[ipos])
                for icycle in cycles:
                    data = results[icycle]['null']
                    line += form2.format(data[comp]['worth'][ipos][0])
                output.write(line+'\n')
            
            output.write(f'\n*** {comp}: reactivity addition\n')
            output.write(header)
            for ipos in candidates:
                line = form1.format(ipos,*ass2rp[ipos])
                for icycle in cycles:
                    data = results[icycle]['null']
                    line += form2.format(data[comp]['worth'][ipos][1])
                output.write(line+'\n')
                
            output.write(f'\n*** {comp}: peak assembly power (W)\n')
            output.write(header)
            for ipos in candidates:
                line = form1.format(ipos,*ass2rp[ipos])
                for icycle in cycles:
                    data = results[icycle]['null']
                    line += form3.format(data[comp]['peak_power'][ipos][1])
                output.write(line+'\n')
            
            output.write(f'\n*** {comp}: peak node-averaged power density (W/cc)\n')
            output.write(header)
            for ipos in candidates:
                line = form1.format(ipos,*ass2rp[ipos])
                for icycle in cycles:
                    data = results[icycle]['null']
                    line += form3.format(data[comp]['peak_pwrden'][ipos][2])
                output.write(line+'\n')
        
        return 0   # return advance 
    
    # ------------------ summary for normal refueling 
    line = '\n reloading positions: (power values are evaluated at BOC of reloaded cycles, peak power density means peak node-averaged power density)\n'\
        ' cycle |  reload position   |  discharge       | peak assembly power (W) with reloaded assembly | peak power density (W/cc) with reloaded assembly\n'\
        ' index |  IASS  IRING  IPOS |  burnup (MWD/MT) |  estimated      calculated     difference (%)  |  estimated      calculated     difference (%)   \n'
    output.write(line)
    form = '{:5d}    {:5d}{:5d}{:5d}     {:13.5E}     {:13.5E}   {:13.5E}       {:5.2f}        {:13.5E}   {:13.5E}       {:5.2f}\n'
    for i in range(1,num_cycle):
        pwr1  = results[i]['reload_pekass'][3]
        pwr2  = results[i+1]['boc_pekass'][3]
        diff1 = (pwr1 / pwr2 - 1.0) * 100
        peak1 = results[i]['reload_peak'][4]
        peak2 = results[i+1]['boc_peak'][4]
        diff2 = (peak1 / peak2 - 1.0) * 100
        output.write(form.format(i, *results[i]['reload_pos'], results[i]['burnup'], pwr1, pwr2, diff1, peak1, peak2, diff2))
    form = '{:5d}    {:5d}{:5d}{:5d}     {:13.5E}     {:13.5E}     N/A                N/A         {:13.5E}     N/A                N/A ! REBUS calculation not performed after last reloading\n'
    i = num_cycle
    output.write(form.format(i,*results[i]['reload_pos'], results[i]['burnup'], results[i]['reload_pekass'][3], results[i]['reload_peak'][4]))
    output.write('note: IRING = hexagon ring number, IPOS = position number\n')

    #output.write('\n depletion history: (from REBUS output)\n cycle       k-effective       peaking factor\n'\
    #' index       boc      eoc       boc      eoc   \n')
    #line = '{:5d}     {:8.5f} {:8.5f}  {:8.5f} {:8.5f}\n'
    #for i in range(num_cycle):
    #    output.write(line.format(i+1, *keff[i], *peakf[i]))

    output.write('\n fuel management history: (negative cycle length means it is not edited)\n'\
    ' burned  cycle        reload     keff   | peak assembly     | peak assembly | peak node position      | peak power    \n'\
    ' cycles  length,day   status            | IASS  IRING  IPOS | power (W)     | IASS  IRING  IPOS  INOD | density (W/cc)\n')
    line = '{:5d}  {:12.4E}   {:6s}  {:8.5f}  {:6d}{:6d}{:6d}  {:13.5E}   {:6d}{:6d}{:6d}{:6d}  {:13.5E}\n'
    i = cycles[0]
    if i > 1:
        output.write(line.format(i-1, 0.0, 'after',keff[0][0], *results[i]['boc_pekass'], *results[i]['boc_peak']))
    else:
        output.write(line.format(i-1, 0.0, '-',keff[0][0], *results[i]['boc_pekass'], *results[i]['boc_peak']))
    for i in range(num_cycle):
        j = cycles[i]
        output.write(line.format(j, cycle_length[i], 'before', keff[i][1], *results[j]['eoc_pekass'], *results[j]['eoc_peak']))
        if i+1 < num_cycle:
            output.write(line.format(j, 0.0, 'after', keff[i+1][0], *results[j+1]['boc_pekass'], *results[j+1]['boc_peak']))
        else:
            output.write(f'{cycles[-1]:5d}  {0.0:12.4E}   after    => REBUS calculation not performed after the last reloading\n')

    if print_discharge:
        output.write(f'\n composition of discharged fuel: (elapsed time starts from BOC of cycle {cycles[0]})\n'\
        ' burned  elapsed      active isotope densities (1/barn-cm)\n'\
        f' cycles  time, day    {"".join([name.ljust(12) for name in isolist])}\n')
        line = '\n{:5d} {:12.4E}  '+'{:12.4E}'*niso
        time = 0.0
        for i in range(num_cycle):
            time += cycle_length[i]
            output.write(line.format(i+1, time, *assden[i]))
    
    output.write('\n\n reloaded assembly power and worth vs burnup:\n'\
        ' burned | reload position | discharg burnup | power of reloaded assembly | assembly worth (keff change in pcm)\n'\
        ' cycles |  IRING    IPOS  |    (MWD/MT)     |     EOC         reloaded   |  discharged    refueled       added\n')
    line = ' {:5d}   {:8d}{:8d}   {:14.5E}   {:13.5E} {:13.5E}  {:12.1f}{:12.1f}{:12.1f}\n'
    for i in range(num_cycle):
        j = cycles[i]
        data = results[j]
        values = data['reload_pos'][1:]
        values.append(data['burnup'])
        values.extend(data['reload_pwr'])
        pos = data['reload_pos'][0]
        worth = data['worth'][pos]
        worth = [worth[0]-worth[1], worth[0], worth[1]]
        values.extend(np.asarray(worth)*1.E+5)
        output.write(line.format(j,*values))

# ----------------------------------------------------------------------
def reload_post(workdir, cycles=None, power_edit=None):
    ''' post processing reload.out of each cycle
        workdir - working directory for REBUS/reload job
        power_edit - None, all, or list of cycle indices for which the assembly power distributions are edited
    '''
    output = open('reload_post.out','w')
    root_dir = os.getcwd()
    os.chdir(workdir)
    if cycles is None: cycles = sorted([int(f.split('_')[1]) for f in os.listdir() if 'reload.out_' == f[:11]])
    print('Discharged assembly compositions are not edited, if needed, turn on print_discharge flag in function reload_post')
    print_reload_summary(output, cycles=cycles, print_discharge=True)
    maxcycle = max(cycles)
    print(f"maxcycle = {maxcycle}")
    if power_edit:
        if power_edit == 'all': power_edit = cycles
        if (not maxcycle in power_edit): power_edit.append(maxcycle)
        output.write('\n\n evolution of assembly power distributions:\n'\
        ' Fuel assembly position    |Assembly powers at BOL and EOC of each burned cycle -->\n'\
        '  No.   IRING   IPOS   ')
        for cycle in power_edit:
            if (cycle>maxcycle): 
                print(f'Skip {cycle} in power_edit (larger than max iteration)')
                continue
            if (cycle==maxcycle):
                output.write(''.join([f'{cycle:8d}(EOL)']))
            else:
                output.write(''.join([f'{" "*5}{cycle:8d}']))
        asspwr = []
        asspos = []
        if not os.path.exists('reload.out_1'):
            print(f'Output reload.out_1 is not found in the working directory, BOL assembly powers cannot be read')
            sys.exit(1)
        with open('reload.out_1','r') as f:
            skip_to(f,'*** Assembly powers at the beginning of cycle (BOC) ***')
            skip_lines(f,3)
            while 1:
                line = f.readline()
                if len(line.strip()) == 0: break
                line = line.split()
                asspos.append([int(val) for val in line[:3]])
                asspwr.append([float(line[3])])
        nass = len(asspos)
        nedit= 0
        for icycle in power_edit:
            if (icycle>maxcycle): 
                # print(f'Skip {icycle} in power_edit (larger than max iteration)')
                continue
            file = f'reload.out_{icycle}'
            if not os.path.exists(file):
                print(f'Output {file} is not found in the working directory')
                sys.exit(1)
            print(f'read file{file}')
            with open(file,'r') as f:
                skip_to(f,'*** Assembly powers at the end of cycle (EOC) ***')
                skip_lines(f,2)
                for iass in range(nass):
                    line = f.readline().split()
                    asspwr[iass].append(int(float(line[3])))
                nedit += 1
        # nedit = len(power_edit)+1
        line = '\n{:5d} {:6d} {:6d}  ' + '{:13.5E}'*nedit
        for iass in range(nass):
            output.write(line.format(*asspos[iass], *asspwr[iass]))
        output.write('\n')

    output.close()
    os.chdir(root_dir)

# ----------------------------------------------------------------------
if __name__ == "__main__":
    nargs = len(sys.argv)
    if nargs < 2:
        print('Going to post-process ReloadPy outputs. Default working directory is the current directory, if not use:')
        print(' ./reload_files.py  rebus_job_directory [cycle range for which assembly powers are summarized]')
        print(' cycle range specification: all or [start, end (included), step].\n'\
        '   if only one positive number is specified, then cycle 1 to that cycle is included.\n')
        workdir = '.'
    else:
        workdir = sys.argv[1]
    if nargs <= 2:
        cycles_for_power_edit = None
    elif nargs == 3:            
        if sys.argv[2].lower() == 'all':
            cycles_for_power_edit = 'all'
        else:
            end = int(sys.argv[2])
            cycles_for_power_edit =[]
            for i in range(1,end+1): cycles_for_power_edit.append(i)
    elif nargs == 4:
        start, end = [int(val) for val in sys.argv[2:]]
        cycles_for_power_edit =[]
        for i in range(start,end+1): cycles_for_power_edit.append(i)
    else:
        start, end, step = [int(val) for val in sys.argv[2:]]
        # cycles_for_power_edit = range(start, end+1, step)
        cycles_for_power_edit =[]
        for i in range(start, end+1, step): cycles_for_power_edit.append(i)
        if (not end in cycles_for_power_edit):
            cycles_for_power_edit.append(end)
    # cycles_for_power_edit = itertools.chain(cycles_for_power_edit, [end+1])
    # print(cycles_for_power_edit)
    # print(cycles_for_power_edit)
    reload_post(workdir, power_edit=cycles_for_power_edit)
