#!/usr/bin/env python3
import sys
import os
import shutil
import numpy as np
from time import monotonic
from text_io import *
from rebus_files import *
from reload_files import *

# ------------------------------------------------------------------------------
def run_shell(command, stop_on_error=True):
    ''' execute shell command '''
    import subprocess
    args = command.split()
    output = subprocess.run(args, capture_output=True)
    if output.returncode != 0:
        print('shell returncode:',output.returncode)
        print('error message:',output.stdout.decode('utf-8'))
        if stop_on_error: raise RuntimeError('shell command failed')
    return (output.returncode, output.stdout.decode('utf-8')) 

def print_timestamp(output):
    from datetime import datetime
    from socket import gethostname
    now = datetime.now()
    date = now.date()
    now = now.strftime("%H:%M:%S")
    output.write(f'Simulation initiated on {gethostname()} at {date} {now}\n')

def get_random_composition(options, rnum):
    # import numpy as np
    # rnum = np.random.normal(options['ufrac_mean'], options['ufrac_std'], 1)
    fi = np.array(options['fi_0'])*(1.-rnum) + np.array(options['fi_1'])*rnum
    # fi = fi * (1.-rnum[0])
    mass_sum  = np.dot(fi, options['amass_salt'])
    conv_f = (options['salt_frac']*options['salt_dens']/mass_sum)
    random_comp = []
    for nd in options['fixednd']:
        random_comp.append(nd)
    for i in range(0, options['niso_salt']):
        random_comp[i] += conv_f*fi[i]

    return(random_comp)

# ------------------------------------------------------------------------------
total_elapsed = monotonic()
if len(sys.argv) > 1:
    input_file = sys.argv[1]
else:
    input_file  = 'reload.inp'      # default input file name
if not os.path.exists(input_file):
    raise RuntimeError(f'input file {input_file} not found')

# default working directory
root_dir  = os.getcwd()         # all inputs and executables are put here by default
dif3d_dir = 'work_dif3d'        # working directory for DIF3D run (inside main working directory)
persent_dir = 'work_persent'    # working directory for PERSENT run (inside main working directory)

# read user input
print(f'read input options from {input_file} ...')
options, reload_composition = reload_input_scan(input_file)
# set up main working directory
work_dir = 'run_'+options['case_name']
# total number of fuel cycles to be modeled, the first cycle starts as it is, the other cycles start with reloaded fuel
num_cycle = options['num_cycles']
first_cycle = options['num_past_cycles'] + 1
isRestart = first_cycle > 1
if first_cycle > num_cycle:
    print(f'required number of depletion cycles = {num_cycle}\n'\
        f'number of previous depletion cycles = {first_cycle-1}')
    print('done with all depletion cycles already, restart execution stopped')
    quit()
if isRestart:
    output = open(f'reload_{options["case_name"]}.out_restart','w')
    if not os.path.exists(work_dir):
        raise ValueError(f'old working directory {work_dir} does not exist')
else:
    output = open(f'reload_{options["case_name"]}.out','w')
    if os.path.exists(work_dir): run_shell(f'rm -r {work_dir}')
    os.mkdir(work_dir)
    if os.path.exists('rebus.burnup_0'): shutil.move('rebus.burnup_0',work_dir)
    if os.path.exists('reload.fluence_0'): shutil.move('reload.fluence_0',work_dir)
    # user provided reference power distribution
    if options['reload_refpwr'] == 'user':
        if os.path.exists('reload.refpwr'): 
            shutil.copy2('reload.refpwr',f'{work_dir}/reload.refpwr')
        else:
            raise ValueError('user provided reference power distribution file reload.refpwr not found')
print_timestamp(output)
print(f'check input options and REBUS model ...')
reload_input_check(input_file, options, output)
output.flush()

if (options['random_enri']=='T'):
    if (len(options['fixednd'])!=len(reload_composition['RFUEL'])):
        raise ValueError('User provided different number of NDs in RFUEL and fixednd')
    base_composition = get_random_composition(options, options['ufrac_mean'])

# executable and input files
rebus_exe, dif3d_exe, persent_exe = options['rebus_exe'], options['dif3d_exe'], options['persent_exe']
reload_exe = options['reload_utility']  # utility code to determine assembly position for fuel reloading
mk1415_exe = options['type14_utility']  # utility code to make type 14 and 15 cards of A.NIP3 using REBUS interface datasets
rebus_input = options['rebus_input']
xsdata_file = options['isotxs_file']

# global options
adjust_cycle = options['adjust_cycle_length'] == 'T'
adjust_keoc = adjust_cycle and options['deprate_ratio'] is not None

# relative tolerance for cycle length adjustment
eps_cycle = options['eps_cycle']
# absolute tolerance for EOC keff comparison
eps_keoc = options['eps_keoc']
manual_sequence = options['initial_sequence']
if (options['initial_sequence']):
    if (first_cycle > 0):
        manual_sequence = manual_sequence[(first_cycle-1)*2:]
manual_cycle = options['initial_cycle']
null_reload = options['null_refueling'] == 'T'

# --------------------------------------------------------------------
#   work starts here
print('start REBUS calculation')
os.chdir(work_dir)
if '/' not in rebus_exe:  rebus_exe  = f'../{rebus_exe}'
if '/' not in reload_exe: reload_exe = f'../{reload_exe}'
if '/' not in mk1415_exe: mk1415_exe = f'../{mk1415_exe}'
if options['use_persent'] == 'T':
    if '/' not in dif3d_exe: dif3d_exe = f'../{dif3d_exe}'
    if '/' not in persent_exe: persent_exe = f'../{persent_exe}'
    if not os.path.exists('dif3d.x'): run_shell(f'ln -s {dif3d_exe} dif3d.x')
    if not os.path.exists('persent.x'): run_shell(f'ln -s {persent_exe} persent.x')
elif options['dif3d_check'] == 'T':
    if '/' not in dif3d_exe: dif3d_exe = f'../{dif3d_exe}'
    run_shell(f'ln -s {dif3d_exe} dif3d.x')

# initialize storage of results
results = {}
burnup = {}
reload_regions = None
cycle_length = []

# in case of restart run, load in previous results
if isRestart:
    # read reload.out for previous cycles
    for i in range(1,first_cycle):
        results[i] = reload_output_scan(f'reload.out_{i}')
        print(f'{i}')
        cycle_length.append(results[i]['this_cycle_length'])
    if not null_reload:
        # mark previous refueling positions
        if options['disperse_reload']=='T':
            icycle = first_cycle-1
            with open('reload.previous','w') as f:
                i = options['disperse_numpast']
                if i == 0: i = icycle           # check all previous positions
                j = icycle if i > icycle else i
                f.write(f'{j}  {options["disperse_mindist"]:13.5E}   ! number of previously reloaded assemblies and minimum distance\n')
                i = 1 if i > icycle else icycle+1-i
                line = ' '.join([f'{results[j]["reload_pos"][0]} ' for j in range(i,icycle+1)])
                f.write(line+'  ! positions of previously reloaded assemblies')
        # get reloaded regions in last reloading
        reload_regions = []
        reload_map = results[first_cycle-1]['reload_map']
        for comp in reload_map:
            reload_regions.extend(reload_map[comp])
    # and assembly burnups so far
    restart_burnup = f'rebus.burnup_{options["num_past_cycles"]}'
    if not os.path.exists(restart_burnup):
        raise ValueError(f'burnup file {restart_burnup} generated from previous depletion cycles is not found')
    burnup = import_history_burnup(restart_burnup)
    # reset fast fluence history to the current step
    shutil.copy2(f'reload.fluence_{options["num_past_cycles"]}', 'reload.fluence')
    # clean old links
    for item in ['user.NHFLUX','user.NAFLUX','init.NAFLUX']:
        if os.path.exists(item): os.remove(item)
    if os.path.exists(persent_dir): shutil.rmtree(persent_dir)

else:   # check initial burnup/fluence files (may be provided if the job not started from fresh core)
    if os.path.exists('rebus.burnup_0'): burnup = import_history_burnup('rebus.burnup_0')
    if os.path.exists('reload.fluence_0'): os.rename('reload.fluence_0','reload.fluence')

if rebus_input:  # override saved input in working directory in restart case
    shutil.copy2(f'{root_dir}/{rebus_input}','rebus.in')
elif isRestart:  # using saved input from last run
    shutil.copy2(f'rebus.in_{first_cycle}','rebus.in')
else:
    raise ValueError('required REBUS input file is not provided')
shutil.copy2(f'{root_dir}/{xsdata_file}','ISOTXS')


output.write('\n *** program execution footprints:\n\n')
# initial run of rebus
output.write(f'--- cycle {first_cycle:4d}:\n------ perform REBUS calculation for depletion')
elapsed = monotonic()
if adjust_cycle:
    if (manual_cycle):
        target_keff = 1.0
    else:
        target_keff = options['target_eoc_keff']
    min_cycle_length = options['min_cycle_length']
    this_cycle_length = rebus_iterate_cycle_length(output, 'rebus.in', 'rebus.out', target_keff, min_cycle_length, eps_keoc, eps_cycle, rebus_exe)
    if not this_cycle_length: quit()
else:
    #run_shell(f'{rebus_exe} < rebus.in > rebus.out')   # why it is not working?
    os.system(f'{rebus_exe} < rebus.in > rebus.out')
elapsed = monotonic() - elapsed
output.write(f'\n       finished in{" "*54}{elapsed:13.5e} seconds\n')
if adjust_cycle:
    cycle_length.append(this_cycle_length)
else:
    cycle_length.append(read_cycle_length('ABURN'))

#if fetch_anip3_info('ANIP3')['ndim'] == 3:
#    print('Warning: deprate_ratio is neglected when 3D core model is used')
#    output.write('Warning: deprate_ratio is neglected when 3D core model is used. target_eoc_keff is not adjusted\n')
#    options['deprate_ratio'] = None
#    adjust_keoc = False
active_regions = fetch_active_regions('ANIP3','ABURN')
change_aburn35('ABURN',active_regions)

# prepare input for utility code
if os.path.exists('RTFLUX'): 
    options['flux_data'] = 'RTFLUX'
elif os.path.exists('NHFLUX'):
    options['flux_data'] = 'NHFLUX'
else:
    output.flush()
    raise RuntimeError('neither RTFLUX nor NHFLUX is generated by DIF3D/REBUS, flux solution cannot be obtained')
options['fuel_regions'] = active_regions
options['num_fuel_region'] = len(active_regions)
shutil.copy2(f'{root_dir}/{input_file}','reload.comp')
prepare_utility_input('reload',options)

# tags of intermediate files
tobe_renamed = ['ZNATDN_0000_','NDXSRF_0000_','LABELS_0000_','PWDINT_0000_','RTFLUX_0000_', 'NHFLX0_0000_']
tobe_removed = ['GEODST_','LABELS_','RFILES_','ZNATDN_','PWDINT_','RTFLUX_','RZFLUX_','NHFLX0_']
# successive cycles with assembly reloading
finished_all_cycles = False
for icycle in range(first_cycle, num_cycle+1):
    print(f'finished depletion cycle {icycle}')
    # save BOC data and clean useless interface files, EOC data are saved by default
    datasets = next(os.walk('.'), (None, None, []))[2]
    files1 = [file for file in datasets if file[:12] in tobe_renamed]
    files2 = [file for file in datasets if file not in files1 and file[:7] in tobe_removed]
    files2.extend(['EOCRF','POINTR','STACK'])
    for file in files1: os.rename(file, f'boc_{file[:6]}')
    for file in files2: os.remove(file)
    if options['use_VARIANT'] and options['VARIANT_Pn_order']==1:   # NHFLUX not saved at EOC
        os.remove('NHFLUX')     # this is the NHFLUX at BOC
        run_shell('ln -s NHFLX0  NHFLUX')   # in P1 case, use NHFLX0 as NHFLUX at EOC
    # update region burnups by reading rebus output
    output.write('------ read REBUS output for burnup\n')
    output.flush()
    burnup = update_burnup('rebus.out', burnup, reload_regions)
    export_burnup_for_reload(options['burnup_data'], burnup)
    print_burnup(f'rebus.burnup_{icycle}', burnup)
    # use manual reloading sequence if any
    manual_sequence = set_manual_reload(manual_sequence,'reload.manual')
    
    # perform reloading search
    output.write('------ find assembly to be reloaded\n')
    output.write('--------- run reload search utility program in')
    if os.path.exists(dif3d_dir): shutil.rmtree(dif3d_dir)   # avoid warning in reload.x
    if os.path.exists(persent_dir): shutil.rmtree(persent_dir)  # avoid warning in reload.x
    elapsed = monotonic()
    ierr, errmsg = run_shell(reload_exe, stop_on_error=False)
    if ierr != 0:
        output.write(f'\n{" "*10}utility program aborted\n')
        break
    elapsed = monotonic() - elapsed
    output.write(f'{" "*26}{elapsed:13.5e} seconds\n')
    # read reload.out for editing
    output.write('--------- read utility output in')
    elapsed = monotonic()
    results[icycle] = reload_output_scan('reload.out')
    if results[icycle] == {}:
        output.write('\n--------- failed to read reload output, probably due to termination of utility program. Check reload.out in working directory for error message')
        break
    if not results[icycle]['search_success']:
        output.write('\n--------- search for reloading position failed, program terminated\n')
        break
    os.rename('reload.out',f'reload.out_{icycle}')
    shutil.copy2('reload.fluence', f'reload.fluence_{icycle}')
    # get peak_power_limit from first cycle if necessary
    if icycle == 1:
        if not null_reload and options['fixed_peak_power'] == 'T':
            if options['peak_power_limit'] < 0.0 or options['assembly_power_limit'] < 0.0:
                options['peak_power_limit'] = results[icycle]['peak_power_limit']
                options['assembly_power_limit'] = results[icycle]['assembly_power_limit']
                prepare_utility_input('reload',options)     # put peak powers in first cycle as power limits for later cycles
    elapsed = monotonic() - elapsed
    output.write(f'{" "*40}{elapsed:13.5e} seconds\n')
    # save rebus input/output files of this cycle
    save_rebus_input('rebus.in',icycle)
    os.rename('rebus.out',f'rebus.out_{icycle}')
    
    # update positions of recently reloaded assemblies
    if options['disperse_reload']=='T':
        with open('reload.previous','w') as f:
            i = options['disperse_numpast']
            if i == 0: i = icycle           # check all previous positions
            j = icycle if i > icycle else i
            f.write(f'{j}  {options["disperse_mindist"]:13.5E}   ! number of previously reloaded assemblies and minimum distance\n')
            i = 1 if i > icycle else icycle+1-i
            line = ' '.join([f'{results[j]["reload_pos"][0]} ' for j in range(i,icycle+1)])
            f.write(line+'  ! positions of previously reloaded assemblies')
    # mark just reloaded regions for burnup accumulation
    reload_regions = []
    if not null_reload:
        reload_map = results[icycle]['reload_map']
        for comp in reload_map:
            reload_regions.extend(reload_map[comp])
    # make new rebus input
    if (options['random_enri']=='T'):
        # update composition before generating rebus input
        rnum = np.random.normal(options['ufrac_mean'], options['ufrac_std'], 1)
        output.write(f'U fraction sampling {rnum[0]*100:.6f} % at {icycle} - th cycles\n')
        temp_composition = get_random_composition(options, rnum[0])
        for i in range(0, len(temp_composition)):
            reload_composition['RFUEL'][i][1] = temp_composition[i]

    output.write('------ prepare new REBUS input with reloaded assembly in')
    next_cycle_length = results[icycle]['next_cycle_length'] if adjust_cycle else None
    if (manual_cycle):
        if (icycle <= len(options['initial_cycle'])):
            next_cycle_length = options['initial_cycle'][icycle-1]
    elapsed = monotonic()
    flux_guess = None
    if options['use_VARIANT']: flux_guess = 'NHFLUX'
    if null_reload:
        done = make_rebus_input('rebus.in', False, mk1415_exe=mk1415_exe, reload_comp=active_regions[0], \
            null_refueling=True, initial_flux=flux_guess)
    else:
        done = make_rebus_input('rebus.in', adjust_cycle, cycle_length=next_cycle_length, mk1415_exe=mk1415_exe, \
            reload_map=reload_map, reload_comp=reload_composition, initial_flux=flux_guess)
    if (options['random_enri']=='T'):
        # retreive composition after generating rebus input
        for i in range(0, len(base_composition)):
            reload_composition['RFUEL'][i][1] = base_composition[i]
    if not done:
        output.flush()
        raise RuntimeError('failed to make new rebus input')
    elapsed = monotonic() - elapsed
    output.write(f'{" "*16}{elapsed:13.5e} seconds\n')
    output.flush()

    if icycle == num_cycle:
        save_rebus_input('rebus.in',icycle+1)
        finished_all_cycles = True
        break
    # run rebus for next cycle
    output.write(f'--- cycle {icycle+1:4d}:\n------ perform REBUS calculation for depletion')
    elapsed = monotonic()
    
    if adjust_cycle and not(manual_cycle):
        new_cycle_length = rebus_iterate_cycle_length(output, 'rebus.in', 'rebus.out', target_keff, min_cycle_length, \
            eps_keoc, eps_cycle, rebus_exe, deprate_ratio=options['deprate_ratio'])
        if not new_cycle_length: break
        cycle_length.append(new_cycle_length)
    else:
        #run_shell(f'{rebus_exe} < rebus.in > rebus.out')
        os.system(f'{rebus_exe} < rebus.in > rebus.out')
        if (manual_cycle):
            cycle_length.append(next_cycle_length)
        else:
            cycle_length.append(cycle_length[icycle-1])
    elapsed = monotonic() - elapsed
    output.write(f'\n       finished in{" "*54}{elapsed:13.5e} seconds\n')

# cumulative summary edit over all cycles
if finished_all_cycles:
    output.write('--- finished all fuel cycles\n')
else:
    output.write(f'\n--- !!! job terminated in cycle {icycle} !!!\n')
    num_cycle = icycle-1
total_elapsed = monotonic() - total_elapsed
output.write(f'--- total execution time (s) = {total_elapsed:13.5e}\n')

output.write('\n *** summary edits (for more information, refer to reload.out_xxx output files in the working directory)\n')

cycles = range(1,num_cycle+1)
if null_reload:
    print_reload_summary(output, cycles=cycles, results=results)
else:
    print_reload_summary(output, cycles=cycles, cycle_length=cycle_length, results=results)

output.close()

