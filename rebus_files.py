import os
from text_io import skip_lines, skip_to, read_next_line, copy_lines

# ----------------------------------------------------------------------
def scan_file_header(f):
    data = f.readline().split()
    maxcard = int(data[1])
    freeform = int(data[3])
    nline = maxcard // 16
    if nline*16 < maxcard: nline += 1
    ncards = []
    for line in range(nline):
        ncards.extend([int(i) for i in f.readline().split()])
    return (data[0].strip(), maxcard, freeform, ncards)

# ----------------------------------------------------------------------
def read_isotope_mapping(aburn):
    ''' read active isotope mapping, return {isotope:[list of equivalent labels]}
        aburn = A.BURN file name
    '''
    with open(aburn,'r') as f:
        setname, maxcard, freeform, ncards = scan_file_header(f)
        # incomplete
        
# ----------------------------------------------------------------------
def read_burnup(file, info=None, fluence=False):
    ''' burnup = read_burnup(file, info)
        info = (num_burn_step, num_active_isotope)
        if fluence = False:
            burnup = {region_name:[heavy_metal_mass, burnup, burned_cycles], ....}
        else:
            burnup = {region_name:[heavy_metal_mass, burnup, burned_cycles, fast_fluence], ....}
        here fast_fluence is peak fast fluence in each region
        Note that peak fast fluence cannot be acumulated through multiple cycles directly
    '''
    if info:
        num_burn_step, num_active_isotope = info
    else: # scan input copy to get basic info
        with open('ABURN','r') as fin:
            setname, maxcard, freeform, ncards = scan_file_header(fin)
            nline = sum(ncards[:2])
            skip_lines(fin,nline)
            if freeform:
                line = fin.readline().split()
                num_burn_step = int(line[6])
            else:
                line = fin.readline()
                num_burn_step = int(line[60:66])
            nline = sum(ncards[3:8])
            skip_lines(fin,nline)
            active_isotopes = []
            for i in range(ncards[8]):
                line = fin.readline()
                name = line.split()[1] if freeform else line[6:12].strip()
                if name not in active_isotopes:
                    active_isotopes.extend([name])
            num_active_isotope = len(active_isotopes)
    # start to read rebus output
    with open(file,'r') as f:
        skip_to(f,f'REACTOR CONDITIONS AFTER')
        skip_to(f,'MASSES (IN KG) OF ACTIVE ISOTOPES IN EACH REGION')
        nskip = (num_active_isotope+8)//9
        # read initial heavy metal mass of each region
        burnup = {}
        while 1:
            line = read_next_line(f)
            if 'INITIAL TOTAL MASS OF FISSIONABLE ISOTOPES IN THIS REGION =' in line:
                mass = float(line.split('=')[1])*1e-3       # convert to MT unit
                i = 0
                while i < nskip:
                    line = read_next_line(f)
                    if 'FCC004' in line and ' PAGE' in line:
                        continue
                    else:
                        i += 1
                regnam = read_next_line(f).split()[0]
                burnup[regnam] = [mass]
            elif 'TOTAL REACTOR LOADING (IN KG) OF ACTIVE ISOTOPES AT TIME NODE' in line:
                break
        num_active_region = len(burnup)
        # read burnups
        skip_to(f,'START OF CUMULATIVE EDITS')
        for i in range(num_burn_step+1):
            skip_to(f,'CUMULATIVE BURNUP AFTER')
        skip_to(f,'AVERAGE BURNUP (MWD/MT) OF EACH REGION')
        nread = (num_active_region+8)//9
        i = 0
        while i < nread:
            line = read_next_line(f)
            if 'FCC004' in line and ' PAGE' in line:
                skip_to(f,'AVERAGE BURNUP (MWD/MT) OF EACH REGION')
                continue
            elif 'AVERAGE BURNUP (MWD/MT) OF EACH REGION' in line:
                continue
            else:
                line = line.split()
                if line[0] == 'REGION':
                    regs = line[1:]
                    continue
                else:
                    i += 1
                    for j in range(len(regs)):
                        burnup[regs[j]].extend([float(line[j]), 1])

        # read peak burnups and peak fast fluence
        skip_to(f,'CUMULATIVE PEAK BURNUP AND FAST FLUENCE AFTER')
        if fluence:
            skip_to(f,'PEAK FAST FLUENCE (N/CM**2) OF EACH REGION')
            i = 0
            while i < nread:
                line = read_next_line(f)
                if 'FCC004' in line and ' PAGE' in line:
                    skip_to(f,'PEAK FAST FLUENCE (N/CM**2) OF EACH REGION')
                    continue
                elif 'PEAK FAST FLUENCE (N/CM**2) OF EACH REGION' in line:
                    continue
                else:
                    line = line.split()
                    if line[0] == 'REGION':
                        regs = line[1:]
                        continue
                    else:
                        i += 1
                        for j in range(len(regs)):
                            burnup[regs[j]].extend([float(line[j])])
    return burnup

# ----------------------------------------------------------------------
def update_burnup(file, burn0, reload_regions=None, info=None, readfluence=False):
    ''' burn0 - burnup data of previous cycles {region_name:[heavy_metal_mass, burnup, burned_cycles, fast_fluence], ...}
        burnup unit: MWD/MT,  mass unit: MT
        reload_regions: list of previously reloaded fuel region in the last cycle
        info = [num_burn_step, num_active_isotopes]
    '''
    burnup = read_burnup(file, info=info, fluence=readfluence)    # burnups in this cycle
    if readfluence: print('Warning: peak fast fluence read from REBUS one-cycle calculation will not be accumulated')
    if burn0:   # not fresh core
        for reg in burnup:
            if reload_regions and reg in reload_regions: continue
            # add burnups in previous cycles for un-refueled regions
            if reg not in burn0:
                raise ValueError(f'active region {reg} is not presented in previous cycle')
            total_mwd = burn0[reg][0] * burn0[reg][1] + burnup[reg][0] * burnup[reg][1]
            burnup[reg][0] = burn0[reg][0]      # keep the same initial heavy metal mass unless it is reloaded
            burnup[reg][1] = total_mwd / burn0[reg][0]    # accumulated burnup in all cycles up to now
            burnup[reg][2] += burn0[reg][2]
            #burnup[reg][3] += burn0[reg][3]
    return burnup

# ----------------------------------------------------------------------
def print_burnup(file, burnup, editfluence=False):
    with open(file,'w') as f:
        if editfluence:
            line = ' edit of accumulated burnup and peak fast fluence by region:\n'\
                ' region   burned      burnup     heavy metal    peak fast fluence\n'\
                '  name    cycles     (MWD/MT)     mass (MT)      in this cycle (n/cm^2)\n\n'
            form = ' {:6s}   {:4d}  {:14.5E}{:14.5E} {:14.5E}\n'
            for reg in burnup:
                line += form.format(reg, burnup[reg][2], burnup[reg][1], burnup[reg][0], burnup[reg][3])
        else:
            line = ' edit of accumulated burnup by region:\n'\
                ' region   burned      burnup     heavy metal\n'\
                '  name    cycles     (MWD/MT)     mass (MT)\n\n'
            form = ' {:6s}   {:4d}  {:14.5E}{:14.5E}\n'
            for reg in burnup:
                line += form.format(reg, burnup[reg][2], burnup[reg][1], burnup[reg][0])
        f.write(line)

# ----------------------------------------------------------------------
def write_anip14(output, comp_name, comp_data, freeform):
    ''' comp_data is a two-level list [[isotope_name, atom_density], ...]'''
    if freeform:
        head = f'14      {comp_name:6s}'
    else:
        head = f'14          {comp_name:6s}'
    i = 0
    line = head
    if freeform:
        for iso in comp_data:
            if i == 3:
                line += '\n'+head+f' {iso[0]:6s}{iso[1]:12.5E}'
                i = 1
            else:
                line += f' {iso[0]:6s}{iso[1]:12.5E}'
                i += 1
    else:
        for iso in comp_data:
            if i == 3:
                line += '\n'+head+f'{iso[0]:6s}{iso[1]:12.5E}'
                i = 1
            else:
                line += f'{iso[0]:6s}{iso[1]:12.5E}'
                i += 1
    output.write(line+'\n')

# ----------------------------------------------------------------------
def modify_anip1415(nt14out, output, reload_map, reload_composition, freeform):
    ''' replace composition (type 14) in reloaded regions with reloaded composition
        note that composition name is same as region name in ANIP1415 generated with mk1415.x 
        fixed format is used
        PS: since REBUS requires subzones to be defined, we just define a subzone for each modified zone
            because subzones are not really used in the targeted non-equilibrium depletion problem
    '''
    with open(nt14out,'r') as f:
        lines = f.read().splitlines()
    modified = []
    reg2comp = {}
    for comp in reload_map:
        for reg in reload_map[comp]:
            if reg in reg2comp.values():
                raise ValueError(f'multiple compositions are reloaded into region {reg}')
            reg2comp[reg] = comp
        # make subzones using reloaded compositions
        write_anip14(output, comp, reload_composition[comp], freeform)       
    
    for line in lines:
        card = line[:2]
        if card == '14':
            if freeform:
                reg = line.split()[1]
            else:
                reg = line[12:18].strip()       # composition name is same as region name
            if reg in reg2comp:
                if reg not in modified:
                    modified.extend([reg])
                    comp = reg2comp[reg]    # composition name of reloaded assembly
                    write_anip14(output, reg, [[f'{comp}', 1.0]], freeform)
            else:
                output.write(line+'\n')
        else:
            output.write(line+'\n')

# ----------------------------------------------------------------------
def parse_anip14_card(line, freeform):
    ''' reads one A.NIP3 type 14 card and return [['atom_name', num_density],...]
    '''
    if freeform:
        card = line.split()[2:]
        comp = []
        for i in range(0, len(card), 2):
            comp.append([card[i], float(card[i+1])])
    else:
        card = line[18:]
        comp = []
        while len(card.strip()) > 0:
            comp.append([card[:6], float(card[6:18])])
            card = card[18:]
    return comp
    
# ----------------------------------------------------------------------
def copy_anip1415(nt14out, output, a_fuel_zone, freeform):
    ''' copy nt14out (ANIP1415) to rebus input. Since REBUS requires subzones to be defined, 
        we just define a subzone for the given fuel zone. Subzones are not really used in the 
        targeted non-equilibrium depletion problem
    '''
    with open(nt14out,'r') as f:
        lines = f.read().splitlines()
    write_anip14(output, a_fuel_zone, [['SUBZN', 1.0]], freeform)
    for line in lines:
        card = line[:2]
        if card == '14':
            if freeform:
                reg = line.split()[1]
            else:
                reg = line[12:18].strip()       # composition name is same as region name
            if reg == a_fuel_zone:
                comp = parse_anip14_card(line, freeform)
                write_anip14(output, 'SUBZN', comp, freeform)
            else:
                output.write(line+'\n')
        else:
            output.write(line+'\n')

# ----------------------------------------------------------------------
def write_aburn35(output,regions,freeform):
    ''' regions contains all the fuel regions whose composition name is same as region name
    '''
    if len(regions) > 99999:
        raise ValueError('too many active zones')
    i = 0
    if freeform:
        line = '35    P{:05d} {:6s} {:6s}   1   1\n'
    else:
        line = '35    P{:05d}{:6s}{:6s}     1     1\n'
    for reg in regions:
        i += 1
        output.write(line.format(i,reg,reg))

# ----------------------------------------------------------------------
def change_aburn35(filename, active_regions):
    ''' change type 35 cards of A.BURN such that fuel paths are specified with active regions
    '''
    tmpfile = open('tempfile','w')
    with open(filename,'r') as f:
        line = f.readline()
        line = line[:13] + '    0' + line[18:]
        tmpfile.write(line)
        data = line.split()
        maxcard = int(data[1])
        freeform = int(data[3])
        if maxcard < 35:
            raise ValueError('change_aburn35 error: no type 35 card is found in A.BURN')
        nhead = maxcard // 16
        if nhead*16 < maxcard: nhead += 1
        ncards = []
        for k in range(nhead):
            ncards.extend([int(i) for i in f.readline().split()])
        if ncards[34] == 0:
            raise ValueError('change_aburn35 error: only type 35 card is supported for fuel management path')
        num35 = ncards[34]
        ncards[34] = len(active_regions)
        for k in range(nhead):
            i = k*16; j = i+16
            if j > maxcard: j = maxcard
            line = ''.join([f'{i:9d}' for i in ncards[i:j]])
            tmpfile.write(line+'\n')
        ncards[34] = num35

        nline = sum(ncards[:34])
        copy_lines(f, tmpfile, nline)
        # make type 35 cards with region wise paths
        skip_lines(f,ncards[34])
        write_aburn35(tmpfile, active_regions, freeform)
        # other inputs
        nline = sum(ncards[35:])
        copy_lines(f, tmpfile, nline)
    tmpfile.close()
    os.rename('tempfile', filename)

# ----------------------------------------------------------------------
def scan_rebus_output(filename, keff=True, peakf=False, discharged=False):
    ''' keff  = True / False, whether to read keff
        peakf = True / False, whether to read peaking factor
        discharged = None / False / list of discharged compositions whose atom densities need to be retrieved
    '''
    results = {}
    with open(filename,'r') as f:
        if discharged:
            results['discharged'] = {}
            regvol = {}
            skip_to(f,'X PROBLEM CHARACTERISTICS X')
            line = skip_to(f,'NUMBER OF SUB-INTERVALS FOR THE BURN CYCLE')
            nsub = int(line.split()[-1])  # number of subintervals
            # read region volumes 
            skip_to(f,'K-EFFECTIVE =')
            skip_to(f,'REGION     ZONE  ZONE        VOLUME    INTEGRATION(1)')
            skip_lines(f,2)
            while 1:
                line = read_next_line(f)
                if 'DIF3D' in line and 'PAGE' in line:
                    skip_to(f,'REGION     ZONE  ZONE        VOLUME    INTEGRATION(1)')
                    skip_lines(f,2)
                elif '(1) INTEGRATION WEIGHT FACTOR' in line or '(2) THE PEAK POWER DENSITY' in line:
                    continue
                else:
                    line = line.split()
                    if line[0] == 'TOTALS': break
                    reg = line[1]
                    if reg in discharged:
                        regvol[reg] = float(line[4])
            results['discharged']['region_vol'] = regvol
            # move to the last step
            while 1:
                skip_to(f,'***** FCC004 HAS BEEN ENTERED *****')
                line = skip_to(f,'REACTOR CONDITIONS AFTER')
                step = int(line.split('AFTER')[1].split('BURNUP')[0])
                if step == nsub: break
            # read atom densities over regions
            niso = 0
            isolist = []
            nleft = len(discharged)       # number of requested regions
            while 1:
                if niso == 0:
                    skip_to(f,'ATOM DENSITIES (IN ATOMS/BARN-CM.) OF ACTIVE ISOTOPES IN EACH REGION')
                    skip_lines(f,2)
                    while 1:
                        line = f.readline()
                        if 'FCC004' in line and 'PAGE' in line:
                            skip_to(f,'ATOM DENSITIES (IN ATOMS/BARN-CM.) OF ACTIVE ISOTOPES IN EACH REGION')
                            skip_lines(f,2)
                            continue
                        if len(line.strip()) == 0: break
                        isolist.extend(line.split())
                    niso = len(isolist)
                    nskip = (niso+8)//9
                    results['discharged']['active_isotopes'] = isolist

                line = read_next_line(f)
                if 'FCC004' in line and 'PAGE' in line:
                    skip_to(f,'ATOM DENSITIES (IN ATOMS/BARN-CM.) OF ACTIVE ISOTOPES IN EACH REGION')
                    continue
                elif line.strip() == 'REGION':
                    skip_lines(f,nskip)
                    continue 
                else:
                    line = line.split()
                    region = line[0]
                    if region in discharged:
                        data = [float(val) for val in line[1:]]
                        for i in range(nskip-1):
                            data.extend([float(val) for val in f.readline().split()])
                        if len(data) == niso: 
                            nleft -= 1
                            results['discharged'][region] = data
                        if nleft == 0: break
                    else:
                        skip_lines(f,nskip-1)
                        continue
        if keff:
            skip_to(f,'REACTOR CHARACTERISTICS SUMMARY')
            line = skip_to(f,'KEFF')
            results['keff'] = [float(val) for val in line.split()[1:]]
        if peakf:
            if not keff:
                skip_to(f,'REACTOR CHARACTERISTICS SUMMARY')
            line = skip_to(f,'PEAKING FACTOR')
            results['peakf'] = [float(val) for val in line.split()[2:]]
    return results

# ----------------------------------------------------------------------
def read_cycle_length(aburn):
    with open(aburn,'r') as f:
        setname, maxcard, freeform, ncards = scan_file_header(f)
        nline = sum(ncards[:2])
        skip_lines(f,nline)
        line = f.readline()
        if freeform:
            cycle_length = float(line.split()[4])
        else:
            cycle_length = float(line[36:48].strip())
    return cycle_length

# ----------------------------------------------------------------------
def modify_cycle_length(aburn03, cycle_length, freeform):
    ''' rewrite type 03 card of a.burn with new cycle_length
    '''
    if freeform:
        newcard = aburn03.split()       # will strip linefeed
        newcard[4] = f'{cycle_length:11.4E}'.strip()
        newcard = ' '.join(newcard)+'\n'
        if(len(newcard) > 72):
            raise ValueError('card length exceeded 72 column after modifying cycle length')
    else:
        newcard = aburn03
        newcard = newcard[:36] + f'{cycle_length:12.4E}' + newcard[48:]
    return newcard

# ----------------------------------------------------------------------
def fetch_active_regions(anip, aburn):
    ''' fetch active regions from original ANIP3 and ABURN
    '''
    active_regions = []
    with open(aburn,'r') as f:
        fpath = []
        setname, maxcard, freeform, ncards = scan_file_header(f)
        nline = sum(ncards[:34])
        skip_lines(f,nline)
        if freeform:
            for i in range(ncards[34]):
                fpath.extend([ f.readline().split()[2] ])
        else:
            for i in range(ncards[34]):
                fpath.extend([ f.readline()[12:18].strip() ])
    with open(anip,'r') as f:
        setname, maxcard, freeform, ncards = scan_file_header(f)
        nline = sum(ncards[:13])
        skip_lines(f,nline)
        marker = []
        if freeform:
            for i in range(ncards[13]):
                data = f.readline().split()
                if data[1] in fpath or data[2] in fpath: marker.extend([data[1]])
            for i in range(ncards[14]):
                data = f.readline().split()
                if data[1] in marker: active_regions.extend(data[2:])
        else:
            for i in range(ncards[13]):
                line = f.readline()
                zone = line[12:18].strip()
                subz = line[18:24].strip()
                if zone in fpath or subz in fpath: marker.extend([zone])
            for i in range(ncards[14]):
                line = f.readline().strip()
                zone = line[6:12].strip()
                if zone in marker:
                    line = line[12:]
                    leng = len(line)
                    active_regions.extend([line[j:j+6].strip() for j in range(0, leng, 6)])
    return active_regions

# ----------------------------------------------------------------------
def fetch_anip3_info(file):
    ''' get more information from file (ANIP3) generated by DIF3D/REBUS
    '''
    info = {}
    with open(file,'r') as f:
        setname, maxcard, freeform, ncards = scan_file_header(f)
        nline = sum(ncards[:2])
        skip_lines(f,nline)
        igeom = int(f.readline().split()[1])
        info['ndim'] = 2 if igeom in [110,114,116] else 3
        if igeom == 110 or igeom == 120:
            info['nsect'] = 6
        elif igeom == 114 or igeom == 124:
            info['nsect'] = 1
        elif igeom == 116 or igeom == 126:
            info['nsect'] = 2
        else:
            info['nsect'] = 0       # not considered

        if ncards[28] > 0:
            nline = sum(ncards[3:28])
            skip_lines(f,nline)
            info['hex_flat2flat'] = float(f.readline().split()[1])
    return info

# ----------------------------------------------------------------------
def make_rebus_input(filename, change_cycle, cycle_length=None, isRestart=False, mk1415_exe=None, \
    reload_map=None, reload_comp=None, null_refueling=False, initial_flux=None):
    ''' prepare complete rebus input using existing interface files listed in _rebus_files
    filename        - filename for rebus input
    change_cycle    - True or False, whether to change cycle length
    cycle_length    - new cycle length, needed if change_cycle = True
    isRestart       - True or False (default). 
                      True means to make rebus input using BOC composition and new cycle length (during cycle length iteration)
                      False means to make rebus input using EOC and reloaded composition and new cycle length
    mk1415_exe      - path to mk1415 utility code, needed if isRestart = False
    reload_map      - composition to region map for reloaded fuels
    reload_comp     - reloaded fuel compositions / a fuel region name when null_refueling = True
    null_refueling  - if True, will not modify composition
    initial_flux    - if True, will use NHFLUX/RTFLUX as initial flux guess
    '''
    _rebus_files = ['ASUMMAR', 'ASTP027', 'ADIF3D', 'AHMG4C', 'ANIP3', 'ABURN']
    _optional_files = ['ASUMMAR', 'AHMG4C']
    subname = 'make_rebus_input'

    if change_cycle and not cycle_length:
        print(f'{subname} error: cycle_length is not specified')
        return False
    if not isRestart:
        if not mk1415_exe:
            print(f'{subname} error: path to mk1415 utility (mk1415_exe) is not specified')
            return False
        if not null_refueling:
            if not reload_comp or (not reload_map):
                print(f'{subname} error: unspecified reloaded composition')
                return False
    rebusin = open(filename,'w')
    if initial_flux is None:
        rebusin.write('BLOCK=OLD\nDATASET=ISOTXS\nDATASET=NHFLUX\nBLOCK=STP027\n')
    else:
        initial_flux = initial_flux.upper()
        if initial_flux in ['RTFLUX', 'NHFLUX']:
            rebusin.write(f'BLOCK=OLD\nDATASET=ISOTXS\nDATASET={initial_flux}\nBLOCK=STP027\n')
        else:
            print(f'{subname} error: unknow dataset {initial_flux} for initial flux guess')
            return False

    for file in _rebus_files:
        if os.path.exists(file):
            with open(file,'r') as f:
                setname, maxcard, freeform, ncards = scan_file_header(f)
                if freeform:
                    rebusin.write(f'UNFORM={setname}\n')
                else:   # fixed format
                    rebusin.write(f'DATASET={setname}\n')

                if file == 'ANIP3' and (not isRestart):  # need to modify compositions
                    nline = sum(ncards[:12])
                    copy_lines(f, rebusin, nline)
                    # skip 13, 14, 15 cards
                    nline = sum(ncards[12:15])
                    skip_lines(f,nline)
                    # get atom densities after depletion and make new type 14 cards
                    for file in ['NDXSRF', 'ZNATDN', 'LABELS']:
                        if not os.path.exists(file):
                            print(f'{subname} error: required dataset {file} does not exist in the working directory')
                            return False
                    if freeform:
                        os.system(f'{mk1415_exe} y')        # will produce new type 14, 15 cards in ANIP1415
                    else:
                        os.system(mk1415_exe)
                    if null_refueling:
                        copy_anip1415('ANIP1415', rebusin, reload_comp, freeform)    # add one subzone 
                    else:
                        modify_anip1415('ANIP1415', rebusin, reload_map, reload_comp, freeform) # modifiy type 14 cards
                    # other inputs
                    nline = sum(ncards[15:])
                    copy_lines(f, rebusin, nline)
                    
                elif file == 'ABURN' and change_cycle:  # need to modify cycle length
                    if ncards[2] > 1:
                        print(f'{subname} error: more than one type 03 cards present in A.BURN')
                        return False
                    nline = sum(ncards[:2])
                    copy_lines(f, rebusin, nline)
                    line = f.readline()
                    line = modify_cycle_length(line, cycle_length, freeform)
                    rebusin.write(line)
                    nline = sum(ncards[3:])
                    copy_lines(f, rebusin, nline)
                else:
                    nline = sum(ncards)
                    copy_lines(f, rebusin, nline)
        else:
            if file not in _optional_files:
                print(f'{subname} error: file {file} not found')
                return False
    rebusin.close()
    return True

# ----------------------------------------------------------------------
def save_rebus_input(inpfile, icycle):
    ''' save rebus input file for restart run, need to remove DATASET=NHFLUX in OLD block '''
    with open(inpfile,'r') as f:
        lines = [line for line in f if line.strip() != 'DATASET=NHFLUX']
    with open(f'{inpfile}_{icycle}','w') as f:
        for line in lines: f.write(line)
    
# ----------------------------------------------------------------------
def rebus_iterate_cycle_length(output, inpfile, outfile, user_eoc_keff, min_cycle_length, eps_keoc, eps_cycle, rebus_exe, deprate_ratio=None):
    errmsg = '\n       Error in <iterate_cycle_length>:{:s}'
    new_cycle_length = None
    low_cycle_limit = min_cycle_length * (1.0 - eps_cycle)
    while True:
        os.system(f'{rebus_exe} < {inpfile} > {outfile}')
        # read keffs from rebus output 
        boc_keff, eoc_keff = scan_rebus_output(outfile,keff=True)['keff']
        # get cycle length
        if new_cycle_length:    # update with new estimate
            cycle_length = new_cycle_length
        else:
            if deprate_ratio:
                target_keff = boc_keff - deprate_ratio * (boc_keff - user_eoc_keff)
            else:
                target_keff = user_eoc_keff
            cycle_length = read_cycle_length('ABURN')
            if target_keff - boc_keff > 1.0e-5:   # check this only once
                output.write(errmsg.format('targeted keff at eoc > calculated keff at boc'))
                return False
        # check convergence
        if cycle_length == low_cycle_limit and eoc_keff - target_keff < -eps_keoc:
            output.write(errmsg.format('cannot achieve target eoc keff with minimum cycle length because of too fast reactivity loss'))
            return False
        speed = (boc_keff-eoc_keff) / cycle_length
        if eps_keoc < 0.0: eps_keoc = eps_cycle * min_cycle_length * speed
        if abs(eoc_keff - target_keff) > eps_keoc or cycle_length < low_cycle_limit:  # restart with a new cycle length
            # we roll back even when the guessed cycle length is too short because it does not increase effort 
            # and save time on making input with new compositions and accumulating burnups
            new_cycle_length = (boc_keff - target_keff) / speed
            if new_cycle_length < low_cycle_limit:
                new_cycle_length = low_cycle_limit  # if estimated cycle_length < minimum cycle length, set it to minimum cycle length
            output.write(f'\n       target_keff ={target_keff:8.5f}  eoc_keff ={eoc_keff:8.5f}   diff ={eoc_keff-target_keff:8.5f}   tolerance ={eps_keoc:8.5f}')
            output.write(f'\n       min_cycle_length ={min_cycle_length:8.2f}  cycle_length ={cycle_length:8.2f}  tolerance ={eps_cycle*min_cycle_length:8.2f}')
            output.write(f'\n       restart REBUS job with a new cycle length of {new_cycle_length:13.5E} days')
            done = make_rebus_input(inpfile, True, cycle_length=new_cycle_length, isRestart=True, initial_flux='NHFLUX')
            if not done:
                output.write(errmsg.format('failed to modify cycle length of rebus input'))
                return False
        else:   # converged 
            break
    return cycle_length







