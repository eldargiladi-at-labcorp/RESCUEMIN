from csv import DictReader
from collections import defaultdict
import os
import argparse
import shutil

def filter_high_coverage(minimize_segs,loc2coverage,minimum_coverage,roi):
    seg2remove = set()
    removed_segs = dict()

    for seg in minimize_segs:
        chrom,start,end = minimize_segs[seg]
        minCov = get_min_coverage(loc2coverage,chrom,start,end)
        
        if minCov >= minimum_coverage:
            seg2remove.add(seg)
            
    for seg in seg2remove:
        removed_segs[seg] = minimize_segs[seg]
        del minimize_segs[seg]

    return removed_segs


def eval_seg_info(minimize_segments,coverage_segments,coverage_segs,minimize_segs,removed_segs):

    total_seg_len = 0
    bad_seg_len = 0
    num_fail = 0

    total_segs = len(coverage_segs)+len(minimize_segs)+len(removed_segs)
    num_fail = len(coverage_segments)+len(minimize_segments)
    num_pass = total_segs-num_fail

    for seg in coverage_segs:
        chrom,a,b = coverage_segs[seg]
        total_seg_len += b-a

    for seg in minimize_segs:
        chrom,a,b = minimize_segs[seg]
        total_seg_len += b-a

    for seg in removed_segs:
        chrom,a,b = removed_segs[seg]
        total_seg_len += b-a        

    for seg,mm in coverage_segments:
        chrom,a,b = coverage_segs[seg]
        bad_seg_len += b-a

    for seg,efn in minimize_segments:
        chrom,a,b = minimize_segs[seg]
        bad_seg_len += b-a

    good_seg_len = total_seg_len - bad_seg_len
        
    return good_seg_len,bad_seg_len,num_pass,num_fail



def dump(prog_version,minimize_segments,coverage_segments,config_path,config,outfile,coverage_segs,minimize_segs,total1,loc2coverage):

    status = "FAIL"    
    if total1 < config['maxEFN']:
        status = "SUCCESS"

    with open(outfile,"w") as outfp:
        print("# rescueMinimization version=%s" % prog_version,file=outfp)
        print("# config file path=%s" % config_path,file=outfp)
        for attribute in config:
            print("# %s=%s" % (attribute,config[attribute]),file=outfp)
        print("# exit status for sample: %s" % status,file=outfp)
        if status == "FAIL":
            return 0
        
        for name,efn in minimize_segments:
            chrom,start,end = minimize_segs[name]
            mm = get_min_coverage(loc2coverage,chrom,start,end)

            print("%s\t%d\t%d\t%s\t%d\t%.4e" % (minimize_segs[name][0],minimize_segs[name][1],minimize_segs[name][2],name,mm,efn),file=outfp)
        for name,mm in coverage_segments:
            print("%s\t%d\t%d\t%s\t%d\tNA" % (coverage_segs[name][0],coverage_segs[name][1],coverage_segs[name][2],name,mm),file=outfp)
    return 1


def load_sample(input_dir,coverage_file):
    ret = dict()
    sample = os.path.join(input_dir,coverage_file)

    with open(sample,encoding='windows-1252') as infp:
        for row in DictReader(infp):
            loc = "%s:%d" % (row['Mapping'],int(row['Reference position']))
            cov = int(row['Coverage'])
            ret[loc] = cov
    return ret

def get_args():
    parser = argparse.ArgumentParser(
                    prog = 'rescueMinimizer',
                    description = 'Identifies segments that need to be rescued')

    parser.add_argument('input_directory', type=str,help='input directory of sample')
    parser.add_argument('config_file_path', type=str,help='path to config file')
    parser.add_argument('outfile', type=str,help='Output file for results')
    parser.add_argument('--minimization',type=int,help='1 exome fn negative minimization is used, 0 exome coverage only')
    parser.add_argument('--ROI',type=str,help='bed file for ROI segments')
    parser.add_argument('--detectionCurvesTable',type=str,help='Path to table of detection curves')
    parser.add_argument('--defaultAlleleFrequencies',type=float,help='Default allele frequency')
    parser.add_argument('--VEPpositions',type=str,help='Path to vep positions bed file')
    parser.add_argument('--PLPpositions',type=str,help='Path to plp positions bed file')
    parser.add_argument('--gnomadVariants',type=str,help = 'Path to gnomad variants csv')
    parser.add_argument('--maxEFN',type=float,help = 'Maximum expected number of FN')
    parser.add_argument('--maxRescue',type=int,help = 'Maximum number of rescues per sample')
    parser.add_argument('--minCoverage',type=int,help = 'Minimum coverage criteria')
    parser.add_argument('--removeHighCoverage',type=int,help='1: ignore segments with min coverage >= minCoverage from rescue minimization, 0: otherwise') 
    parser.add_argument('--coverageFile',type=str,default="coverage_table.Targetedv1.0_NotMerged_Coverage_Segment.all.csv",help='Path to csv file with coverage information')
    parser.add_argument('--no_dnr', action='store_true',default=False,help='When set make DNR as coverage')

    args = parser.parse_args()

    # Load the config file


    config = load_config_file(args.config_file_path)

    # Overide attributes

    attributes = ['minimization','ROI','detectionCurvesTable','defaultAlleleFrequencies','VEPpositions','PLPpositions','gnomadVariants','maxEFN','maxRescue','minCoverage','coverageFile','removeHighCoverage','no_dnr']

    for attribute in attributes:
        val = getattr(args,attribute)
        if val is None:
            continue
        config[attribute]=val
    
    if config['no_dnr']:
        config['no_dnr']=1
    else:
        config['no_dnr']=0
    
    # Add the default argument to config

    config['input_directory'] = args.input_directory
    config['outfile'] = args.outfile

    return config,args.config_file_path



def load_variants(name):
    retval = dict()
    with open(name,encoding='windows-1252') as infp:
        for row in DictReader(infp,delimiter='\t'):
            variant = "%s:%s:%s:%s" % (row['CHR'],row['POS'],row['REF'],row['ALT'])
            retval[variant] = [row['GROUP']]
    return retval

def get_loc2segments(segs):
    retval = defaultdict(list)
    for name in segs:
        a=segs[name][1]
        b=segs[name][2]
        chrom = segs[name][0]

        for i in range(a+1,b+1):
            vv="%s:%d" % (chrom,i)
            retval[vv].append(name)
            
    return retval
    

def load_config_file(path):
    ret=dict()
    with open(path,encoding='windows-1252') as infp:
        directory = os.path.dirname(path)
        for line in infp:
            line = line.strip()
            if len(line) == 0:
                continue
            if line[0] == "#":
                continue
            if line.find("#") != -1:
                raise RuntimeError("load_config_file: config file can only have comments in seprate lines, and comment line should begin with #: %s" % line)
            
            ff = line.split("=")
            if len(ff) != 2:
                raise RuntimeError("load_config_file: config file entries should be of the form variable=value: %s" % line)
            ff[0]=ff[0].strip()
            ff[1]=ff[1].strip()

            data_files = ['ROI','detectionCurvesTable','VEPpositions','gnomadVariants','PLPpositions']

            if ff[0] in data_files:
                ret[ff[0]] = os.path.join(directory,ff[1])
            else:
                ret[ff[0]] = ff[1]

    attributes = ['minimization','ROI','detectionCurvesTable','defaultAlleleFrequencies','VEPpositions','PLPpositions','gnomadVariants','maxEFN','maxRescue','minCoverage','removeHighCoverage','no_dnr']

    missing = list()
    for x in attributes:
        if not x in ret:
            missing.append(x)
    if len(missing) > 0:
        raise RuntimeError('load_config_file: the following attributes are missing from config file: %s' % ",".join(missing))

    ret['minimization'] = int(ret['minimization'])
    ret['maxEFN'] = float(ret['maxEFN'])
    ret['maxRescue'] = int(ret['maxRescue'])
    ret['minCoverage'] = int(ret['minCoverage'])
    ret['defaultAlleleFrequencies'] = float(ret['defaultAlleleFrequencies'])
    ret['removeHighCoverage'] = int(ret['removeHighCoverage'])
    ret['no_dnr'] = int(ret['no_dnr'])
            
    return ret


def rescue2(total,lll,max_res,epsilon):
    total1 = total
    num = 0
    segments = set()
    while total1 > epsilon and num < max_res:
        lll = sorted(lll,key=lambda entry: sum([x[1] for x in entry[1]]),reverse=True)

        set_first = lll[0][1]

        weight_first = sum([x[1] for x in set_first])
        
        total1 -= weight_first
        segments.add((lll[0][0],weight_first))
        num+=1
                           
        if total1 < epsilon:
            break

        lll = lll[1:]
        for i in range(len(lll)):
            lll[i][1] -= set_first
    
    if total1 < epsilon:
        return num,total1,segments
    return 0,total,set([])


def process(loc2coverage,location2variants,group_cov2prob,roi,loc2segments,epsilon,max_res):
    ll = set(list(location2variants.keys()))
    
    segment2loc_total = defaultdict(set)
    total = 0
    
    for loc in loc2coverage:
        cov = loc2coverage[loc]
        if cov > 999:
            cov = 999

        if not loc in roi:
            continue

        if not loc in ll:
            continue
        
        if not loc in loc2segments:
            continue
        
        loc_total = 0
        for group,af in location2variants[loc]:
            pp = 1-group_cov2prob[group][cov]
            loc_total += af*pp
        total += loc_total
        
        segments = loc2segments[loc]
        for segment in segments:
            segment2loc_total[segment].add((loc,loc_total))

    lll = list(segment2loc_total.items())
    lll = [[x[0],x[1]] for x in lll]
    lll = sorted(lll,key=lambda entry: sum([x[1] for x in entry[1]]),reverse=True)

    num=0
    total1=total
    segments_to_remove = []
    
    if total > epsilon:
        num,total1,segments_to_remove =rescue2(total,lll,max_res,epsilon)

    return segments_to_remove, total,total1,num

def get_min_coverage(loc2coverage,chrom,a,b):
    mm = 100000000  # Very large number

    for i in range(a+1,b+1):  
        loc = "%s:%d" % (chrom,i)
        if not loc in loc2coverage:
            raise RuntimeError("get_min_coverage: location %s not in sample: inconsistent data" % loc)
            ### AD HOC TEMPORARY FIX !!!!! ###
            #print("WARNING: get_min_coverage: location %s not in sample: inconsistent data" % loc)
            continue
        cov = loc2coverage[loc]
        mm=min(mm,cov)

    return mm

def process1(coverage_segs,loc2coverage,min_coverage):
    segments = list()
    for seg in coverage_segs:
        chrom=coverage_segs[seg][0]
        a=coverage_segs[seg][1]
        b=coverage_segs[seg][2]

        mm = get_min_coverage(loc2coverage,chrom,a,b)

        if mm < min_coverage:
            segments.append((seg,mm))

    return segments

def loc2vars(variant2info):
    ret = defaultdict(list)
    for var in variant2info:
        ff = var.split(":")
        loc = "%s:%s" %(ff[0],ff[1])
        ret[loc].append(variant2info[var])
    return ret

def merge(plp2type,vep2type):
    for var in vep2type:
        if var in plp2type:
            continue
        plp2type[var] = vep2type[var]
    return plp2type

def add_af(var2type,variant2gnomad,default_af):
    for vv in var2type:
        if vv in variant2gnomad:
            var2type[vv].append(variant2gnomad[vv])
        else:
            var2type[vv].append(default_af)

def add_af1(var2type,variant2gnomad,default_af):
    for vv in var2type:
        ff = vv.split(":")
        if ff[0] == "X":
            if vv in variant2gnomad:
                var2type[vv].append(variant2gnomad[vv]*2)
            else:
                var2type[vv].append(default_af*2)
        else:
            if vv in variant2gnomad:
                var2type[vv].append(variant2gnomad[vv])
            else:
                var2type[vv].append(default_af)
            

def load_gnomad(name):
    ret = dict()
    with open(name,encoding='windows-1252') as infp:
        for row in DictReader(infp):
            var = row['ID']
            af = float(row['AF_ALL'])
            if var in ret:
                if af!=ret[var]:
                    raise RuntimeError('load_gnomad: inconsistent value for AF')
            ret[var] = af
    return ret

        

def load_group_cov_prob(name):
    retval = defaultdict(list)

    with open(name,encoding='windows-1252') as infp:
        for row in DictReader(infp):
            group = row['group']
            coverage = int(row['coverage'])
            prob = float(row['prob'])

            retval[group].append(prob)

    return retval


def load_roi(name,no_dnr):
    roi=set() # Position that are treated with rescue minimization (not coverage)
    bad_pos=set() # Positions to ignore in rescue minimization
    coverage_segs = dict()
    minimize_segs = dict()
    with open(name,encoding='windows-1252') as infp:
        for line in infp:
            ff = line.strip().split()
            a=int(ff[1])
            b=int(ff[2])
            chrom = ff[0].replace("chr","")
            name=ff[3]
            status = ff[4]

            if (status == 'Ignore' and no_dnr == 0) or (status == 'AlwaysIgnore'):
                for i in range(a+1,b+1):
                    vv="%s:%d" % (chrom,i)
                    bad_pos.add(vv)
            elif status == 'Minimize':
                for i in range(a+1,b+1):
                    vv="%s:%d" % (chrom,i)
                    roi.add(vv)
                minimize_segs[name]=(chrom,a,b)
            elif status == 'Coverage' or (status == 'Ignore' and no_dnr == 1):
                for i in range(a+1,b+1):
                    vv="%s:%d" % (chrom,i)
                    bad_pos.add(vv)
                coverage_segs[name]=(chrom,a,b)
            else:
                raise RuntimeError('load_roi: unexpected status for segment: %s' % status)

    roi = roi.difference(bad_pos)

    return roi,coverage_segs,minimize_segs

def load_roi_no_minimize(name,no_dnr):
    roi=set() # Position that are treated with rescue minimization (not coverage)
    coverage_segs = dict()
    minimize_segs = dict()
    with open(name,encoding='windows-1252') as infp:
        for line in infp:
            ff = line.strip().split()
            a=int(ff[1])
            b=int(ff[2])
            chrom = ff[0].replace("chr","")
            name=ff[3]
            status = ff[4]

            if (status == 'Ignore' and no_dnr == 0) or (status == 'AlwaysIgnore'):
                pass
            elif status == 'Minimize' or status == 'Coverage' or (status == 'Ignore' and no_dnr == 1):
                coverage_segs[name]=(chrom,a,b)
            else:
                raise RuntimeError('load_roi: unexpected status for segment: %s' % status)

    return roi,coverage_segs,minimize_segs


def main():
    prog_version="2.5.2"
    config,config_path=get_args()
    
    if config['minimization'] == 1:
        roi,coverage_segs,minimize_segs = load_roi(config['ROI'],config['no_dnr'])
    elif config['minimization'] == 0:
        roi,coverage_segs,minimize_segs = load_roi_no_minimize(config['ROI'],config['no_dnr'])
    else:
        raise RuntimeError('--minimization should be 0,1 but it is %d' % config['minimization'])

    default_af=config['defaultAlleleFrequencies']

    # get the coverage information

    input_directory = config['input_directory']
    loc2coverage = load_sample(input_directory,config['coverageFile'])

    # Filter minimize segments with enough coverage

    removed_segs = dict()
    
    if config['removeHighCoverage'] == 1:
        removed_segs = filter_high_coverage(minimize_segs,loc2coverage,config['minCoverage'],roi)
    
    loc2minimize_segments = get_loc2segments(minimize_segs)

    # load group coverage 2 prob

    group_cov2prob = load_group_cov_prob(config['detectionCurvesTable'])

    plp2type = load_variants(config['PLPpositions'])
    vep2type = load_variants(config['VEPpositions'])

    # load plp gnomad and vep gnomad

    variant2gnomad = load_gnomad(config['gnomadVariants'])

    # add allele frequency

    add_af1(plp2type,variant2gnomad,default_af)
    add_af1(vep2type,variant2gnomad,default_af)

    # merge the two groups

    variant2info = merge(plp2type,vep2type)

    # generate location to variant structure

    location2variants = loc2vars(variant2info)

    epsilon = config['maxEFN']
    max_res = config['maxRescue']

    # Do the FN minimization
    minimize_segments,total,total1,num = process(loc2coverage,location2variants,group_cov2prob,roi,loc2minimize_segments,epsilon,max_res)

    # Do the min coverage analysis
    coverage_segments = process1(coverage_segs,loc2coverage,config['minCoverage'])

    # Output the segments that require rescue in a bed file

    dump(prog_version,minimize_segments,coverage_segments,config_path,config,config["outfile"],coverage_segs,minimize_segs,total1,loc2coverage)


if __name__ == "__main__":
    try:
        main()
    except OSError as error:
        print("FAILURE: rescueMinimization terminated abnormally with an OS exception")
        print(error)
        exit(1)
    except RuntimeError as error:
        print("FAILURE: rescueMinimization terminated abnormally with a runtime exception")
        print(error)
        exit(1)
    except Exception as error:
        print("FAILURE: rescueMinimization terminated abnormally with a unknown exception")
        print(error)
        exit(1)

    print("SUCCESS: rescueMinimization completed successfully")
    exit(0)
