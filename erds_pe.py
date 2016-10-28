from __future__ import division
import os
import argparse
import glob
import jxcnv_functions as jf
from DataManager import *
from hmm.Model import *
from hmm.ModelParams import *
import operator 
import numpy as np
from VCFReader import *
from ParameterEstimation import *
import fileinput
import pdb
import sys
import time

def bamlist2RPKM(args):
    MAQ = 20 #TODO: set the MAQ as a input parameter.
    print "MAQ threshold:",MAQ

    try:
        import pysam
    except:
        print 'Cannot load pysam module!'
        sys.exit(0)
    try:
        # read target
	target_fn = str(args.target)
	targets = jf.loadTargets(target_fn)
	num_target = len(targets)
    except IOError:
        print 'Cannot read target file: ', target_fn
        sys.exit(0)

    if not os.path.exists(args.output):
        os.mkdir(args.output)

    try:
        bamlist_f = open(args.input)
    except IOError:
        sys.exit(0)

    for line in bamlist_f.readlines():
        line = line.strip('\n')
        temp = line.split('\t')
        sample_name = temp[0] + '.' + temp[2]
        bam_file = line.split('\t')[1]
        
        f = pysam.AlignmentFile(bam_file, 'rb')

        if not f.has_index():
            print 'No index found for ', bam_file
            sys.exit(0)
    
        readcount = np.zeros(num_target)
        exon_bp = np.zeros(num_target)
        targetIDs = np.zeros(num_target)

        # detect contig naming scheme here # TODO, add an optional "contigs.txt" file or automatically handle contig naming
        bam_contigs = f.references
        targets_contigs = [str(t) for t in set(map(operator.itemgetter("chr"), targets))]
                    
        targets2contigmap = {}
        for targets_contig in targets_contigs:
            if targets_contig in bam_contigs:
                targets2contigmap[targets_contig] = targets_contig
            elif jf.chrInt2Str(targets_contig) in bam_contigs:
                targets2contigmap[targets_contig] = jf.chrInt2Str(targets_contig)
            elif jf.chrInt2Str(targets_contig).replace("chr","") in bam_contigs:
                targets2contigmap[targets_contig] = jf.chrInt2Str(targets_contig).replace("chr","")
            else:
                print "[ERROR] Could not find contig '%s' from %s in bam file! \n[ERROR] Perhaps the contig names for the targets are incompatible with the bam file ('chr1' vs. '1'), or unsupported contig naming is used?" % (targets_contig, target_fn)
                sys.exit(0)

        print 'Calculating RC and RPKM values...'
        
        total_reads = 0
        counter = 0
        for t in targets:
            t_chr = targets2contigmap[str(t['chr'])]

            t_start = t['start']
            t_stop = t['stop']
            
            try:
                iter = f.fetch(t_chr, t_start, t_stop)
            except:
                print "[ERROR] Could not retrieve mappings for region %s:%d-%d. Check that contigs are named correctly and the bam file is properly indexed" % (t_chr,t_start,t_stop)
                sys.exit(0)

            for i in iter:
                if i.pos+1 >= t_start and i.mapq >= MAQ:
                    readcount[counter] += 1
                    total_reads += 1

            exon_bp[counter] = t_stop - t_start + 1
            targetIDs[counter] = counter + 1
            counter += 1

        print 'Found %d reads in the target regions of bam file with MAQ >= %d: ' %(total_reads, MAQ), bam_file
        # calculate RPKM values for all targets 
        rpkm = (10**9*(readcount)/(exon_bp))/(total_reads)
        rpkm_f = open(args.output+'/'+sample_name+'.rc.rpkm', 'w')
        
        rpkm_f.write('chr\tstart\tstop\tRC\tRPKM\n')
        for i in range(len(rpkm)):
            rpkm_f.write(targets[i]['chr'] + '\t' + str(targets[i]['start']) + '\t' + str(targets[i]['stop']) + '\t' + str(readcount[i]) + '\t' + str(rpkm[i]) + '\n')
        rpkm_f.close()

    bamlist_f.close()
    
def RPKM2Matrix(args):
    #Renjie modified
    rpkm_dir = str(args.rpkm_dir)
    rpkm_files = glob.glob(rpkm_dir + "/*")
    if len(rpkm_files) == 0:
        print 'Cannot find any rpkm files'
        sys.exit(0)
    if not os.path.exists(args.output):
        os.mkdir(args.output)
        print 'Output dir created: ',args.output

    try:
        # read target
        target_fn = str(args.target)
        targets = jf.loadTargetsStr(target_fn)
        num_target = len(targets)
    except IOError:
        print 'Cannot read target file: ', target_fn
        sys.exit(0)
   
    samples = {}
    for f in rpkm_files:
        s = '.'.join(f.split('/')[-1].split('.')[0:-1])
        samples[s] = f

    RPKM_matrix = np.zeros([num_target, len(samples)], dtype=np.float)
    RC_matrix = np.zeros([num_target, len(samples)], dtype=np.float)
    for i,s in enumerate(samples.keys()):
        rc = np.loadtxt(samples[s], dtype=np.float, delimiter="\t", skiprows=1, usecols=[3])
        rpkm = np.loadtxt(samples[s], dtype=np.float, delimiter="\t", skiprows=1, usecols=[4])
        RC_matrix[:,i] = rc
        RPKM_matrix[:,i] = rpkm
        print "Successfully read RC and RPKM for " + s
    
    output_rpkm = '/RPKM_matrix.rpkm.raw'
    output_rc = '/RC_matrix.rc.raw'

    if args.output:
        output_rpkm = args.output+'/RPKM_matrix.raw'
        output_rc = args.output+'/RC_matrix.raw'
    output_rc_f = open(output_rc, 'w')
    output_rpkm_f = open(output_rpkm, 'w')

    output_rc_f.write('Targets\t' + '\t'.join(samples.keys()) + '\n')
    for i in range(len(RC_matrix)):
        output_rc_f.write(targets[i] + '\t' + '\t'.join(str(r) for r in RC_matrix[i]) + '\n')
    output_rc_f.close()

    output_rpkm_f.write('Targets\t' + '\t'.join(samples.keys()) + '\n')
    for i in range(len(RPKM_matrix)):
        output_rpkm_f.write(targets[i] + '\t' + '\t'.join(str(r) for r in RPKM_matrix[i]) + '\n')
    output_rpkm_f.close()
                   
def svd(args):
    filename = args.rpkm_matrix 
    f_dir = os.path.dirname(filename)
    if f_dir != '':
        f_dir = f_dir + '/'

    output = filename
    if args.output:
        output = f_dir + str(args.output)
   
    # count the columns number of the data file
    f = open(filename)
    temp = f.readline().strip().split('\t')
    colsnum = len(temp)
    
    # skip 1st row and 4 columns
    print 'Loading file...'
    data = np.loadtxt(filename, dtype=np.float, delimiter='\t', skiprows=1, usecols=range(4, colsnum)) 
    # loading targets str
    targets = jf.loadTargetsStrFromFirstCol(filename)
    # names of samples
    samples = temp[4:] 

    print 'SVD...'
    U, S, Vt = np.linalg.svd(data, full_matrices=False)
    
    index = S < 0.7 * np.mean(S)
    new_S = np.diag(S * index)
    
    # reconstruct data matrix
    data_new = np.dot(U, np.dot(new_S, Vt))

    # save to files
    file_u = open(output + '.U', 'w')
    file_s = open(output + '.S', 'w')
    file_vt = open(output + '.Vt', 'w')
    
    print 'Saving SVD files...'
    np.savetxt(file_u, U, delimiter='\t')
    np.savetxt(file_s, S, delimiter='\t')
    np.savetxt(file_vt, Vt, delimiter='\t')
    file_u.close()
    file_s.close()
    file_vt.close()

    print 'Saving matrix..'
    jf.saveRPKMMatrix(output + '.SVD', samples, targets, np.transpose(data_new))

def discover(args) :
    paramsfile = args.params
    sample_req = args.sample
    hetsnp = args.hetsnp
    tagsnp = args.tagsnp
    vcf_file = args.vcf
    mode = 'SVD'
    
    if hetsnp == 'True' or hetsnp == 'TRUE':
        hetsnp = True
    else:
        hetsnp = False
    
    if tagsnp == 'True' or tagsnp == 'TRUE':
        tagsnp = True
    else:
        tagsnp = False

    datafile = args.datafile 
    f_dir = os.path.dirname(datafile)
    if f_dir != '':
        f_dir = f_dir + '/'

    if args.output:
		outputfile = f_dir + str(args.output)

    tagsnp_file = args.tagsnp_file

    sample_flag = False #used to check whether sample_req exists

    print 'Loading data file...',
    dataloader = DataManager(datafile)
    print 'Done!'
    print 'Loading paramters...',
    params = dataloader.getParams(paramsfile)
    print 'Done!'
    dataloader.skipHeadline()
    sample = dataloader.getNextSample()

    targets_list = dataloader.getTargetsList()
    output_aux = file(outputfile+'.aux', 'w')
    output_aux.write('SAMPLE_ID\tCNV_TYPE\tFULL_INTERVAL\tINDEX\tINTERVAL\tREAD_DEPTH\n')
    output = file(outputfile,'w')
    output.write('SAMPLE_ID\tCNV_TYPE\tINTERVAL\tCHROMOSOME\tSTART\tSTOP\tLENGTH\n')

    if (hetsnp or tagsnp) and vcf_file == '':
        print 'Error: please indicate a vcf file!'
        system.exit(0)

    if vcf_file != '':
        vcf_reader = VCFReader(vcf_file)
    else:
	vcf_reader = False

    if tagsnp:
        print 'Loading tagSNP information ...',
        cnp_dict = vcf_reader.loadTagSNP(tagsnp_file)
        print 'Done!'

    while sample :
        if sample_req == '' or (sample_req != '' and sample['sample_id'] == sample_req):
            sample_flag = True
            print time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) ,sample_req,'......'

            #Renjie added: To check whether the VCF contains sample_req.
            if vcf_file != '':
	            vcf_checker = vcf.Reader(open(vcf_file,'r'))
	            if sample['sample_id'] in vcf_checker.samples:
	                sample_in_VCF = True
	            elif sample_req in vcf_checker.samples:
	                sample_in_VCF = True
	            else:
	                print 'No sample %s in VCF file.'%sample_req
	                sample_in_VCF = False

	            if hetsnp and sample_in_VCF :
	                print 'Parsing SNV information from VCF file for: ' + sample['sample_id']
	                snp_info = vcf_reader.getSNPInfo(sample['sample_id'], targets_list)

	            if tagsnp and sample_in_VCF:
	                print 'Analysing tagSNP information from tagSNP database for: ' + sample['sample_id'],
	                cnp_list = vcf_reader.findTagSNPForSample(sample['sample_pop'], sample['sample_id'], cnp_dict)
	                tagsnp_info_list = vcf_reader.findExonWithTagSNP(cnp_list, targets_list, overlap_threshold=0.5)
	                print len(tagsnp_info_list)


            #estimate NB paramters from sample['observations']  
            sample_observations = []
            remove_list = []
            sample['observations'] = [ float(x) for x in sample['observations']]
            
            #slicing: target_index is used to split observations sequence
            target_index_begin = 0
            target_index_end = 0
            temp = 1

            sample_observations_list = []
            snp_info_list = []

            for i, targets in enumerate(targets_list):
                target_index_end = target_index_begin + len(targets)
                if hetsnp and sample_in_VCF:
                    snp_info_list.append(snp_info[target_index_begin:target_index_end])
                sample_observations_list.append(sample['observations'][target_index_begin:target_index_end])

                target_index_begin = target_index_end

            for i in range(len(sample_observations_list)):
                sample_observations_list[i] = ndarray.tolist(stats.zscore(sample_observations_list[i]))
                   
            for i, targets in enumerate(targets_list):
                print 'Running HMM for sample[' + sample['sample_id'] + ']: ',
                print 'chr' + targets[0]._chr + ' [' + str(temp) + '|' + str(len(targets_list)) + ']'
                temp += 1
		
                #Run the HMM 
                if not hetsnp and not tagsnp:
                    modelParams = ModelParams(mode, params, targets, het_nums=0, tagsnp=0)
                elif sample_in_VCF and hetsnp and not tagsnp:
                	modelParams = ModelParams(mode, params, targets, snp_info_list[i], tagsnp=0)
                elif sample_in_VCF and not hetsnp and tagsnp:
                	modelParams = ModelParams(mode, params, targets, het_nums=0, tagsnp=tagsnp_info_list[i])
                elif sample_in_VCF and hetsnp and tagsnp:
                	modelParams = ModelParams(mode, params, targets, snp_info_list[i], tagsnp_info_list[i])
                elif not sample_in_VCF and hetsnp and tagsnp:
                    modelParams = ModelParams(mode, params, targets, het_nums=0, tagsnp=0)
                else:
                    pdb.set_trace()
	
                model = Model(mode, modelParams, sample_observations_list[i])
                pathlist = list()
                
                if vcf_reader and sample_in_VCF:
                    pathlist = model.forwardBackward_Viterbi(mode, if_snp = True)
                else:
                    pathlist = model.forwardBackward_Viterbi(mode, if_snp = False)
                dataloader.outputCNVaux(output_aux, sample['sample_id'], targets, pathlist, sample_observations_list[i])
                dataloader.outputCNV(output, sample['sample_id'], targets, pathlist, sample_observations_list[i])

        sample = dataloader.getNextSample()

    output.close()
    output_aux.close()
    dataloader.closeFile()

    if not sample_flag:
        print 'Could not find the sample_id specified.'

parser = argparse.ArgumentParser(prog='ERDS-pe', description='Designed by RT.')
subparsers = parser.add_subparsers()

#BAM List -> RPKM
svd_parser = subparsers.add_parser('rpkm', help="Create RPKM matrix from a BAM list")
svd_parser.add_argument('--target', required=True, help='Target definition file')
svd_parser.add_argument('--input', required=True, help='BAM file list, each line for each sample')
svd_parser.add_argument('--output', required=True, help='Directory for RPKM files')
svd_parser.set_defaults(func=bamlist2RPKM)

#RPKM files -> Matrix
svd_parser = subparsers.add_parser('merge_rpkm', help="Merge RPKM files to a matrix")
svd_parser.add_argument('--rpkm_dir', required=True, help='RPKM files')
svd_parser.add_argument('--target', required=True, help='Target definition file')
svd_parser.add_argument('--output', required=False, help='Matrix file')
svd_parser.set_defaults(func=RPKM2Matrix)

#SVD
svd_parser = subparsers.add_parser('svd', help="SVD")
svd_parser.add_argument('--rpkm_matrix', required=True, help='')
svd_parser.add_argument('--output', required=False, help='')
svd_parser.set_defaults(func=svd)

#CNV discover
cnv_parser = subparsers.add_parser('discover', help="Run HMM to discover CNVs")
cnv_parser.add_argument('--params', required=True, help='Parameters used by HMM')
cnv_parser.add_argument('--datafile', required=True, help='Read depth file.')
cnv_parser.add_argument('--output', required=True, help='Output file.')
cnv_parser.add_argument('--sample', required=False, default='', help='Optionally, users can choose one sample to run.')
cnv_parser.add_argument('--vcf', required=False, default='', help='Optionally, users can input snp information by specifing a vcf file')
cnv_parser.add_argument('--hetsnp', required=False, default=False)
cnv_parser.add_argument('--tagsnp', required=False, default=False)
cnv_parser.add_argument('--tagsnp_file',required = False, help='TagSNP file location.')
cnv_parser.set_defaults(func=discover)

args = parser.parse_args()
args.func(args)
