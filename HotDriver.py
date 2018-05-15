import os,subprocess,sys,getopt,scipy,math
from scipy import stats

def usage():
    print 'Usage:\n'
    print 'Please check the /example_data directory for individual file formats...\n'
    
    print ' -m, --mutation (required)        Mutation data\n'
    print ' -c, --cpgIsland (required)       CpG island region bed file on hg19\n'
    print ' -l, --gene_length (required)     Genes with corresponding amino acid length\n'
    print ' -o, --output (required)          Hotspot mutations that are predicted\n'
    print ' -g, --gene_list (optional)       Genes to investigate. Default is to include all available genes\n'
    print ' -p, --pvalue (optional)          Adjusted pvalue cutoff that is set to identify hotspot mutations. Default: 0.01\n'

class ArgumentError(Exception):
	pass

def poisson_probability(actual, mean):
    # naive:   math.exp(-mean) * mean**actual / factorial(actual)

    # iterative, to keep the components from getting too large or small:
    p = math.exp(-mean)
    for i in xrange(actual):
        p *= mean
        p /= i+1
    return p

def parse_arguments(argv):
    mutation,cpgIsland,gene_length,output = None,None,None,None
    gene_list = ''
    pValCutoff = 0.01
    
    try:
        opts, args = getopt.getopt(argv[1:], "m:g:c:l:o:g:p:h", ["mutation=", "gene_list=", "cpgIsland=", "gene_length=", "output=", "gene_list", "pvalue=", "help"])
    except getopt.GetoptError as err:
        # print help information and exit:
        print str(err)  # will print something like "option -m not recognized"
        usage()

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
            
        elif opt in ("-m", "--mutation"):
            mutation = arg
            if not os.path.isfile(mutation):
                raise ArgumentError("Argument Error: Invalid mutation data passed.")
        
        elif opt in ("-g", "--gene_list"):
            gene_list = arg
            if not os.path.isfile(gene_list):
                raise ArgumentError("Argument Error: Invalid gene list passed.")
        
        
        elif opt in ("-c", "--cpgIsland"):
            cpgIsland = arg
            if not os.path.isfile(cpgIsland):
                raise ArgumentError("Argument Error: Invalid cpgIsland region data passed.")
        
            
        elif opt in ("-l", "--gene_length"):
            gene_length = arg
            if not os.path.isfile(gene_length):
                raise ArgumentError("Argument Error: Invalid gene amino acid length data passed.")
        
        elif opt in ("-p", "--pvalue"):
            pValCutoff = float(arg)
            
        elif opt in ("-o", "--output"):
            output = arg
            
        else:
            raise ArgumentError("Bad argument: I don't know what %s is" % opt)

    if mutation is None or cpgIsland is None or gene_length is None:
        raise ArgumentError("You need to supply mutation file, cpgIsland region, and gene amino acid length information all!") 
    
    # return argument values
    return mutation,gene_list,cpgIsland,gene_length,pValCutoff,output

def run (cmd):
    subprocess.call(cmd, shell = True)
    return

def consistentHigh(aaMutCounts,geneMutCounts,geneLen):
    x=0
    for count in aaMutCounts:
        if int(count)>1:
            x+=1
    for count in geneMutCounts:
        if float(count)/geneLen>=1:
            x+=1
    if x>=1:
        return True
    else:
        return False

def sum(s):
    sum=0
    for x in s:
        sum+=int(x)
    return sum

def is_num(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def pvalue_correction(num_Test,uncorrected_Pval):
    corrected_Pval=dict()
    prev_bh_value = 0
    i=0
    for test_id,pval in sorted(uncorrected_Pval.items(), key=lambda x: x[1]):
        bh_value = pval * num_Test/(i+1)
        bh_value = min(bh_value, 1)
        bh_value = max(bh_value, prev_bh_value)
        prev_bh_value = bh_value
        corrected_Pval[test_id]=bh_value
        i+=1    
    return corrected_Pval

def mutation_attribute(gene,mutType,ref,alt,mutPos,CpGisland):
    if mutType == 'Missense':
        if mutPos not in CpGisland:
            if ref[0]=='A':
                if alt[0]=='G':
                    return 'Missense_ATts'
                else:
                    return 'Missense_ATtv'
            elif ref[0]=='T':
                if alt[0]=='C':
                    return 'Missense_ATts'
                else:
                    return 'Missense_ATtv'
            
            elif ref[0]=='G':
                if alt[0] == 'A':
                    return 'Missense_NoCpG_CGts'
                else:
                    return 'Missense_NoCpG_CGtv'
            elif ref[0]=='C':
                if alt[0] == 'T':
                    return 'Missense_NoCpG_CGts'
                else:
                    return 'Missense_NoCpG_CGtv'

        else:
            if ref[0]=='A':
                if alt[0]=='G':
                    return 'Missense_ATts'
                else:
                    return 'Missense_ATtv'
            elif ref[0]=='T':
                if alt[0]=='C':
                    return 'Missense_ATts'
                else:
                    return 'Missense_ATtv'
            
            elif ref[0]=='G':
                if alt[0] == 'A':
                    return 'Missense_CpG_CGts'
                else:
                    return 'Missense_CpG_CGtv'
            elif ref[0]=='C':
                if alt[0] == 'T':
                    return 'Missense_CpG_CGts'
                else:
                    return 'Missense_CpG_CGtv'
                        
    elif mutType == 'Nonsense':
        if mutPos not in CpGisland:
            if ref[0]=='A':
                if alt[0]=='G':
                    return 'Nonsense_ATts'
                else:
                    return 'Nonsense_ATtv'
            elif ref[0]=='T':
                if alt[0]=='C':
                    return 'Nonsense_ATts'
                else:
                    return 'Nonsense_ATtv'
            
            elif ref[0]=='G':
                if alt[0] == 'A':
                    return 'Nonsense_NoCpG_CGts'
                else:
                    return 'Nonsense_NoCpG_CGtv'
            elif ref[0]=='C':
                if alt[0] == 'T':
                    return 'Nonsense_NoCpG_CGts'
                else:
                    return 'Nonsense_NoCpG_CGtv'

        else:
            if ref[0]=='A':
                if alt[0]=='G':
                    return 'Nonsense_ATts'
                else:
                    return 'Nonsense_ATtv'
            elif ref[0]=='T':
                if alt[0]=='C':
                    return 'Nonsense_ATts'
                else:
                    return 'Nonsense_ATtv'
            
            elif ref[0]=='G':
                if alt[0] == 'A':
                    return 'Nonsense_CpG_CGts'
                else:
                    return 'Nonsense_CpG_CGtv'
            elif ref[0]=='C':
                if alt[0] == 'T':
                    return 'Nonsense_CpG_CGts'
                else:
                    return 'Nonsense_CpG_CGtv'
                
    elif mutType == 'Silent':
        if mutPos not in CpGisland:
            if ref[0]=='A':
                if alt[0]=='G':
                    return 'Silent_ATts'
                else:
                    return 'Silent_ATtv'
            elif ref[0]=='T':
                if alt[0]=='C':
                    return 'Silent_ATts'
                else:
                    return 'Silent_ATtv'
            
            elif ref[0]=='G':
                if alt[0] == 'A':
                    return 'Silent_NoCpG_CGts'
                else:
                    return 'Silent_NoCpG_CGtv'
            elif ref[0]=='C':
                if alt[0] == 'T':
                    return 'Silent_NoCpG_CGts'
                else:
                    return 'Silent_NoCpG_CGtv'

        else:
            if ref[0]=='A':
                if alt[0]=='G':
                    return 'Silent_ATts'
                else:
                    return 'Silent_ATtv'
            elif ref[0]=='T':
                if alt[0]=='C':
                    return 'Silent_ATts'
                else:
                    return 'Silent_ATtv'
            
            elif ref[0]=='G':
                if alt[0] == 'A':
                    return 'Silent_CpG_CGts'
                else:
                    return 'Silent_CpG_CGtv'
            elif ref[0]=='C':
                if alt[0] == 'T':
                    return 'Silent_CpG_CGts'
                else:
                    return 'Silent_CpG_CGtv'
                        
    elif mutType == 'Insertion':
        return 'Insertion'
    elif mutType == 'Deletion':
        return 'Deletion'

mutation,gene_list,cpgIsland,gene_length,pValCutoff,output = parse_arguments(sys.argv)
cwd = os.getcwd()
bedtools = cwd+'/tools/bedtools' ## the user can also define the location of bedtools

print "Check whether the output directory exist..."
if not os.path.isdir('/'.join(output.split('/')[:-1])):
    print "please make sure the output directory exists...\n"
    sys.exit()

## mutation subtypes to include in the hotspot mutation identification
subtypeList = ['Missense_ATts','Missense_ATtv','Missense_CpG_CGts','Missense_CpG_CGtv','Missense_NoCpG_CGts','Missense_NoCpG_CGtv',
               'Nonsense_ATts','Nonsense_ATtv','Nonsense_CpG_CGts','Nonsense_CpG_CGtv','Nonsense_NoCpG_CGts','Nonsense_NoCpG_CGtv',
               'Silent_ATts','Silent_ATtv','Silent_CpG_CGts','Silent_CpG_CGtv','Silent_NoCpG_CGts','Silent_NoCpG_CGtv',
               'Insertion',
               'Deletion']

if gene_list!='':
    print "Record the genes to investigate..."
    InvestigateGenes = set()
    for data in open(gene_list):
        data = data.rstrip().split("\t")
        gene = data[0]
        if gene not in InvestigateGenes:
            InvestigateGenes.add(gene)
else:
    print "Investigate all available genes..."
        
print "Record the CDS amino acid length for selected cancer genes..."
geneLen=dict()
for data in open(gene_length):
    data = data.rstrip().split("\t")
    gene,length = data
    if gene_list!='':
        if gene in InvestigateGenes and length!='NA':
            if gene not in geneLen:
                geneLen[gene]=int(length)
    else:
        if length!='NA':
            if gene not in geneLen:
                geneLen[gene]=int(length)

print "Identify mutations that locate in CpG region..."
mutationBed = output+'.tmpmutation.bed'
intersect =  output+'.CpG.mutation.bed'
cmd = "sed -e '1d' "+mutation+" | cut -f3,4,5 | sort | uniq | sed -e 's/^/chr/' > "+mutationBed+"\n"
cmd += bedtools+' intersect -a '+cpgIsland+' -b '+mutationBed+' > '+intersect+'\n'
cmd += 'rm '+mutationBed+'\n'
run(cmd)

print "Record the CpG island regions..."
CpGisland = dict()
cpgBed = open(intersect)
for data in cpgBed:
    data = data.rstrip().split('\t')
    chrom,start,end=data
    if chrom+':'+start+'-'+end not in CpGisland:
        CpGisland[chrom+':'+start+'-'+end]=""
cmd = 'rm '+intersect+'\n'
run(cmd)

print "Compute the number of each mutation subtypes on individual genes and on individual amino acid positions ..."
geneMuts,geneSubMuts,geneChrMuts=dict(),dict(),dict()
overallGeneMuts=dict()
mutData = open(mutation)
header = mutData.readline()
for data in mutData:
    data = data.rstrip().split("\t")
    sample,gene,chrom,start,end,ref,alt,aaPos,mutType,strand = data
    gene = gene.split('_')[0]
    mutPos = 'chr'+chrom+':'+start+'-'+end
    mutant = chrom+':'+start+'-'+end+ref+'/'+alt
    if gene in geneLen:
        subType=mutation_attribute(gene,mutType,ref,alt,mutPos,CpGisland)
        if subType in subtypeList:
            if gene not in geneMuts:
                geneMuts[gene]=[0]*len(subtypeList)
                geneSubMuts[gene]=dict()
		geneChrMuts[gene]=dict()
        
            loc = subtypeList.index(subType)
            if gene in geneMuts:
                geneMuts[gene][loc]+=1
            if gene in geneSubMuts:
                if aaPos not in geneSubMuts[gene]:
                    geneSubMuts[gene][aaPos]=[0]*len(subtypeList)
                if aaPos in geneSubMuts[gene]:
                    geneSubMuts[gene][aaPos][loc]+=1
	    if gene in geneChrMuts:
		if aaPos not in geneChrMuts[gene]:
		    geneChrMuts[gene][aaPos]=dict()
		if aaPos in geneChrMuts[gene]:
		    if mutant not in geneChrMuts[gene][aaPos]:
			geneChrMuts[gene][aaPos][mutant]=0
		    if mutant in geneChrMuts[gene][aaPos]:
			geneChrMuts[gene][aaPos][mutant]+=1
	
	if mutType == 'Silent':		
	    if gene not in overallGeneMuts:
		overallGeneMuts[gene]=0
	    if gene in overallGeneMuts:
		overallGeneMuts[gene]+=1
                
print "Identify the hotspot mutations ..."
numTest=0
aaTruth,uncorrected_Pval=dict(),dict()
for gene in geneSubMuts:
    geneMutCounts = geneMuts[gene]

    for aaPos in geneSubMuts[gene]:
        aaMutCounts = geneSubMuts[gene][aaPos]
        statistic=0.0
	
	n=0
        gene_aaPos = gene+'_'+aaPos
        for i in range(len(subtypeList)):
            if geneMutCounts[i]!=0:
		n+=1
		#if gene in overallGeneMuts:
		    #expectation=float(overallGeneMuts[gene])/(geneLen[gene]*6)
		#else:
		    #expectation=1.0/(geneLen[gene]*6)
	    #else:
                expectation  = float(geneMutCounts[i])/geneLen[gene]
		observation = aaMutCounts[i]
		p = max(1e-200,poisson_probability(observation,expectation))
		statistic += -2*(math.log(p,10))
        
	uncorrectPval = scipy.stats.chisqprob(statistic,2*n)
        numTest+=1
       
        if gene_aaPos not in uncorrected_Pval:
            uncorrected_Pval[gene_aaPos]=uncorrectPval

        if consistentHigh(aaMutCounts,geneMutCounts,geneLen[gene]):
            aaTruth[gene_aaPos]=""

print "Perform false discovery rate correction ..."            
corrected_Pval = pvalue_correction(numTest,uncorrected_Pval)

print "Output the hotspot mutations ..."
idenHotspot = open(output,'w')
idenHotspot.write('Gene\taaPos\tMutCount\tMutSubtypes\tMutPositions\tPvalue\tadj.Pvalue\n')
for gene_aaPos,adjPval in sorted(corrected_Pval.items(),key=lambda kv:kv[1]):
    
    if len(gene_aaPos.split('_'))==2: 
	gene,aaPos = gene_aaPos.split('_')
	mutCount = sum(geneSubMuts[gene][aaPos])
    
	aaSubMutCounts=list()
	for i in range(len(subtypeList)):
	    if geneSubMuts[gene][aaPos][i]!=0:
		aaSubMutCounts.append(subtypeList[i]+'('+str(geneSubMuts[gene][aaPos][i])+')')
	mutSubtypes = ','.join(aaSubMutCounts)
	pVal = uncorrected_Pval[gene_aaPos]
	
	aaMutInfor=list()
	for pos,count in sorted(geneChrMuts[gene][aaPos].items(),key=lambda kv:kv[1],reverse=True):
	    aaMutInfor.append(pos+'('+str(count)+')')
	mutInfor=','.join(aaMutInfor)
	
	if gene_aaPos in aaTruth and adjPval<pValCutoff:
	    idenHotspot.write('%s\t%s\t%d\t%s\t%s\t%.3e\t%.3e\n' % (gene,aaPos,mutCount,mutSubtypes,mutInfor,pVal,adjPval))

