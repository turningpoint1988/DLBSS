import numpy as np
import sys, random, os, h5py, argparse, pwd, gzip, string
from os.path import dirname, abspath, join, exists

shapeTypes=("MGW", "Roll", "ProT", "HelT")

complement = string.maketrans('ATCGN', 'TAGCN')

#############################################################################
# Used to take the complement of a string.

def reverseComplement(sequence):
  return sequence.upper().translate(complement)[::-1]


#############################################################################
def getKmer(k, alphabet=None, dir=None):

  kmerfilename="%dmer.seq.gz" % k
  if dir is not None:
    kmerfilename = dir + "/" + kmerfilename

  if not os.path.isfile(kmerfilename):
    kmersDict,kmersList = generateKmer(k, alphabet=alphabet, outputfilename=kmerfilename)
  else:
    kmersDict,kmersList = loadKmer(kmerfilename, k=k)
  return kmersDict,kmersList

####generate all possile k-mers and store in a list

def generateKmer(k, alphabet=None, outputfilename=None):
  assert k>0, "paremeter k=%d needs to be positive" % k
  if alphabet is None:
    alphabet="ACGT"
  alphabet_len=len(alphabet)
#  total=pow(alphabet_len,k)
  kmersList=[""]
  kmersDict={}
  for k_ind in xrange(k):
    # pop out strings of len (k_ind-1)
    poptotal=pow(alphabet_len, k_ind)
    kmersNum=0
    for pop_ind in xrange(poptotal):
      oldseq=kmersList.pop(0)
      for char_ind in xrange(alphabet_len):
        newseq = oldseq+alphabet[char_ind]
        kmersList.append(newseq)
        if k_ind==k-1:
          kmersDict[newseq] = kmersNum
          kmersNum += 1

  # write kmer to gzipped file
  if outputfilename is not None:
    if outputfilename[-3:] == ".gz":
      output = gzip.open(outputfilename, "w")
    else:
      output = open(outputfilename, "w")
    for kmer in kmersList:
      output.write("%s\n" % kmer)
    output.close()

  sys.stderr.write("Generate %d %d-mers\n" % (len(kmersDict),k))
  return kmersDict,kmersList

####load all possile k-mers and store in a list

def loadKmer(inputfilename, k=None):

  if inputfilename[-3:] == ".gz":
    inputfile = gzip.open(inputfilename)
  else:
    inputfile = open(inputfilename)

  kmersList = [line.strip() for line in inputfile]
  kmersDict={}
  for kmer_ind in xrange(len(kmersList)):
    kmersDict[kmersList[kmer_ind]] = kmer_ind

  if k is not None:
    assert k==len(kmersList[0]), "paremeter k=%d needs to be match sequences in inputfile." % k
  else:
    k=len(kmersList[0])

  sys.stderr.write("Read %d %d-mers\n" % (len(kmersDict),k))
  return kmersDict,kmersList
#############################################################################
def getUniqKmer(k, alphabet=None, kmersDict=None, kmersList=None, dir=None):

  uniqKmerfilename="%dmer.uniq.gz" % k
  if dir is not None:
    uniqKmerfilename = dir + "/" + uniqKmerfilename

  if not os.path.isfile(uniqKmerfilename):
    if kmersDict is None or kmersList is None:
      kmersDict,kmersList = getKmer(k, alphabet=alphabet, dir=dir)

    uniqKmersDict,uniqKmersList = generateUniqKmer(k, alphabet=alphabet,
                                                   kmersDict=kmersDict, kmersList=kmersList,
                                                   outputfilename=uniqKmerfilename)
  else:
    uniqKmersDict,uniqKmersList = loadUniqKmer(uniqKmerfilename, k=k)

  return uniqKmersDict,uniqKmersList

####generate all unique kmersList and store in dictionary and List

def generateUniqKmer(k, alphabet=None, kmersDict=None, kmersList=None, outputfilename=None):
  assert k>0, "paremeter k=%d needs to be positive" % k
  if kmersDict is None or kmersList is None:
    kmersDict,kmersList = generateKmer(k=k,alphabet=alphabet)

  uniqKmersDict = {} # key is kmer; values are index of seq and rev-comp-seq
  uniqKmersList = []
  for kmer_ind in xrange(len(kmersList)):
    kmer=kmersList[kmer_ind]
    kmerRevComp=reverseComplement(kmer)
#    sys.stderr.write("%s=>%s\n" % (kmer,kmerRevComp))
    if kmerRevComp in uniqKmersDict:
      uniqKmersDict[kmerRevComp].append(kmer_ind)
    else:
      uniqKmersDict[kmer]=[len(uniqKmersList),kmer_ind]
      uniqKmersList.append(kmer)

  # write kmer to gzipped file
  if outputfilename is not None:
    if outputfilename[-3:] == ".gz":
      output = gzip.open(outputfilename, "w")
    else:
      output = open(outputfilename, "w")
    for uniqKmer in uniqKmersList:
      if len(uniqKmersDict[uniqKmer]) == 2:
        output.write("%s %d %d\n" % (uniqKmer, uniqKmersDict[uniqKmer][0],uniqKmersDict[uniqKmer][1]))
      else:
        output.write("%s %s %d %d %d\n" % (uniqKmer,kmersList[uniqKmersDict[uniqKmer][2]],
                                        uniqKmersDict[uniqKmer][0], uniqKmersDict[uniqKmer][1], uniqKmersDict[uniqKmer][2]))
    output.close()

  sys.stderr.write("Generate %d uniq %d-mers\n" % (len(uniqKmersDict),k))
  return uniqKmersDict,uniqKmersList

####read unique kmer file and re-install uniqKmersDict,uniqKmersList,kmersDict,kmersList

def loadUniqKmer(inputfilename, k=None):
  if inputfilename[-3:] == ".gz":
    inputfile = gzip.open(inputfilename)
  else:
    inputfile = open(inputfilename)

  uniqKmersDict = {} # key is kmer; values are index of seq and rev-comp-seq
  for line in inputfile:
    words = line.rstrip().split()
    if len(words) == 3:
      uniqKmer = words[0]
      uniqKmersDict[uniqKmer] = [int(words[1]),int(words[2])]
    elif len(words) == 5:
      uniqKmer = words[0]
      uniqKmerRev = words[1]
      uniqKmersDict[uniqKmer] = [int(words[2]),int(words[3]),int(words[4])]
    else:
      sys.stderr.write("Error! Wrong file format: %s" % line)

  if k is not None:
    assert k==len(uniqKmersDict.keys()[0]), "paremeter k=%d needs to be match sequences in inputfile." % k
  else:
    k=len(uniqKmer)

  sys.stderr.write("Read %d uniq %d-mers\n" % (len(uniqKmersDict),k))

  uniqKmersList = [None] * len(uniqKmersDict)
  for uniqKmer in uniqKmersDict:
    uniqKmersList[uniqKmersDict[uniqKmer][0]] = uniqKmer

  return uniqKmersDict,uniqKmersList
#############################################################################
# load mapper
def Load_mapper(mapperfile):
    mapper = {}
    with open(mapperfile,'r') as f:
         for x in f:
             line = x.strip().split()
             word = line[0]
             vec = [float(item) for item in line[1:]]
             mapper[word] = vec
    return mapper
#############################################################################
# load 5-mer DNA shape look up table
# key is the 5mer, info is the shape vector
# FWD 5mer and REV 5mer are saved as separate keys
# Roll and HelT have 2 shape values
# if flanking = True, appending up to 2N flanking sequences in the table

def loadShapeTable(shapeFilename, flanking=False):

  shapeDict = {}
  shapeFile = open(shapeFilename)
  for line in shapeFile:
    words = line.rstrip().split("\t")
    assert len(words) == 8, "Wrong lookup table format. %s" % shapeFilename
    fwdKmer = words[0]
    if fwdKmer == "#FWD": # first header line
      continue
    revKmer = words[1]
    MGW = float(words[2])
    ProT = float(words[3])
    Roll1 = float(words[4])
    Roll2 = float(words[5])
    HelT1 = float(words[6])
    HelT2 = float(words[7])
    shapeDict[fwdKmer] = (MGW, ProT, Roll1, Roll2, HelT1, HelT2)
    shapeDict[revKmer] = (MGW, ProT, Roll2, Roll1, HelT2, HelT1)
  shapeFile.close()
  sys.stderr.write("Load shape look up table for all 5mers.\n")

  if flanking is True:

    t4mersDict,t4mersList = generateKmer(4, alphabet=None, outputfilename=None)
    for t4mer in t4mersList:
      # N+(k-1)mer
      newKmer = "N"+t4mer
      kmerShape = np.zeros(6)
      for char in "ACGT":
        kmer = char+t4mer
        kmerShape += shapeDict[kmer]
      kmerShape /= 4.0
      shapeDict[newKmer] = kmerShape
      # (k-1)mer+N
      newKmerRevComp=reverseComplement(newKmer)
      shapeDict[newKmerRevComp] = kmerShape
    nFlankingSeqs = 2*len(t4mersList)
    sys.stderr.write("Added %d N-flanking sequences into the DNA shape look up tabe.\n" % nFlankingSeqs)

    t3mersDict,t3mersList = generateKmer(3, alphabet=None, outputfilename=None)
    for t3mer in t3mersList:
      # NN+(k-2)mer
      newKmer = "NN"+t3mer
      kmerShape = np.zeros(6)
      for char1 in "ACGT":
        for char2 in "ACGT":
          kmer = char1+char2+t3mer
          kmerShape += shapeDict[kmer]
      kmerShape /= 16.0
      shapeDict[newKmer] = kmerShape
      # (k-2)mer+NN
      newKmerRevComp=reverseComplement(newKmer)
      shapeDict[newKmerRevComp] = kmerShape
      # N+(k-2)mer+N
      newKmer = "N"+t3mer+"N"
      kmerShape = np.zeros(6)
      for char1 in "ACGT":
        for char2 in "ACGT":
          kmer = char1+t3mer+char2
          kmerShape += shapeDict[kmer]
      kmerShape /= 16.0
      shapeDict[newKmer] = kmerShape
    nFlankingSeqs = 3*len(t3mersList)
    sys.stderr.write("Added %d NN-flanking sequences into the DNA shape look up tabe.\n" % nFlankingSeqs)

    t2mersDict,t2mersList = generateKmer(2, alphabet=None, outputfilename=None)
    for t2mer in t2mersList:
      # NN+(k-3)+N mer
      newKmer = "NN"+t2mer+"N"
      kmerShape = np.zeros(6)
      for char1 in "ACGT":
        for char2 in "ACGT":
            for char3 in "ACGT":
              kmer = char1+char2+t2mer+char3
              kmerShape += shapeDict[kmer]
      kmerShape /= 64.0
      shapeDict[newKmer] = kmerShape
      # N+(k-2)mer+NN
      newKmerRevComp=reverseComplement(newKmer)
      shapeDict[newKmerRevComp] = kmerShape
    nFlankingSeqs = 2*len(t2mersList)
    sys.stderr.write("Added %d NN+N-flanking sequences into the DNA shape look up tabe.\n" % nFlankingSeqs)

    t1mersDict,t1mersList = generateKmer(1, alphabet=None, outputfilename=None)
    for t1mer in t1mersList:
      # NN+(k-1)mer+NN
      newKmer = "NN"+t1mer+"NN"
      kmerShape = np.zeros(6)
      for char1 in "ACGT":
        for char2 in "ACGT":
          for char3 in "ACGT":
            for char4 in "ACGT":
              kmer = char1+char2+t1mer+char3+char4
              kmerShape += shapeDict[kmer]
      kmerShape /= 256.0
      shapeDict[newKmer] = kmerShape
    nFlankingSeqs = 1*len(t1mersList)
    sys.stderr.write("Added %d NN+NN-flanking sequences into the DNA shape look up tabe.\n" % nFlankingSeqs)

  return shapeDict

#############################################################################
# calculate dnashape values of given sequences using dnashape look up table
# if flanking=True, the sequences can contain up to 2 Ns in the flanking sequences 

def calculateDNAShape(seqs, seqlen, seqNum, shapeDict):
  
    assert seqNum == len(seqs), "seqNum=%d doesn't match input sequence number %d" % (seqNum, len(seqs))
    assert seqlen == len(seqs[0]), "seqlen=%d doesn't match input sequence length %d" % (seqlen, len(seqs[0]))
    seqShapeDict={}
    nRows=seqNum
    for shapeType in shapeTypes:
        if shapeType == "MGW" or shapeType == "ProT":
           nCols=seqlen-4
        else:
           nCols=seqlen-3
        seqShapeDict[shapeType] = np.zeros(shape=(nRows, nCols), dtype=np.float32)

    for index in xrange(seqNum):
        seq= seqs[index]
        RollP = 0.0
        HelTP = 0.0
        for seq_index in xrange(seqlen-4):
            seq5mer = seq[seq_index:seq_index+5]
            MGW, ProT, Roll1, Roll2, HelT1, HelT2 = shapeDict[seq5mer]
            seqShapeDict["MGW"][index,seq_index] = MGW
            seqShapeDict["ProT"][index,seq_index] = ProT
            if seq_index == 0:
               seqShapeDict["Roll"][index,seq_index] = Roll1
               seqShapeDict["HelT"][index,seq_index] = HelT1
            else:
               seqShapeDict["Roll"][index,seq_index] = (Roll1+RollP)/2.0
               seqShapeDict["HelT"][index,seq_index] = (HelT1+HelTP)/2.0
            RollP = Roll2
            HelTP = HelT2
        seqShapeDict["Roll"][index,seqlen-4] = RollP
        seqShapeDict["HelT"][index,seqlen-4] = HelTP

    return seqShapeDict
#############################################################################
# build k-mer kernel. rows are all training seqs, cols are all uniq k-mers

def buildKmerFeature(trainSeqs, k, alphabet=None,
                     uniqKmersDict=None, uniqKmersList=None,
                     kmersDict=None, kmersList=None):

  assert len(trainSeqs) > 0, "Error! Training sequence set can not be empty"
  assert k>0, "paremeter k=%d needs to be positive" % k

  if kmersDict is None or kmersList is None:
    kmersDict,kmersList = generateKmer(k=k,alphabet=alphabet)
  if uniqKmersDict is None or uniqKmersList is None:
    uniqKmersDict,uniqKmersList = generateUniqKmer(k=k,alphabet=alphabet,kmersDict=kmersDict, kmersList=kmersList)

  nRows = len(trainSeqs)
  nCols = len(uniqKmersDict)

  kmerFeatureArray = np.zeros(shape=(nRows, nCols), dtype=np.float32)

  for seq_ind in xrange(nRows):
    for seq_kmer_start in xrange(len(trainSeqs[0])-k+1):
      seq_kmer=trainSeqs[seq_ind][seq_kmer_start:(seq_kmer_start+k)]
      if seq_kmer not in kmersDict:
        sys.stderr.write("Error! %s not a valid %d-mer!\n" % (seq_kmer,k))
        sys.exit(1)
      for j in xrange(nCols):
        jKmerFwd=uniqKmersList[j]
        jKmerInds=uniqKmersDict[jKmerFwd]
        if len(jKmerInds) == 2:
          jKmerRev = jKmerFwd
        else:
          jKmerRev=kmersList[jKmerInds[2]]
        if seq_kmer == jKmerFwd or seq_kmer == jKmerRev:
          kmerFeatureArray[seq_ind,j] += 1

  return kmerFeatureArray
#############################################################################
# store data as hdf5 format
def outputHDF5Shape(data, filename, shapename='shape'):
    print 'data shape: {}'.format(data.shape)
    comp_kwargs = {'compression': 'gzip', 'compression_opts': 1}
    with h5py.File(filename, 'w') as f:
        f.create_dataset(shapename, data=data, **comp_kwargs)
#############################################################################
# store data as hdf5 format
def outputHDF5Data(data, label, filename, labelname='intensity', dataname='data'):
    print 'data shape: {}'.format(data.shape)
    print 'label shape: {}'.format(label.shape)
    comp_kwargs = {'compression': 'gzip', 'compression_opts': 1}
    #label = [[x.astype(np.float32)] for x in label]
    with h5py.File(filename, 'w') as f:
        f.create_dataset(dataname, data=data, **comp_kwargs)
        f.create_dataset(labelname, data=label, **comp_kwargs)
#############################################################################
def reduceDim(array):
    rows, cols = array.shape
    array_d = np.zeros(shape=(rows, cols-1), dtype=np.float32)
    for row in range(rows):
        for col in range(cols-1):
            array_d[row, col] = (array[row, col] + array[row, col+1]) / 2
    return array_d
#############################################################################
def convert_shape(seqs, intensity, seqs_flank, mapper, shapeDict, out_path, name):
    seqs_vector = []    
    for seq in seqs:
       mat = []
       for element in seq:
           if element in mapper.keys(): mat.append(mapper[element])
           else: print >> sys.stderr, 'invalid character!'; sys.exit(1)
       seqs_vector.append(mat)
       
    seqs_vector = np.asarray(seqs_vector)
    intensity = np.asarray(intensity)
    intensity = intensity.reshape((intensity.shape[0], 1))
    # store data and label into out_filename
    out_filename = out_path + '/%s.hdf5' % name
    if os.path.exists(out_filename):
       os.remove(out_filename)
    outputHDF5Data(seqs_vector, intensity, out_filename)
    ## calculate the shape of seqences ##
    seqShapeDict = calculateDNAShape(seqs_flank, len(seqs_flank[0]), len(seqs_flank), shapeDict)
    # "MGW"
    MGW = seqShapeDict["MGW"]
    MGW_arr = MGW.reshape((MGW.shape[0], MGW.shape[1], 1))
    # "ProT"
    ProT = seqShapeDict["ProT"]
    ProT_arr = ProT.reshape((ProT.shape[0], ProT.shape[1], 1))
    # "Roll"
    Roll = seqShapeDict["Roll"]
    Roll_d = reduceDim(Roll)
    Roll_d_arr = Roll_d.reshape((Roll_d.shape[0], Roll_d.shape[1], 1))
    # "HelT"
    HelT = seqShapeDict["HelT"]
    HelT_d = reduceDim(HelT)
    HelT_d_arr = HelT_d.reshape((HelT_d.shape[0], HelT_d.shape[1], 1))

    fourMPRH = np.concatenate((MGW_arr, ProT_arr, Roll_d_arr, HelT_d_arr), axis=-1)
    
    # store fourMPRH in hdf5 format
    out_filename = out_path + '/%s_MPRH.hdf5' % name
    if os.path.exists(out_filename):
       os.remove(out_filename)
    outputHDF5Shape(fourMPRH, out_filename)
    # store MGW in hdf5 format
    out_filename = out_path + '/%s_MGW.hdf5' % name
    if os.path.exists(out_filename):
       os.remove(out_filename)
    outputHDF5Shape(MGW_arr, out_filename)
    # store ProT in hdf5 format
    out_filename = out_path + '/%s_ProT.hdf5' % name
    if os.path.exists(out_filename):
       os.remove(out_filename)
    outputHDF5Shape(ProT_arr, out_filename)
    # store Roll in hdf5 format
    out_filename = out_path + '/%s_Roll.hdf5' % name
    if os.path.exists(out_filename):
       os.remove(out_filename)
    outputHDF5Shape(Roll_d_arr, out_filename)
    # store Roll in hdf5 format
    out_filename = out_path + '/%s_HelT.hdf5' % name
    if os.path.exists(out_filename):
       os.remove(out_filename)
    outputHDF5Shape(HelT_d_arr, out_filename)
#############################################################################
def parse_args():
    parser = argparse.ArgumentParser(description="Convert sequence and target for Caffe")
    user = pwd.getpwuid(os.getuid())[0]
    
    parser.add_argument("-f", dest="file", type=str, default="", help="file containing sequences in FASTA format")
    parser.add_argument("-shape", "--shapefile", dest="shapefile", default="", help="a shape file")

    return parser.parse_args()

if __name__ == "__main__":
    
    args = parse_args()
    file_path = args.file
    mapper = {'A':[1.,0.,0.,0.],'C':[0.,1.,0.,0.],'G':[0.,0.,1.,0.],'T':[0.,0.,0.,1.],'N':[0.,0.,0.,0.]}
    
    if args.shapefile == "":
       print >> sys.stderr, "Need shape file";sys.exit(1)
    else:
       shapeDict = loadShapeTable(args.shapefile, True)

    seqs = []; seqs_flank = []; intensity = []
    with open(file_path, 'r') as f:
        lines = f.readlines()
        
    for line in lines:
        line_split = line.strip().split()
        if len(line_split[1]) != 35:
            print line_split[1]
            continue
        else:
            seqs.append(line_split[1])
            seqs_flank.append('NN' + line_split[1] + 'NN')
            intensity.append(float(line_split[0]))
    
    out_dir = join(dirname(abspath(file_path)), 'data')
    if not exists(out_dir):
        os.mkdir(out_dir)
    convert_shape(seqs, intensity, seqs_flank, mapper, shapeDict, out_dir, 'train')
    
