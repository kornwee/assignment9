import seqbio

def gcContent(seq):
    print('Input ' + seq)
    result = seqbio.SeqCal.gcContent(seq)
    print('GC content = ' + str(result))

def countBases(seq, hideInput = False):
    if not hideInput:
        print('Input ' + seq)

    result = seqbio.SeqCal.countBasesDict(seq)
    print('countBases = ' + str(result))

def countBasesRev(seq):
    print('Input ' + seq)
    seq = seqbio.SeqPattern.reverseSeq(seq)
    seq = seqbio.SeqPattern.reverseComplementSeq(seq)
    countBases(seq, True)

def enzTargetsScan(seq, enz, hideInput = False):
    if not hideInput:
        print('Input ' + seq)

    result = seqbio.SeqPattern.enzTargetsScan(seq, enz)
    print('EcoRI sites = ' + str(result))

def enzTargetsScanRev(seq, enz):
    print('Input ' + seq)
    seq = seqbio.SeqPattern.reverseComplementSeq(seq)
    enzTargetsScan(seq, enz, True)

def transcription(seq, hideInput = False):
    if not hideInput:
        print('Input ' + seq)
    
    result = seqbio.dnaconvert.dna2rna(seq)
    print('Transcription = ' + str(result))

def transcriptionRev(seq):
    print('Input ' + seq)
    seq = seqbio.SeqPattern.reverseComplementSeq(seq)
    transcription(seq, True)

def translation(seq, hideInput = False):
    if not hideInput:
        print('Input ' + seq)
    
    result = seqbio.dnaconvert.dna2protein(seq)
    print('Translation = ' + str(result))

def translationRev(seq):
    print('Input ' + seq)
    seq = seqbio.SeqPattern.reverseComplementSeq(seq)
    translation(seq, True)