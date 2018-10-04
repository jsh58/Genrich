#!/usr/bin/python

# JMG 10/2018

# Produce BED file of N homopolymers for a
#   fasta file (e.g. reference genome).

import sys
import gzip

def openRead(filename):
  '''
  Open filename for reading. '-' indicates stdin.
    '.gz' suffix indicates gzip compression.
  '''
  if filename == '-':
    return sys.stdin
  try:
    if filename[-3:] == '.gz':
      f = gzip.open(filename, 'rb')
    else:
      f = open(filename, 'rU')
  except IOError:
    sys.stderr.write('Error! Cannot open %s for reading\n' % filename)
    sys.exit(-1)
  return f

def openWrite(filename):
  '''
  Open filename for writing. '-' indicates stdout.
    '.gz' suffix indicates gzip compression.
  '''
  if filename == '-':
    return sys.stdout
  try:
    if filename[-3:] == '.gz':
      f = gzip.open(filename, 'wb')
    else:
      f = open(filename, 'w')
  except IOError:
    sys.stderr.write('Error! Cannot open %s for writing\n' % filename)
    sys.exit(-1)
  return f

def printNs(fOut, head, seq, minLen):
  '''
  Print BED intervals of Ns.
  '''
  count = 0
  start = -1
  for i in xrange(len(seq)):
    if seq[i].upper() != 'N':
      if start != -1:
        if i - start >= minLen:
          fOut.write('\t'.join([head, str(start), str(i)]) + '\n')
          count += 1
        start = -1

    else:
      if start == -1:
        start = i

  if start != -1:
    if i - start >= minLen:
      fOut.write('\t'.join([head, str(start), str(i)]) + '\n')
      count += 1

  return count

def parseFasta(fIn, fOut, minLen):
  '''
  Parse fasta file, write output on the fly.
  '''
  count = pureNs = 0
  head = ''    # header (1st space-delim token)
  seq = ''     # sequence is pure Ns
  length = 0   # length of sequence

  # analyze fasta reads
  for line in fIn:
    if line[0] == '>':

      # process previous read
      if head:
        count += 1
        pureNs += printNs(fOut, head, seq, minLen)

      # start new read
      head = line.rstrip().split(' ')[0][1:]
      seq = ''
      length = 0

    elif head:
      # save sequence
      seq += line.rstrip()

  if fIn != sys.stdin:
    fIn.close()

  # process last read
  if head:
    count += 1
    pureNs += printNs(fOut, head, seq, minLen)

  if fOut != sys.stdout:
    fOut.close()

  return count, pureNs

def main():
  '''Main.'''
  args = sys.argv[1:]
  if len(args) < 2:
    sys.stderr.write('Usage: python findNs.py  <input>')
    sys.stderr.write('  <output>  [<minLen>]\n')
    sys.stderr.write('  <input>     Input fasta file\n')
    sys.stderr.write('  <output>    Output BED file of \'N\' homopolymers\n')
    sys.stderr.write('  <minLen>    Minimum length of Ns (def. 100bp)\n')
    sys.exit(-1)

  # get CL args
  fIn = openRead(args[0])
  fOut = openWrite(args[1])
  minLen = 100
  if len(args) > 2:
    minLen = int(args[2])

  # parse fasta
  count, pureNs = parseFasta(fIn, fOut, minLen)

  sys.stderr.write('Total fasta sequences in %s: %d\n' % (args[0], count))
  sys.stderr.write('Intervals of Ns (min. %dbp): %d\n' % (minLen, pureNs))

if __name__ == '__main__':
  main()
