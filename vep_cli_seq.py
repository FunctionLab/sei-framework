#Usage: python sequence_class_seqs.py seq_file1 seq_file2
#This script provides the basic capability to predict sequence class 
#variant effect scores from a pair of two sequences (corresponding to REF and ALT).
#Each input file should only the sequence without line break or header.
#Sequence class diff scores are printed in the order of sequence class
#size as in the README doc.

import sys
import torch
import numpy as np
from torch import nn
from selene_sdk.sequences import Genome
from selene_sdk.utils import NonStrandSpecific
from model.sei import Sei


if __name__ == "__main__":
    sei_model = nn.DataParallel(NonStrandSpecific(Sei()))
    sei_model.load_state_dict(torch.load('./model/sei.pth'))

    clustervfeat = np.load('./model/projvec_targets.npy')
    histone_inds = np.load('./model/histone_inds.npy')

    refseqstr = str(np.loadtxt(sys.argv[1], dtype='str'))
    altseqstr = str(np.loadtxt(sys.argv[2], dtype='str'))

    sei_model.cuda()
    chromatin_profile_ref = sei_model(torch.FloatTensor(Genome.sequence_to_encoding(refseqstr)[None,:,:]).transpose(1,2).cuda()).detach().cpu().numpy()
    chromatin_profile_alt = sei_model(torch.FloatTensor(Genome.sequence_to_encoding(altseqstr)[None,:,:]).transpose(1,2).cuda()).detach().cpu().numpy()

    chromatin_profile_ref_adjust = chromatin_profile_ref.copy()
    chromatin_profile_ref_adjust[:, histone_inds] = \
        chromatin_profile_ref_adjust[:, histone_inds] * (
        (np.sum(chromatin_profile_ref[:, histone_inds], axis=1)*0.5 +
         np.sum(chromatin_profile_alt[:, histone_inds], axis=1)*0.5) /
        np.sum(chromatin_profile_ref[:, histone_inds], axis=1))[:, None]

    chromatin_profile_alt_adjust = chromatin_profile_alt.copy()
    chromatin_profile_alt_adjust[:, histone_inds] = \
        chromatin_profile_alt_adjust[:, histone_inds] * (
        (np.sum(chromatin_profile_ref[:, histone_inds], axis=1)*0.5 +
         np.sum(chromatin_profile_alt[:, histone_inds], axis=1)*0.5) /
        np.sum(chromatin_profile_alt[:, histone_inds], axis=1))[:, None]

    refproj = (np.dot(chromatin_profile_ref_adjust,clustervfeat.T) /
               np.linalg.norm(clustervfeat,axis=1))
    altproj = (np.dot(chromatin_profile_alt_adjust,clustervfeat.T) /
               np.linalg.norm(clustervfeat,axis=1))
    diffproj = altproj[:,:40] - refproj[:,:40]

    print('\t'.join(diffproj[0].astype(str).tolist()))