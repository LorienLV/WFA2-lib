/*
 *                             The MIT License
 *
 * Wavefront Alignment Algorithms
 * Copyright (c) 2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 * This file is part of Wavefront Alignment Algorithms.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * PROJECT: Wavefront Alignment Algorithms
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: WaveFront-Alignment module for backtracing alignments
 */

#include "utils/commons.h"
#include "wavefront_backtrace.h"

/*
 * Wavefront type
 */
#define BACKTRACE_TYPE_BITS                   4 // 4-bits for piggyback
#define BACKTRACE_TYPE_MASK 0x000000000000000Fl // Extract mask

#define BACKTRACE_PIGGYBACK_SET(offset,backtrace_type) \
  (( ((int64_t)(offset)) << BACKTRACE_TYPE_BITS) | backtrace_type)

#define BACKTRACE_PIGGYBACK_GET_TYPE(offset) \
  ((offset) & BACKTRACE_TYPE_MASK)
#define BACKTRACE_PIGGYBACK_GET_OFFSET(offset) \
  ((offset) >> BACKTRACE_TYPE_BITS)

typedef enum {
  backtrace_M       = 9,
  backtrace_D2_ext  = 8,
  backtrace_D2_open = 7,
  backtrace_D1_ext  = 6,
  backtrace_D1_open = 5,
  backtrace_I2_ext  = 4,
  backtrace_I2_open = 3,
  backtrace_I1_ext  = 2,
  backtrace_I1_open = 1,
} backtrace_type;

/*
 * Backtrace Trace Patch Match/Mismsmatch
 */
int64_t wavefront_backtrace_misms(
    wavefront_aligner_t* const wf_aligner,
    const int score,
    const int k) {
  if (score < 0) return WAVEFRONT_OFFSET_NULL;
  wavefront_t* const mwavefront = wf_aligner->wf_components.mwavefronts[score];
  if (mwavefront != NULL &&
      mwavefront->lo <= k &&
      k <= mwavefront->hi) {
    return BACKTRACE_PIGGYBACK_SET(mwavefront->offsets[k]+1,backtrace_M);
  } else {
    return WAVEFRONT_OFFSET_NULL;
  }
}
void wavefront_backtrace_matches(
    wavefront_aligner_t* const wf_aligner,
    const int k,
    wf_offset_t offset,
    int num_matches,
    cigar_t* const cigar) {
  // Parameters
  const uint64_t matches_lut = 0x4D4D4D4D4D4D4D4Dul; // Matches LUT = "MMMMMMMM"
  char* operations = cigar->operations + cigar->begin_offset;
  // Update offset first
  cigar->begin_offset -= num_matches;
  // Blocks of 8-matches
  while (num_matches >= 8) {
    operations -= 8;
    *((uint64_t*)(operations+1)) = matches_lut;
    num_matches -= 8;
  }
  // Remaining matches
  int i;
  for (i=0;i<num_matches;++i) {
    *operations = 'M';
    --operations;
  }
}
/*
 * Backtrace Trace Patch Deletion
 */
int64_t wavefront_backtrace_del1_open(
    wavefront_aligner_t* const wf_aligner,
    const int score,
    const int k) {
  if (score < 0) return WAVEFRONT_OFFSET_NULL;
  wavefront_t* const mwavefront = wf_aligner->wf_components.mwavefronts[score];
  if (mwavefront != NULL &&
      mwavefront->lo <= k+1 &&
      k+1 <= mwavefront->hi) {
    return BACKTRACE_PIGGYBACK_SET(mwavefront->offsets[k+1],backtrace_D1_open);
  } else {
    return WAVEFRONT_OFFSET_NULL;
  }
}
int64_t wavefront_backtrace_del2_open(
    wavefront_aligner_t* const wf_aligner,
    const int score,
    const int k) {
  if (score < 0) return WAVEFRONT_OFFSET_NULL;
  wavefront_t* const mwavefront = wf_aligner->wf_components.mwavefronts[score];
  if (mwavefront != NULL &&
      mwavefront->lo <= k+1 &&
      k+1 <= mwavefront->hi) {
    return BACKTRACE_PIGGYBACK_SET(mwavefront->offsets[k+1],backtrace_D2_open);
  } else {
    return WAVEFRONT_OFFSET_NULL;
  }
}
int64_t wavefront_backtrace_del1_ext(
    wavefront_aligner_t* const wf_aligner,
    const int score,
    const int k) {
  if (score < 0) return WAVEFRONT_OFFSET_NULL;
  wavefront_t* const d1wavefront = wf_aligner->wf_components.d1wavefronts[score];
  if (d1wavefront != NULL &&
      d1wavefront->lo <= k+1 &&
      k+1 <= d1wavefront->hi) {
    return BACKTRACE_PIGGYBACK_SET(d1wavefront->offsets[k+1],backtrace_D1_ext);
  } else {
    return WAVEFRONT_OFFSET_NULL;
  }
}
int64_t wavefront_backtrace_del2_ext(
    wavefront_aligner_t* const wf_aligner,
    const int score,
    const int k) {
  if (score < 0) return WAVEFRONT_OFFSET_NULL;
  wavefront_t* const d2wavefront = wf_aligner->wf_components.d2wavefronts[score];
  if (d2wavefront != NULL &&
      d2wavefront->lo <= k+1 &&
      k+1 <= d2wavefront->hi) {
    return BACKTRACE_PIGGYBACK_SET(d2wavefront->offsets[k+1],backtrace_D2_ext);
  } else {
    return WAVEFRONT_OFFSET_NULL;
  }
}
/*
 * Backtrace Trace Patch Insertion
 */
int64_t wavefront_backtrace_ins1_open(
    wavefront_aligner_t* const wf_aligner,
    const int score,
    const int k) {
  if (score < 0) return WAVEFRONT_OFFSET_NULL;
  wavefront_t* const mwavefront = wf_aligner->wf_components.mwavefronts[score];
  if (mwavefront != NULL &&
      mwavefront->lo <= k-1 &&
      k-1 <= mwavefront->hi) {
    return BACKTRACE_PIGGYBACK_SET(mwavefront->offsets[k-1]+1,backtrace_I1_open);
  } else {
    return WAVEFRONT_OFFSET_NULL;
  }
}
int64_t wavefront_backtrace_ins2_open(
    wavefront_aligner_t* const wf_aligner,
    const int score,
    const int k) {
  if (score < 0) return WAVEFRONT_OFFSET_NULL;
  wavefront_t* const mwavefront = wf_aligner->wf_components.mwavefronts[score];
  if (mwavefront != NULL &&
      mwavefront->lo <= k-1 &&
      k-1 <= mwavefront->hi) {
    return BACKTRACE_PIGGYBACK_SET(mwavefront->offsets[k-1]+1,backtrace_I2_open);
  } else {
    return WAVEFRONT_OFFSET_NULL;
  }
}
int64_t wavefront_backtrace_ins1_ext(
    wavefront_aligner_t* const wf_aligner,
    const int score,
    const int k) {
  if (score < 0) return WAVEFRONT_OFFSET_NULL;
  wavefront_t* const i1wavefront = wf_aligner->wf_components.i1wavefronts[score];
  if (i1wavefront != NULL &&
      i1wavefront->lo <= k-1 &&
      k-1 <= i1wavefront->hi) {
    return BACKTRACE_PIGGYBACK_SET(i1wavefront->offsets[k-1]+1,backtrace_I1_ext);
  } else {
    return WAVEFRONT_OFFSET_NULL;
  }
}
int64_t wavefront_backtrace_ins2_ext(
    wavefront_aligner_t* const wf_aligner,
    const int score,
    const int k) {
  if (score < 0) return WAVEFRONT_OFFSET_NULL;
  wavefront_t* const i2wavefront = wf_aligner->wf_components.i2wavefronts[score];
  if (i2wavefront != NULL &&
      i2wavefront->lo <= k-1 &&
      k-1 <= i2wavefront->hi) {
    return BACKTRACE_PIGGYBACK_SET(i2wavefront->offsets[k-1]+1,backtrace_I2_ext);
  } else {
    return WAVEFRONT_OFFSET_NULL;
  }
}
/*
 * Backtrace wavefronts
 */
void wavefront_backtrace_linear(
    wavefront_aligner_t* const wf_aligner,
    const int alignment_score,
    const int alignment_k,
    const wf_offset_t alignment_offset) {
  // Parameters
  wavefront_sequences_t* const sequences = &wf_aligner->sequences;
  const int pattern_length = sequences->pattern_length;
  const int text_length = sequences->text_length;
  const wavefront_penalties_t* const penalties = &wf_aligner->penalties;
  const distance_metric_t distance_metric = penalties->distance_metric;
  // Prepare cigar
  cigar_t* const cigar = wf_aligner->cigar;
  cigar_clear(cigar);
  cigar->end_offset = cigar->max_operations - 1;
  cigar->begin_offset = cigar->max_operations - 2;
  cigar->operations[cigar->end_offset] = '\0';
  // Compute starting location
  int score = alignment_score;
  int k = alignment_k;
  int h = WAVEFRONT_H(alignment_k,alignment_offset);
  int v = WAVEFRONT_V(alignment_k,alignment_offset);
  wf_offset_t offset = alignment_offset;
  // Account for ending insertions/deletions
  if (v < pattern_length) {
    int i = pattern_length - v;
    while (i > 0) {cigar->operations[(cigar->begin_offset)--] = 'D'; --i;};
  }
  if (h < text_length) {
    int i = text_length - h;
    while (i > 0) {cigar->operations[(cigar->begin_offset)--] = 'I'; --i;};
  }
  // Trace the alignment back
  while (v > 0 && h > 0 && score > 0) {
    // Compute scores
    const int mismatch = score - penalties->mismatch;
    const int gap_open1 = score - penalties->gap_opening1;
    // Compute source offsets
    const int64_t misms = (distance_metric != indel) ?
        wavefront_backtrace_misms(wf_aligner,mismatch,k) :
        WAVEFRONT_OFFSET_NULL;
    const int64_t ins = wavefront_backtrace_ins1_open(wf_aligner,gap_open1,k);
    const int64_t del = wavefront_backtrace_del1_open(wf_aligner,gap_open1,k);
    const int64_t max_all = MAX(misms,MAX(ins,del));
    // Check source score
    if (max_all < 0) break; // No source
    // Traceback Matches
    const int max_offset = BACKTRACE_PIGGYBACK_GET_OFFSET(max_all);
    const int num_matches = offset - max_offset;
    wavefront_backtrace_matches(wf_aligner,k,offset,num_matches,cigar);
    offset = max_offset;
    // Update coordinates
    v = WAVEFRONT_V(k,offset);
    h = WAVEFRONT_H(k,offset);
    if (v <= 0 || h <= 0) break;
    // Traceback Operation
    const backtrace_type backtrace_type = BACKTRACE_PIGGYBACK_GET_TYPE(max_all);
    switch (backtrace_type) {
      case backtrace_M:
        score = mismatch;
        cigar->operations[(cigar->begin_offset)--] = 'X';
        --offset;
        break;
      case backtrace_I1_open:
        score = gap_open1;
        cigar->operations[(cigar->begin_offset)--] = 'I';
        --k; --offset;
        break;
      case backtrace_D1_open:
        score = gap_open1;
        cigar->operations[(cigar->begin_offset)--] = 'D';
        ++k;
        break;
      default:
        fprintf(stderr,"[WFA::Backtrace] Wrong type trace.4\n");
        exit(1);
        break;
    }
    // Update coordinates
    v = WAVEFRONT_V(k,offset);
    h = WAVEFRONT_H(k,offset);
  }
  // Account for last operations
  if (v > 0 && h > 0) {
    // Account for beginning series of matches
    const int num_matches = MIN(v,h);
    wavefront_backtrace_matches(wf_aligner,k,offset,num_matches,cigar);
    v -= num_matches;
    h -= num_matches;
  }
  // Account for beginning insertions/deletions
  while (v > 0) {cigar->operations[(cigar->begin_offset)--] = 'D'; --v;};
  while (h > 0) {cigar->operations[(cigar->begin_offset)--] = 'I'; --h;};
  // Set CIGAR
  ++(cigar->begin_offset);
  cigar->score = alignment_score;
}
void wavefront_backtrace_affine(
    wavefront_aligner_t* const wf_aligner,
    const affine2p_matrix_type component_begin,
    const affine2p_matrix_type component_end,
    const int alignment_score,
    const int alignment_k,
    const wf_offset_t alignment_offset) {
  // Parameters
  wavefront_sequences_t* const sequences = &wf_aligner->sequences;
  const int pattern_length = sequences->pattern_length;
  const int text_length = sequences->text_length;
  const wavefront_penalties_t* const penalties = &wf_aligner->penalties;
  const distance_metric_t distance_metric = penalties->distance_metric;
  // Prepare cigar
  cigar_t* const cigar = wf_aligner->cigar;
  cigar_clear(cigar);
  cigar->end_offset = cigar->max_operations - 1;
  cigar->begin_offset = cigar->max_operations - 2;
  cigar->operations[cigar->end_offset] = '\0';
  // Compute starting location
  affine2p_matrix_type matrix_type = component_end;
  int score = alignment_score;
  int k = alignment_k;
  int h = WAVEFRONT_H(alignment_k,alignment_offset);
  int v = WAVEFRONT_V(alignment_k,alignment_offset);
  wf_offset_t offset = alignment_offset;
  // Account for ending insertions/deletions
  if (component_end == affine2p_matrix_M) { // ends-free
    if (v < pattern_length) {
      int i = pattern_length - v;
      while (i > 0) {cigar->operations[(cigar->begin_offset)--] = 'D'; --i;};
    }
    if (h < text_length) {
      int i = text_length - h;
      while (i > 0) {cigar->operations[(cigar->begin_offset)--] = 'I'; --i;};
    }
  }
  // Trace the alignment back
  while (v > 0 && h > 0 && score > 0) {
    // Compute scores
    const int mismatch = score - penalties->mismatch;
    const int gap_open1 = score - penalties->gap_opening1 - penalties->gap_extension1;
    const int gap_open2 = score - penalties->gap_opening2 - penalties->gap_extension2;
    const int gap_extend1 = score - penalties->gap_extension1;
    const int gap_extend2 = score - penalties->gap_extension2;
    // Compute source offsets
    int64_t max_all;
    switch (matrix_type) {
      case affine2p_matrix_M: {
        const int64_t misms = wavefront_backtrace_misms(wf_aligner,mismatch,k);
        const int64_t ins1_open = wavefront_backtrace_ins1_open(wf_aligner,gap_open1,k);
        const int64_t ins1_ext  = wavefront_backtrace_ins1_ext(wf_aligner,gap_extend1,k);
        const int64_t max_ins1 = MAX(ins1_open,ins1_ext);
        const int64_t del1_open = wavefront_backtrace_del1_open(wf_aligner,gap_open1,k);
        const int64_t del1_ext  = wavefront_backtrace_del1_ext(wf_aligner,gap_extend1,k);
        const int64_t max_del1 = MAX(del1_open,del1_ext);
        if (distance_metric == gap_affine) {
          max_all = MAX(misms,MAX(max_ins1,max_del1));
          break;
        }
        const int64_t ins2_open = wavefront_backtrace_ins2_open(wf_aligner,gap_open2,k);
        const int64_t ins2_ext  = wavefront_backtrace_ins2_ext(wf_aligner,gap_extend2,k);
        const int64_t max_ins2 = MAX(ins2_open,ins2_ext);
        const int64_t del2_open = wavefront_backtrace_del2_open(wf_aligner,gap_open2,k);
        const int64_t del2_ext  = wavefront_backtrace_del2_ext(wf_aligner,gap_extend2,k);
        const int64_t max_del2 = MAX(del2_open,del2_ext);
        const int64_t max_ins = MAX(max_ins1,max_ins2);
        const int64_t max_del = MAX(max_del1,max_del2);
        max_all = MAX(misms,MAX(max_ins,max_del));
        break;
      }
      case affine2p_matrix_I1: {
        const int64_t ins1_open = wavefront_backtrace_ins1_open(wf_aligner,gap_open1,k);
        const int64_t ins1_ext  = wavefront_backtrace_ins1_ext(wf_aligner,gap_extend1,k);
        max_all = MAX(ins1_open,ins1_ext);
        break;
      }
      case affine2p_matrix_I2: {
        const int64_t ins2_open = wavefront_backtrace_ins2_open(wf_aligner,gap_open2,k);
        const int64_t ins2_ext  = wavefront_backtrace_ins2_ext(wf_aligner,gap_extend2,k);
        max_all = MAX(ins2_open,ins2_ext);
        break;
      }
      case affine2p_matrix_D1: {
        const int64_t del1_open = wavefront_backtrace_del1_open(wf_aligner,gap_open1,k);
        const int64_t del1_ext  = wavefront_backtrace_del1_ext(wf_aligner,gap_extend1,k);
        max_all = MAX(del1_open,del1_ext);
        break;
      }
      case affine2p_matrix_D2: {
        const int64_t del2_open = wavefront_backtrace_del2_open(wf_aligner,gap_open2,k);
        const int64_t del2_ext  = wavefront_backtrace_del2_ext(wf_aligner,gap_extend2,k);
        max_all = MAX(del2_open,del2_ext);
        break;
      }
      default:
        fprintf(stderr,"[WFA::Backtrace] Wrong type trace.1\n");
        exit(1);
        break;
    }
    // Check source score
    if (max_all < 0) break; // No source
    // Traceback matches
    if (matrix_type == affine2p_matrix_M) {
      const int max_offset = BACKTRACE_PIGGYBACK_GET_OFFSET(max_all);
      const int num_matches = offset - max_offset;
      wavefront_backtrace_matches(wf_aligner,k,offset,num_matches,cigar);
      offset = max_offset;
      // Update coordinates
      v = WAVEFRONT_V(k,offset);
      h = WAVEFRONT_H(k,offset);
      if (v <= 0 || h <= 0) break;
    }
    // Traceback operation
    const backtrace_type backtrace_type = BACKTRACE_PIGGYBACK_GET_TYPE(max_all);
    switch (backtrace_type) {
      case backtrace_M:
        score = mismatch;
        matrix_type = affine2p_matrix_M;
        break;
      case backtrace_I1_open:
        score = gap_open1;
        matrix_type = affine2p_matrix_M;
        break;
      case backtrace_I1_ext:
        score = gap_extend1;
        matrix_type = affine2p_matrix_I1;
        break;
      case backtrace_I2_open:
        score = gap_open2;
        matrix_type = affine2p_matrix_M;
        break;
      case backtrace_I2_ext:
        score = gap_extend2;
        matrix_type = affine2p_matrix_I2;
        break;
      case backtrace_D1_open:
        score = gap_open1;
        matrix_type = affine2p_matrix_M;
        break;
      case backtrace_D1_ext:
        score = gap_extend1;
        matrix_type = affine2p_matrix_D1;
        break;
      case backtrace_D2_open:
        score = gap_open2;
        matrix_type = affine2p_matrix_M;
        break;
      case backtrace_D2_ext:
        score = gap_extend2;
        matrix_type = affine2p_matrix_D2;
        break;
      default:
        fprintf(stderr,"[WFA::Backtrace] Wrong type trace.2\n");
        exit(1);
        break;
    }
    switch (backtrace_type) {
      case backtrace_M:
        cigar->operations[(cigar->begin_offset)--] = 'X';
        --offset;
        break;
      case backtrace_I1_open:
      case backtrace_I1_ext:
      case backtrace_I2_open:
      case backtrace_I2_ext:
        cigar->operations[(cigar->begin_offset)--] = 'I';
        --k; --offset;
        break;
      case backtrace_D1_open:
      case backtrace_D1_ext:
      case backtrace_D2_open:
      case backtrace_D2_ext:
        cigar->operations[(cigar->begin_offset)--] = 'D';
        ++k;
        break;
      default:
        fprintf(stderr,"[WFA::Backtrace] Wrong type trace.3\n");
        exit(1);
        break;
    }
    // Update coordinates
    v = WAVEFRONT_V(k,offset);
    h = WAVEFRONT_H(k,offset);
  }
  // Account for last operations
  if (matrix_type == affine2p_matrix_M) {
    if (v > 0 && h > 0) {
      // Account for beginning series of matches
      const int num_matches = MIN(v,h);
      wavefront_backtrace_matches(wf_aligner,k,offset,num_matches,cigar);
      v -= num_matches;
      h -= num_matches;
    }
    // Account for beginning insertions/deletions
    while (v > 0) {cigar->operations[(cigar->begin_offset)--] = 'D'; --v;};
    while (h > 0) {cigar->operations[(cigar->begin_offset)--] = 'I'; --h;};
  } else {
    // DEBUG
    if (v != 0 || h != 0 || (score != 0 && penalties->match == 0)) {
      fprintf(stderr,"[WFA::Backtrace] I?/D?-Beginning backtrace error\n");
      fprintf(stderr,">%.*s\n",pattern_length,sequences->pattern);
      fprintf(stderr,"<%.*s\n",text_length,sequences->text);
      exit(-1);
    }
  }
  // Set CIGAR
  ++(cigar->begin_offset);
  cigar->score = alignment_score;
}

/**
 * Tests whether mwavefront[k - 1] contains an offset coherent with the chain of
 * insertions that we are following. In case it does, we have found a path
 * back to M. If a path back to M is found, updates the cigar with the
 * insertions between the current positions and the original position in the M
 * matrix. Also updates v and h with the current position in the M matrix.
 *
 *
 * @param cigar The cigar.
 * @param mwavefront A wavefront in the M matrix.
 * @param k The diagonal of the current cell. We are going to access
 * mwavefront[k - 1].
 * @param v_lo The lowest allowed v coordinate so the path is considered valid
 * (inclusive).
 * @param v_hi The highest allowed v coordinate so the path is considered valid
 * (inclusive).
 * @param v The original v coordinate when we started following the chain of
 * insertions.
 * @param h The original h coordinate when we started following the chain of
 * insertions.
 * @return True if a path back to M is found, false otherwise.
 */
static bool backtrace_ins_m_only(cigar_t* const cigar,
                                 const wavefront_t* const mwavefront,
                                 const int k,
                                 const int v_lo,
                                 const int v_hi,
                                 int* const v,
                                 int* const h) {

  if (mwavefront == NULL || mwavefront->lo > k - 1 || k - 1 > mwavefront->hi) {
    return false;
  }

  const int new_offset = mwavefront->offsets[k - 1];
  const int new_v = WAVEFRONT_V(k - 1, new_offset);
  const int new_h = WAVEFRONT_H(k - 1, new_offset);

  if (new_v < v_lo || new_v > v_hi) {
    return false;
  }

  // We have found a path back to M.

  // Remove the unwanted matches previously added to the cigar.
  const int pos_in_range = (new_v - v_lo);
  cigar->begin_offset += pos_in_range;

  const int nins = (*h + pos_in_range)- new_h;
  for (int i = 0; i < nins; ++i) {
    cigar->operations[(cigar->begin_offset)--] = 'I';
  }

  *h = new_h;
  *v = new_v;

  return true;
}

/**
 * Tests whether mwavefront[k + 1] contains an offset coherent with the chain of
 * deletions that we are following. In case it does, we have found a path
 * back to M. If a path back to M is found, updates the cigar with the
 * deletions between the current positions and the original position in the M
 * matrix. Also updates v and h with the current position in the M matrix.
 *
 *
 * @param cigar The cigar.
 * @param mwavefront A wavefront in the M matrix.
 * @param k The diagonal of the current cell. We are going to access
 * mwavefront[k + 1].
 * @param h_lo The lowest allowed h coordinate so the path is considered valid
 * (inclusive).
 * @param h_hi The highest allowed h coordinate so the path is considered valid
 * (inclusive).
 * @param v The original v coordinate when we started following the chain of
 * insertions.
 * @param h The original h coordinate when we started following the chain of
 * insertions.
 * @return True if a path back to M is found, false otherwise.
 */
static bool backtrace_del_m_only(cigar_t* const cigar,
                                 const wavefront_t* const mwavefront,
                                 const int k,
                                 const int h_lo,
                                 const int h_hi,
                                 int* const v,
                                 int* const h) {

  if (mwavefront == NULL || mwavefront->lo > k + 1 || k + 1 > mwavefront->hi) {
    return false;
  }

  const int new_offset = mwavefront->offsets[k + 1];
  const int new_v = WAVEFRONT_V(k + 1, new_offset);
  const int new_h = WAVEFRONT_H(k + 1, new_offset);

  if (new_h < h_lo || new_h > h_hi) {
    return false;
  }

  // We have found a path back to M.

  // Remove the unwanted matches previously added to the cigar.
  const int pos_in_range = (new_h - h_lo);
  cigar->begin_offset += pos_in_range;

  const int ndel = (*v + pos_in_range)- new_v;
  for (int i = 0; i < ndel; ++i) {
    cigar->operations[(cigar->begin_offset)--] = 'D';
  }

  *h = new_h;
  *v = new_v;

  return true;
}

#if 1
/**
 * Check that the cigar produced by the M-only backtrace is coherent, i.e.,
 * it has the expected score and matches and mismatches are consistent with
 * the pattern and text.
 * 
 * @param penalties The penalties.
 * @param cigar The cigar produced by the M-only backtrace.
 * @param sequences The sequences.
 * @param expected_score The expected score of the alignment.
 * @param affine2p True if the alignment is affine2p, false otherwise.
 * @return True if the cigar is coherent, false otherwise.
 */
static bool
check_cigar_backtrace_affine_m_only(const wavefront_penalties_t* const penalties,
                                    const cigar_t* const cigar,
                                    const wavefront_sequences_t* const sequences,
                                    const int expected_score,
                                    const bool affine2p) {

    const char* const pattern = sequences->pattern;
    const char* const text = sequences->text;

    int score = 0;
    int score2 = 0;   // For Dual affine.

    int v = 0;
    int h = 0;

    char prev_op = ' ';

    for (int i = cigar->begin_offset; i < cigar->end_offset; ++i) {
      const char op = cigar->operations[i];

      if ((prev_op == 'I' || prev_op == 'D') && op != prev_op && affine2p) {
        // We have finished a chain of gaps, get the best score.
        score = MIN(score, score2);
      }

      if (op == 'M') {
        if (pattern[v] != text[h]) {
          return false;
        }

        score += penalties->match;

        ++v;
        ++h;
      }
      else if (op == 'X') {
        if (pattern[v] == text[h]) {
          return false;
        }

        score += penalties->mismatch;

        ++v;
        ++h;
      }
      else if (op == 'I') {
        if (prev_op != 'I') {
          score2 = score + penalties->gap_opening2;
          score += penalties->gap_opening1;
        }

        score += penalties->gap_extension1;
        score2 += penalties->gap_extension2;

        ++h;
      }
      else if (op == 'D') {
        if (prev_op != 'D') {
          score2 = score + penalties->gap_opening2;
          score += penalties->gap_opening1;
        }

        score += penalties->gap_extension1;
        score2 += penalties->gap_extension2;

        ++v;
      }
      else {
        return false;
      }

      prev_op = op;
    }

    // In case the sequence ends with a gap.
    if ((prev_op == 'I' || prev_op == 'D') && affine2p) {
      score = MIN(score, score2);
    }

    if (score != expected_score) {
      return false;
    }

    return true;
}

#endif

/**
 * Retrieve the cigar of the alignment for gap-affine and dual gap-affine
 * using exclusively the information in the M matrix (M wavefronts). This is
 * slightly slower than the general backtracking, but it enables storing only
 * the M matrix when performing the alignment (I1, I2, D1 and D2 only need a
 * small scope).
 *
 * @param wf_aligner The wavefront aligner.
 * @param component_begin The matrix where the alignment starts. Unused.
 * @param component_end The matrix where the alignment ends. We assume is the M
 * matrix. Unused.
 * @param alignment_score The score of the alignment.
 * @param alignment_k The diagonal that contains the cell (N, M), where the
 * backtracking starts.
 * @param alignment_offset The offset of the cell (N, M), where the backtracking
 * starts.
 */
void wavefront_backtrace_affine_m_only(
    wavefront_aligner_t* const wf_aligner,
    const affine2p_matrix_type component_begin,
    const affine2p_matrix_type component_end,
    const int alignment_score,
    const int alignment_k,
    const wf_offset_t alignment_offset) {

  // Parameters
  wavefront_sequences_t* const sequences = &wf_aligner->sequences;
  const int pattern_length = sequences->pattern_length;
  const int text_length = sequences->text_length;
  const wavefront_penalties_t* const penalties = &wf_aligner->penalties;
  const distance_metric_t distance_metric = penalties->distance_metric;
  // Prepare cigar
  cigar_t* const cigar = wf_aligner->cigar;
  cigar_clear(cigar);
  cigar->end_offset = cigar->max_operations - 1;
  cigar->begin_offset = cigar->max_operations - 2;
  cigar->operations[cigar->end_offset] = '\0';

  // TODO: Be sure that this is correct.
  bool in_mmatrix = true; // In this function, we always start in the M matrix.
  assert(component_end == affine2p_matrix_M);

  int score = alignment_score;
  int h = WAVEFRONT_H(alignment_k,alignment_offset);
  int v = WAVEFRONT_V(alignment_k,alignment_offset);

  // Variables for searching paths back to M from I1, D1, I2 and D2.
  int k_ins = -1; // I1, I2.
  int k_del = -1; // D1, D2.

  int score1 = -1; // I1, D1.
  int score2 = -1; // I2, D2.

  // I1, I2.
  int v_lo = -1;
  int v_hi = -1;
  // D1, D2.
  int h_lo = -1;
  int h_hi = -1;

  // Account for ending insertions/deletions
  if (component_end == affine2p_matrix_M) { // ends-free
    if (v < pattern_length) {
      int i = pattern_length - v;
      while (i > 0) {cigar->operations[(cigar->begin_offset)--] = 'D'; --i;};
    }
    if (h < text_length) {
      int i = text_length - h;
      while (i > 0) {cigar->operations[(cigar->begin_offset)--] = 'I'; --i;};
    }
  }

  // Trace the alignment back
  while (v > 0 || h > 0) {
    if (in_mmatrix) {
      // Extend (find matches in the diagonal) backwards.
      int init_v = v;
      int init_h = h;

      const char* const pattern = wf_aligner->sequences.pattern;
      const char* const text = wf_aligner->sequences.text;

      while (v > 0 && h > 0 && pattern[v - 1] == text[h - 1]) {
        --v;
        --h;
        // If we come from ins or del, then this may be overwritten.
        cigar->operations[(cigar->begin_offset)--] = 'M';
      }

      // Here h and v might be zero, but we do not care, the loop will be exited
      // in the next iteration.

      const int mismatch = score - penalties->mismatch;
      const wavefront_t* const mwavefront = (mismatch >= 0) ?
        wf_aligner->wf_components.mwavefronts[mismatch]
        : NULL;

      const int k = DPMATRIX_DIAGONAL(h,v);
      const int offset = DPMATRIX_OFFSET(h,v);

      //  If we come from a mismatch, then the backwards extend is correct.
      if (mwavefront != NULL &&
          mwavefront->lo <= k &&
          k <= mwavefront->hi &&
          mwavefront->offsets[k] + 1 == offset) {

        --v;
        --h;
        score = mismatch;

        cigar->operations[(cigar->begin_offset)--] = 'X';
      } else {
          // Otherwise, we come from either I1, D2, I2 or D2.
          // Freeze v and h and start searching a path back to M.

          // In WFA, for a given diagonal and score, we only store the furthest
          // reaching offset. We do not know which was the offset prior to
          // extending it. To account for that, we must allow the indel to come
          // at any point between the offset previos performing the backwards
          // extension and the current offset. We store the range of allowed v
          // and h coordinates.
          //
          // Example, in the following DP table, where there is a chain of 3
          // matches an insertion can come from any of the positions marked with
          // '>'.
          //
          //      A  A  A
          //  A > M
          //  A    > M
          //  A       > M
          //
          in_mmatrix = false;

          h_lo = h;
          h_hi = init_h;

          v_lo = v;
          v_hi = init_v;

          k_ins = k;
          k_del = k;

          score1 = score;
          score2 = score;
      }
    } else {
      // We are searching a path back to M in I1, D1, I2 and D2.
      // Any path that leads to M is valid.

      // --- Follow the ins1 path ---
      {
        const int ins1 = score1 - penalties->gap_opening1 - penalties->gap_extension1;
        const wavefront_t* const mwavefront = (ins1 >= 0) ?
          wf_aligner->wf_components.mwavefronts[ins1]
          : NULL;

        const bool path_found = backtrace_ins_m_only(cigar, mwavefront, k_ins,
                                                     v_lo, v_hi, &v, &h);

        if (path_found) {
          score = ins1;
          in_mmatrix = true;
          continue;
        }
      }

      // --- Follow the del1 path ---
      {
        const int del1 = score1 - penalties->gap_opening1 - penalties->gap_extension1;
        const wavefront_t* const mwavefront = (del1 >= 0) ?
          wf_aligner->wf_components.mwavefronts[del1]
          : NULL;

        const bool path_found = backtrace_del_m_only(cigar, mwavefront, k_del,
                                                     h_lo, h_hi, &v, &h);
        if (path_found) {
          score = del1;
          in_mmatrix = true;
          continue;
        }
      }

      if (distance_metric == gap_affine_2p) {
        // --- Follow the ins2 path ---
        {
          const int ins2 = score2 - penalties->gap_opening2 - penalties->gap_extension2;
          const wavefront_t* const mwavefront = (ins2 >= 0) ?
            wf_aligner->wf_components.mwavefronts[ins2]
            : NULL;

          const bool path_found = backtrace_ins_m_only(cigar, mwavefront, k_ins,
                                                       v_lo, v_hi, &v, &h);

          if (path_found) {
            score = ins2;
            in_mmatrix = true;
            continue;
          }
        }

        // --- Follow the del2 path ---
        {
          const int del2 = score2 - penalties->gap_opening2 - penalties->gap_extension2;
          const wavefront_t* const mwavefront = (del2 >= 0) ?
            wf_aligner->wf_components.mwavefronts[del2]
            : NULL;

          const bool path_found = backtrace_del_m_only(cigar, mwavefront, k_del,
                                                       h_lo, h_hi, &v, &h);

          if (path_found) {
            score = del2;
            in_mmatrix = true;
            continue;
          }
        }
      }

      // Keep exploring I1, I2, D1 and D2.
      --k_ins;
      ++k_del;
      score1 -= penalties->gap_extension1;
      score2 -= penalties->gap_extension2;
    }
  }

  // DEBUG
  if (v != 0 || h != 0 || (score != 0 && penalties->match == 0)) {
    fprintf(stderr,"[WFA::Backtrace] I?/D?-Beginning backtrace error\n");
    fprintf(stderr,">%.*s\n",pattern_length,sequences->pattern);
    fprintf(stderr,"<%.*s\n",text_length,sequences->text);
    exit(-1);
  }

  // Set CIGAR
  ++(cigar->begin_offset);
  cigar->score = alignment_score;

#if 1
  const bool ok = check_cigar_backtrace_affine_m_only(penalties,
                                                      cigar,
                                                      sequences,
                                                      alignment_score,
                                                      distance_metric == gap_affine_2p);
  if (!ok) {
    fprintf(stderr, "[WFA::Backtrace] M-only backtrace not coherent\n");
    exit(-1);
  }
#endif
}
/*
 * Backtrace from BT-Buffer (pcigar)
 */
void wavefront_backtrace_pcigar(
    wavefront_aligner_t* const wf_aligner,
    const int alignment_k,
    const int alignment_offset,
    const pcigar_t pcigar_last,
    const bt_block_idx_t prev_idx_last) {
  // Parameters
  wf_backtrace_buffer_t* const bt_buffer =  wf_aligner->wf_components.bt_buffer;
  // Traceback pcigar-blocks
  bt_block_t bt_block_last = {
      .pcigar = pcigar_last,
      .prev_idx = prev_idx_last
  };
  bt_block_t* const init_block = wf_backtrace_buffer_traceback_pcigar(bt_buffer,&bt_block_last);
  // Fetch initial coordinate
  const int init_position_offset = init_block->pcigar;
  wf_backtrace_init_pos_t* const backtrace_init_pos =
      vector_get_elm(bt_buffer->alignment_init_pos,init_position_offset,wf_backtrace_init_pos_t);
  // Unpack pcigar blocks (packed alignment)
  const int begin_v = backtrace_init_pos->v;
  const int begin_h = backtrace_init_pos->h;
  const int end_v = WAVEFRONT_V(alignment_k,alignment_offset);
  const int end_h = WAVEFRONT_H(alignment_k,alignment_offset);
  if (wf_aligner->penalties.distance_metric <= gap_linear) {
    wf_backtrace_buffer_unpack_cigar_linear(
        bt_buffer,&wf_aligner->sequences,
        begin_v,begin_h,end_v,end_h,wf_aligner->cigar);
  } else {
    wf_backtrace_buffer_unpack_cigar_affine(
        bt_buffer,&wf_aligner->sequences,
        begin_v,begin_h,end_v,end_h,wf_aligner->cigar);
  }
}
