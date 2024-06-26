// wiggle.cpp for CSEM

#include <cstring>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <algorithm>

#include <stdint.h>
#include "sam/bam.h"
#include "sam/sam.h"

#include "utils.h"
#include "wiggle.h"

bool no_fractional_weight = false;
int fragment_length = -1;
bool only_midpoint = false;

void add_bam_record_to_wiggle(const bam1_t *b, Wiggle& wiggle) {
    float w;

    if (no_fractional_weight) w = 1.0;
    else {
      uint8_t *p_tag = bam_aux_get(b, "ZW");
      if (p_tag == NULL) return;
      w = bam_aux2f(p_tag);
    }

    // sanity check
    int totlen = 0;
    uint32_t *p = bam1_cigar(b);
    for (int i = 0; i < (int)b->core.n_cigar; i++, ++p) {
      int op = *p & BAM_CIGAR_MASK;
      int op_len = *p >> BAM_CIGAR_SHIFT;

      switch (op) {
      case BAM_CMATCH : totlen += op_len; break;
      case BAM_CEQUAL : totlen += op_len; break;
      case BAM_CDIFF  : totlen += op_len; break;
      default : assert(false);
      }
    }
    assert(totlen == b->core.l_qseq);

    if (fragment_length < 0) {
      for (int i = 0; i < b->core.l_qseq; i++) wiggle.read_depth[b->core.pos + i] += w;
    }
    else {
      int start, end;

      start = end = -1;

      if ((b->core.flag & 0x0001) == 0) {
	start = std::max(0, ((b->core.flag & 0x0010) ? b->core.pos + b->core.l_qseq - fragment_length : b->core.pos));
	end = std::min((int)wiggle.length, ((b->core.flag & 0x0010) ? b->core.pos + b->core.l_qseq : b->core.pos + fragment_length));
      }
      else if (b->core.isize > 0) { start = b->core.pos; end = start + b->core.isize; }

      if (start < end && only_midpoint) { start = (start + end - 1) / 2; end = start + 1; }

      for (int i = start; i < end; i++) wiggle.read_depth[i] += w;
    }
}

void build_wiggles(const std::string& bam_filename,
                   WiggleProcessor& processor) {

	samfile_t *bam_in = samopen(bam_filename.c_str(), "rb", NULL);
	if (bam_in == 0) { fprintf(stderr, "Cannot open %s!\n", bam_filename.c_str()); exit(-1); }

	bam_header_t *header = bam_in->header;
	bool *used = new bool[header->n_targets];
	memset(used, 0, sizeof(bool) * header->n_targets);

	int cur_tid = -1; //current tid;
	HIT_INT_TYPE cnt = 0;
	bam1_t *b = bam_init1();
	Wiggle wiggle;
	while (samread(bam_in, b) >= 0) {
		if (b->core.flag & 0x0004) continue;

		if (b->core.tid != cur_tid) {
			if (cur_tid >= 0) { used[cur_tid] = true; processor.process(wiggle); }
			cur_tid = b->core.tid;
			wiggle.name = header->target_name[cur_tid];
			wiggle.length = header->target_len[cur_tid];
			wiggle.read_depth.assign(wiggle.length, 0.0);
		}
		
		add_bam_record_to_wiggle(b, wiggle);
		++cnt;
		if (cnt % 1000000 == 0) std::cout<< cnt<< std::endl;
	}
	if (cur_tid >= 0) { used[cur_tid] = true; processor.process(wiggle); }

	for (int32_t i = 0; i < header->n_targets; i++)
		if (!used[i]) {
			wiggle.name = header->target_name[i];
			wiggle.length = header->target_len[i];
			wiggle.read_depth.clear();
			processor.process(wiggle);
		}

	samclose(bam_in);
	bam_destroy1(b);
	delete[] used;
}

UCSCWiggleTrackWriter::UCSCWiggleTrackWriter(const std::string& output_filename,
                                             const std::string& track_name) {
    fo = fopen(output_filename.c_str(), "w");
    fprintf(fo, "track type=wiggle_0 name=\"%s\" description=\"%s\" visibility=full\n",
            track_name.c_str(),
            track_name.c_str());
}

UCSCWiggleTrackWriter::~UCSCWiggleTrackWriter() {
    fclose(fo);
}

void UCSCWiggleTrackWriter::process(const Wiggle& wiggle) {
    int sp, ep;

    if (wiggle.read_depth.empty()) return;
    
    sp = ep = -1;
    for (size_t i = 0; i < wiggle.length; i++) {
        if (wiggle.read_depth[i] > 0) {
            ep = i;
        }
        else {
            if (sp < ep) {
                ++sp;
                fprintf(fo, "fixedStep chrom=%s start=%d step=1\n", wiggle.name.c_str(), sp + 1);
                for (int j = sp; j <= ep; j++) fprintf(fo, "%.7g\n", wiggle.read_depth[j]);
            }
            sp = i;
        }
    }
    if (sp < ep) {
        ++sp;
        fprintf(fo, "fixedStep chrom=%s start=%d step=1\n", wiggle.name.c_str(), sp + 1);
        for (int j = sp; j <= ep; j++) fprintf(fo, "%.7g\n", wiggle.read_depth[j]);
    }
}
