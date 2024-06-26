// wiggle.h for CSEM

#ifndef WIGGLE_H_
#define WIGGLE_H_

#include <cstdio>
#include <string>
#include <vector>

extern bool no_fractional_weight; // if no_frac_weight == true, each alignment counts as weight 1
extern int fragment_length; // if -1, do not extend reads
extern bool only_midpoint; // if true, represent each fragment by its midpoint

struct Wiggle {
    std::string name;
    std::vector<float> read_depth;
    size_t length;
};

class WiggleProcessor {
public:
    virtual ~WiggleProcessor() {}
    virtual void process(const Wiggle& wiggle) = 0;
};

class UCSCWiggleTrackWriter : public WiggleProcessor {
public:
    UCSCWiggleTrackWriter(const std::string& output_filename,
                          const std::string& track_name);
        
    ~UCSCWiggleTrackWriter();

    void process(const Wiggle& wiggle);

private:
    FILE *fo;
};

void build_wiggles(const std::string& bam_filename,
                   WiggleProcessor& processor);

#endif
