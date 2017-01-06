#include <vector>

// "push_back" on lists is slow, so use deques and vectors 
//#define SET_TYPE std::list
#define SET_TYPE std::vector
#define PLIST_TYPE std::deque
//#define PLIST_TYPE std::list

struct gasP_t {
    double rho;
    double r[3];
    double v[3];
    double uint;
    double h;
    double m;
};

struct gizData_t {
    std::vector<gasP_t> gasP;
    int ng;
};
