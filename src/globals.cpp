#include <blond/globals.h>
using namespace blond;
Slices *Context::Slice = nullptr;
GeneralParameters *Context::GP = nullptr;
Beams *Context::Beam = nullptr;
RfParameters *Context::RfP  = nullptr;

// TODO num of threads is global
// should it?
int Context::n_threads = 1;